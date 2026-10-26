#!/usr/bin/env python3
"""
virphy_bridge.py — 桥接 VirPhyKit 计算核心，提供无 GUI 的 CLI 可调用函数
=============================================================================
直接调用 VirPhyKit 源码中的纯计算逻辑，剥离 PyQt5 依赖。
每个函数对应一个 VirPhyKit 模块，返回结构化结果字典。

VirPhyKit 源码路径: <repo>/biosoft/VirPhyKit/src/（随仓库定位，不写死本机路径）
"""

import csv
import io
import json
import os
import re
import sys
import subprocess
import tempfile
import threading
import logging
import shutil
from pathlib import Path
from typing import Optional, Dict, List, Callable

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import linregress

# (2026-08-25 通用化: 不再依赖 VirPhyKit 源码, MJRM/Rename/TempMig 已自研)

from Bio import AlignIO, Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# 2026-10-12: 统一 SCI 配色/样式（本模块产出 Fig 3 的 I / L / M / O / Q / R 面板）。
# 双路径 import: 模块方式 `import utils.virphy_bridge` 与脚本方式均可。
try:
    from utils.sci_palette import (
        apply_sci_style, OKABE_ITO, OKABE_BRANCH, SEQ_CMAP,
        COL_W_1, COL_W_2, despine)
    # 2026-10-16: 图件双格式 (PNG+PDF 同主名); 统一由 plot_io.save_both 收口
    from utils.plot_io import save_both
except ImportError:                                    # pragma: no cover
    from sci_palette import (
        apply_sci_style, OKABE_ITO, OKABE_BRANCH, SEQ_CMAP,
        COL_W_1, COL_W_2, despine)
    from plot_io import save_both

# 2026-10-16: 第三方工具 (TreeTime/BEAST/treeannotator 等) 直写 PDF 的 PNG 兜底
try:
    from utils.pdf_to_png import stage_finalize as _pdf_png_finalize
except ImportError:                                    # pragma: no cover
    from pdf_to_png import stage_finalize as _pdf_png_finalize

# TreeTime
try:
    from treetime import TreeTime
    from treetime.utils import parse_dates
    HAS_TREETIME = True
except ImportError:
    HAS_TREETIME = False


# ═══════════════════════════════════════════════════════════════════
# 日志回调：替代 PyQt5 pyqtSignal
# ═══════════════════════════════════════════════════════════════════

class LogCollector:
    """收集日志消息，替代 Qt signals"""
    def __init__(self, logger: Optional[logging.Logger] = None):
        self.logger = logger or logging.getLogger(__name__)
        self.messages: List[str] = []

    def emit(self, msg: str):
        self.messages.append(msg)
        self.logger.info(msg)

    def warning(self, msg: str):
        self.messages.append(msg)  # 历史坑: warning/error 不入 messages,
        self.logger.warning(msg)   # 下游 log_messages 丢失告警。

    def error(self, msg: str):
        self.messages.append(msg)
        self.logger.error(msg)


# ═══════════════════════════════════════════════════════════════════
# BSSVS 迁移路线判据 —— 全管线唯一定义 (2026-09-15, 审查 P2-9)
# ═══════════════════════════════════════════════════════════════════
# 历史坑: 三个模块各自写默认值 —— beast_parser=5.0, spread3_viz=5.0,
# phylogeo_bridge=3.0 —— 同一份 log 换个入口调用就会得到不同的"显著路线"集合。
# CLI (--phylogeo_bf) 默认是 5.0, 故统一到 5.0; 需要别的阈值一律显式传参。
#
# BF 本体: BF = (p/0.5)/((1-p)/0.5), p = indicators 后验均值, 先验边际 0.5
# 对应 XML 里 poissonPrior mean="0.6931471805599453" (=ln2), 与 Lemey 2009 自洽。
DEFAULT_BF_THRESHOLD = 5.0

# p -> 1 时 BF 解析上发散。JSON 没有 Infinity 字面量 (json.dump 默认写出的
# Infinity 是非标准扩展, 严格解析器会拒绝), 故统一截断到一个有限哨兵值。
BF_INFINITE_CAP = 999.0


def find_indicator_columns(columns, require_location: bool = True):
    """从 log 列名里挑出 BSSVS 指标列 (全管线唯一定义, 2026-09-15 审查 P2-27 补)。

    历史坑: 三个模块各写一份检测条件 ——
        beast_parser     `'indicator' in c and 'location' in c`
        phylogeo_bridge  `'indicator' in c or  'indicators' in c`  (后者冗余)
        spread3_viz      `'indicator' in c`
    同一份 log 换个入口可能选中不同的列集合。统一为: 优先 location 作用域内的
    indicator 列; 若一个都没有 (老 XML 可能不按 location 命名), 再放宽到全部
    indicator 列, 避免造成回归。
    """
    ind = [c for c in columns if 'indicator' in str(c).lower()]
    if not require_location:
        return ind
    scoped = [c for c in ind if 'location' in str(c).lower()]
    return scoped if scoped else ind


def bayes_factor(posterior: float,
                 prior_expectation: float = 0.5,
                 cap: float = BF_INFINITE_CAP) -> float:
    """BSSVS 指标后验 → Bayes Factor (全管线唯一定义)。

    BF = (p / p0) / ((1 - p) / (1 - p0)),  p0 = 0.5 (poissonPrior mean=ln2)。

    · p 落在 (0, 1) 开区间 → 正常公式
    · p >= 1.0 (全部样本取 1) → 公式发散, 返回有限哨兵 cap (默认 999.0)
    · p <= 0.0 → 返回 0.0
    """
    if 0.0 < posterior < 1.0:
        return (posterior / prior_expectation) / \
               ((1.0 - posterior) / (1.0 - prior_expectation))
    return float(cap) if posterior >= 1.0 else 0.0


# ═══════════════════════════════════════════════════════════════════
# 统一 p-distance 口径
# ═══════════════════════════════════════════════════════════════════

_P_DISTANCE = None


def _fallback_p_distance(s1: str, s2: str) -> float:
    """distance_matrix 不可用时的等价兜底 (逻辑与其 p_distance 保持一致):
    跳过任一方为 gap 的位点; N 等 IUPAC 简并字符按展开集判等。"""
    _IUPAC = {'A': 'A', 'T': 'T', 'C': 'C', 'G': 'G',
              'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
              'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
    comp = diff = 0
    for a, b in zip(s1[:min(len(s1), len(s2))], s2[:min(len(s1), len(s2))]):
        a, b = a.upper(), b.upper()
        if a == '-' or b == '-':
            continue
        ea, eb = _IUPAC.get(a), _IUPAC.get(b)
        if not ea or not eb:
            continue
        comp += 1
        if not any(x == y for x in ea for y in eb):
            diff += 1
    return diff / comp if comp > 0 else 0.0


def p_distance_pair(s1: str, s2: str) -> float:
    """两序列 p-distance (统一口径)。

    2026-09-15 修复: 本模块原有三处手写 `diffs = sum(a != b ...)`, 在 MSA 上把
    gap 与 N 当真实差异 → within/between 距离与 Mantel/置换 p 值系统性偏高;
    且与 utils/distance_matrix.p_distance (跳过 gap、N 按 IUPAC 简并展开) 口径
    不一致 —— 同一份数据两个模块给出不同 p-distance。
    现统一复用 distance_matrix.p_distance, 保证全管线单一口径。
    """
    global _P_DISTANCE
    if _P_DISTANCE is None:
        try:
            try:
                from distance_matrix import p_distance as _pd
            except ImportError:
                from utils.distance_matrix import p_distance as _pd
            _P_DISTANCE = _pd
        except ImportError:
            _P_DISTANCE = _fallback_p_distance
    return float(_P_DISTANCE(s1, s2))


# ═══════════════════════════════════════════════════════════════════
# Module 1: TreeTime-RTT  (Root-to-tip Regression)
# ═══════════════════════════════════════════════════════════════════

# 捕获 TreeTime stdout 用的全局锁 (见 _run_treetime_capturing 注释)
_STDOUT_CAPTURE_LOCK = threading.Lock()


def _run_treetime_capturing(tt, **run_kwargs) -> str:
    """跑 tt.run(**run_kwargs), 尽量捕获其 stdout 供日志对比。

    2026-08-27 修复 (第三个 bug: 非线程安全, 是 DRT 并行静默失败的真正原因)
    -------------------------------------------------------------------
    旧版写法:
        old_stdout = sys.stdout
        sys.stdout = io.StringIO()
        try:
            tt.run(...)
            treetime_output = sys.stdout.getvalue()
        finally:
            sys.stdout = old_stdout

    这是直接替换**全局** sys.stdout。DRT fast 模式用 ThreadPoolExecutor 并跑
    多个 run_treetime_rtt, 两个线程会同时替换同一个全局变量:

        线程 A: old = real;  sys.stdout = bufA
        线程 B: old = bufA;  sys.stdout = bufB
        线程 A: bufB.getvalue()  -> 能跑, 但拿到的是 B 的缓冲
        线程 A: finally -> sys.stdout = real
        线程 B: sys.stdout.getvalue() -> 此时 sys.stdout 是 real
                -> AttributeError: '_io.TextIOWrapper' object has no attribute 'getvalue'
        线程 B: finally -> sys.stdout = bufA   <- 全局状态被永久弄坏

    最后一步是致命的: sys.stdout 被还原成一个已被丢弃的 StringIO, 之后本进程
    内**所有 print 静默消失**。实测 DRT round1 的日志整体丢失, 脚本看起来像
    "静默退出", 实际是输出全进了废弃缓冲。这直接导致 T4 验收连续两次失败。

    现方案: 模块级锁 + 非阻塞获取。
      - 拿到锁 (串行路径): 行为与旧版一致, 保留 stdout 交叉校验
      - 拿不到锁 (并行路径): 直接跑, 不碰全局 sys.stdout
    这样既不污染全局状态, 也不会把 TreeTime 串行化 (并行加速得以保留)。
    代价仅是并行时少一份 TreeTime 原始日志, 而 R²/β 早已改为读对象属性
    (见 P0 修复), 该日志仅用于交叉校验。
    """
    if not _STDOUT_CAPTURE_LOCK.acquire(blocking=False):
        # 有其它线程正在捕获, 本轮放弃捕获 (不碰全局 stdout)
        tt.run(**run_kwargs)
        return ""
    try:
        old_stdout = sys.stdout
        buf = io.StringIO()
        sys.stdout = buf
        try:
            tt.run(**run_kwargs)
            return buf.getvalue()
        finally:
            sys.stdout = old_stdout
    finally:
        _STDOUT_CAPTURE_LOCK.release()


def run_treetime_rtt(
    fasta_file: str,
    tree_file: str,
    dates_file: str,
    output_dir: str,
    mapping_file: Optional[str] = None,
    rng_seed: Optional[int] = 42,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    运行 TreeTime root-to-tip 回归分析。

    Parameters
    ----------
    fasta_file : str
        FASTA 多序列比对文件
    tree_file : str
        Newick 格式系统发育树
    dates_file : str
        CSV 格式日期文件，列: name, date
    output_dir : str
        输出目录
    mapping_file : str, optional
        性状映射文件 (region<TAB>location)，用于 color-coded 散点图
    rng_seed : int or None, default 42
        TreeTime 内部 RNG 的固定种子。必须给值才能复现结果。
        2026-08-27 修复: TreeAnc.__init__ 里 self.rng = np.random.default_rng(seed=rng_seed),
        rng_seed 默认 None 意味着每次从 OS 熵取新种子; 而 treetime.py:214 的 run()
        把 sample_from_profile 硬编码成 'root', 根序列永远从 profile 抽样
        (seq_utils.prof2seq -> rng.random), 无法经 kwargs 关闭。实测同一输入
        三次跑出 R² 从 0.0001 到 0.0667。故这里在构造后覆写 tt.rng 以复现。
        传 None 会退回不确定行为, 仅用于显式做种子敏感性分析。
    log : LogCollector, optional

    Returns
    -------
    dict with keys: success, beta, r_squared, p_value, r2_source, rng_seed,
        clock_rate, output_tree, plots, log_messages
    """
    if not HAS_TREETIME:
        return {"success": False, "error": "treetime not installed. Run: pip install phylo-treetime"}

    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "beta": None, "r_squared": None, "p_value": None,
              "r2_source": None, "rng_seed": rng_seed, "output_tree": None,
              "plots": {}, "log_messages": []}

    try:
        # ── 1. Load alignment ──
        alignment = AlignIO.read(fasta_file, "fasta")
        for record in alignment:
            record.seq = Seq(str(record.seq).replace('U', 'T'))
        log.emit(f"Loaded {len(alignment)} sequences from {fasta_file}")

        # ── 2. Load tree ──
        tree = Phylo.read(tree_file, "newick")
        log.emit(f"Loaded newick tree from {tree_file}")

        # ── 3. Load dates ──
        dates = parse_dates(dates_file)
        metadata = pd.read_csv(dates_file)
        has_location = 'location' in metadata.columns
        log.emit(f"Loaded {len(dates)} dates from {dates_file}")

        # ── 4. Handle tips without dates ──
        tree_tips = {clade.name for clade in tree.get_terminals()}
        missing = tree_tips - set(dates.keys())
        # Set outgroup dates to None so TreeTime infers them from tree structure
        if missing:
            log.warning(f"Tree tips without dates (will be inferred by TreeTime): {len(missing)} tips")
            for tip in missing:
                dates[tip] = None  # TreeTime will infer
        valid_tips = [t for t in tree_tips if dates.get(t) is not None]
        if len(valid_tips) < 3:
            result["error"] = f"Only {len(valid_tips)} tips have actual dates"
            result["log_messages"] = log.messages
            return result
        log.emit(f"{len(valid_tips)} tips with dates, {len(missing)} without (to be inferred)")

        # ── 5. Run TreeTime ──
        tt = TreeTime(
            aln=alignment, tree=tree, dates=dates,
            gtr='Jukes-Cantor', seq_len=len(alignment[0]), verbose=4
        )

        # ── 5b. 固定 RNG 种子 (2026-08-27 修复: 结果不确定性) ──
        # 根因: TreeAnc.__init__ 里 self.rng = np.random.default_rng(seed=rng_seed),
        #       rng_seed 默认 None -> 每次从 OS 熵取新种子。
        #       treetime.py:214 的 run() 把 sample_from_profile 硬编码为 'root'，
        #       根序列永远从 profile 抽样 (seq_utils.py:268 rng.random) 且无法经
        #       kwargs 关闭，所以只能在构造后覆写 tt.rng。
        #       实测: 同一输入连跑三次 R² = 0.0001 / 0.0527 / 0.0528 (极差 0.0527);
        #       固定种子后极差 0.0000000000。
        import numpy as _np
        if rng_seed is not None:
            tt.rng = _np.random.default_rng(int(rng_seed))
            log.emit(f"TreeTime RNG seed fixed to {rng_seed} (reproducible)")
        else:
            log.warning("TreeTime RNG seed is None: results are NOT reproducible")

        # 与 VirPhyKit 原版一致: 先最小二乘定根再 run(root='best')
        # (2026-08-27 修复: 之前直接 root='oldest', 与上游语义不一致)
        tt.reroot(root='least-squares')
        log.emit("Tree rerooted using least-squares method (as VirPhyKit)")

        # 2026-10-12 修复(真 bug): 在 run() **之前**存一份 reroot 后的树（枝长 = 替换数）。
        # 原因见下方第 8 步的说明：tt.run() 之后 tt.tree 会被换成 timetree。
        # 注意用模块级 `io`，**不能**用第 8 步才 `import io as _io` 的别名 —— 那会
        # NameError 并被下面的 except 静默吞掉，退回 timetree（实测踩过这个坑）。
        _rtt_nwk = None
        try:
            _b = io.StringIO()
            Phylo.write(tt.tree, _b, "newick")
            _b.seek(0)
            _rtt_nwk = _b.read()
        except Exception as _e:
            log.warning(f"无法保存 run() 前的替换数树, RTT 将退回 timetree: {_e}")

        # Capture stdout for β/R² extraction (log-only; 结果改读 TreeTime 对象属性)
        # 2026-08-27 修复: 旧版直接替换全局 sys.stdout, 并行时非线程安全,
        # 会抛 AttributeError 并把 sys.stdout 永久还原成废弃的 StringIO。
        # 详见 _run_treetime_capturing 的注释。
        #
        # raise_uncaught_exceptions=True 是必需的 (2026-08-27 修复第二个静默故障):
        # treetime.py:55-81 的 TreeTime.run() 把内部所有异常都吞掉后 sys.exit(2):
        #     def run(self, raise_uncaught_exceptions=False, **kwargs):
        #         try:
        #             return self._run(**kwargs)
        #         except TreeTimeError as err:
        #             ...
        #             sys.exit(2)
        #         except BaseException as err:
        #             ...
        #             sys.exit(2)
        # sys.exit 抛的是 SystemExit(BaseException), 不会被本函数的 except Exception
        # 拦住, 会直接穿透并终止**整个管线进程**。批处理下一个坏基因能干掉整轮;
        # 更隐蔽的是在 ThreadPoolExecutor worker 里, SystemExit 被
        # future.set_exception 存下, 无人 result() 时就彻底静默。
        # 传 True 后 TreeTimeError 原样抛出、其余包成 TreeTimeUnknownError
        # (两者都是 Exception 子类, from err 保留原因), 回到本函数的错误处理。
        treetime_output = _run_treetime_capturing(
            tt, root="best", branch_length_mode='joint',
            raise_uncaught_exceptions=True,
            infer_gtr=False, infer_clock=True,
            resolve_polytomies=True, time_marginal=True,
            max_iter=3, reconstruct_tip_states=True)

        log.emit("Time-tree inference completed")

        # ── 6. Extract β and R² ──
        #
        # 2026-08-27 修复 (P0): 旧版照抄 VirPhyKit function_treetime.py, 用 re.search 解析
        # Treetime stdout 取 R²/β。re.search 返回首个匹配, 而 tt.run() 内部有多次
        # get_clock_model 迭代, 每次都会打印一行 "R^2=", 所以拿到的是迭代中间值。
        #
        # 服务器实测 (G 基因, /tmp/_verify_r2_context.py):
        #   1: rate=2.610e-04  R^2=0.0351   ← 旧版抓到这个
        #   7: rate=3.272e-04  R^2=0.0526   ← 最终值 (tt.date2dist.r_val**2)
        # R² 偏低 33%, β 偏低 20%。
        #
        # 另外首支正则 r"R\^2 of root-to-tip regression:" 实测从未命中 (返回 None),
        # 实际一直靠 fallback r"R\^2\s*=". 该正则不是任何 treetime 版本的输出格式。
        #
        # 改为直接读 TreeTime 对象属性 (与 treedater/shinyTempSignal 的做法一致,
        # 即读统计对象而非解析日志)。stdout 解析降级为仅日志用途。
        beta = None
        r_squared = None
        r2_api_source = None

        date2dist = getattr(tt, 'date2dist', None)
        if date2dist is not None:
            # r_val 是 Pearson r (treeregression.py:~317: clock_model['r_val'] = explained_variance())
            # 主流程 covariation=False (ClockTree.__init__ 默认 use_covariation=False),
            # 此时加权退化, r_val**2 即普通决定系数, 与 OLS 的 r.squared 可通约。
            r_val = getattr(date2dist, 'r_val', None)
            if r_val is not None:
                r_squared = float(r_val) ** 2
                r2_api_source = 'date2dist.r_val**2'
            clock_rate = getattr(date2dist, 'clock_rate', None)
            if clock_rate is not None:
                beta = float(clock_rate)

        # stdout 解析: 只用于日志, 不再写入 result
        log_r2 = None
        log_beta = None
        r2_match = re.search(r"R\^2 of root-to-tip regression:\s*([\d\.]+)", treetime_output)
        if not r2_match:
            r2_match = re.search(r"R\^2\s*=\s*([\d\.]+)", treetime_output)
        if r2_match:
            log_r2 = float(r2_match.group(1))
        beta_match = re.search(r"ClockTree\.date2dist: Setting new molecular clock\. rate\s*([\d\.eE+-]+)", treetime_output)
        if not beta_match:
            beta_match = re.search(r"rate\s*=\s*([\d\.eE+-]+)", treetime_output)
        if beta_match:
            log_beta = float(beta_match.group(1))

        if r_squared is None:
            # 对象属性不可用时的兑底 (正常路径不该走到这里)
            r_squared = log_r2
            r2_api_source = 'stdout_fallback' if log_r2 is not None else None
            log.warning(f"R² taken from stdout fallback (date2dist.r_val unavailable): {r_squared}")
        if beta is None:
            beta = log_beta
            log.warning(f"β taken from stdout fallback (date2dist.clock_rate unavailable): {beta}")

        if log_r2 is not None and r_squared is not None and abs(log_r2 - r_squared) > 1e-9:
            log.emit(f"  [note] stdout R² ({log_r2}) 与最终 R² ({r_squared:.6f}) 不同, "
                     f"前者为迭代中间值, 已采用最终值")

        result["r2_source"] = r2_api_source
        clock_rate = getattr(date2dist, 'clock_rate', None)
        log.emit(f"Clock rate: {clock_rate}, β: {beta}, R²: {r_squared} "
                 f"(source: {r2_api_source})")

        # ── 7. Save time-tree ──
        output_tree_file = os.path.join(output_dir, "timetree_inferred.nwk")
        Phylo.write(tt.tree, output_tree_file, "newick")
        result["output_tree"] = output_tree_file

        result["beta"] = beta
        result["r_squared"] = r_squared
        # 2026-09-15 (审查 P2-8): clock_rate 原先只写日志不落 result, 导致
        # temporal_signal.date_randomization_test 承诺的 real_clock_rate / rate_overlap
        # 永远是 None。这里把它作为正式字段返回。
        result["clock_rate"] = float(clock_rate) if clock_rate is not None else None

        # ── 8. 根到端回归 p 值 (核心指标) ──────────────────────────────────
        # 2026-09-15 (审查 P2-12) 修复: 旧版把 linregress + p_value 的计算整个塞在
        # 下面"绘图 best-effort"的 try 里。绘图一失败 (matplotlib/后端/递归深度),
        # except 就把异常吞掉, p_value 静默变回 None —— 而日志还写着
        # "metrics preserved"。现在拆成独立块: 只包住计算, 绘图失败不影响它。
        # NOTE: TreeTime modifies Biopython tree node internals (adds __dict__ attrs)
        # which breaks Biopython's recursive traversal (find_any, get_path, depths).
        # Workaround: Newick round-trip → clean Biopython tree.
        # If you hit RecursionError from find_any() on a TreeTime tree,
        # use: clean_tree = Phylo.read(StringIO(newick_str), 'newick')
        import sys as _sys
        import io as _io
        _old_limit = _sys.getrecursionlimit()
        clean_tree = None
        root_to_tip_dists = []
        sampling_dates_list = []
        reg = None
        try:
            _sys.setrecursionlimit(max(_old_limit, 5000))
            # timetree 快照（第 10 步画 TimeTree 用它）
            nwk_buf = _io.StringIO()
            Phylo.write(tt.tree, nwk_buf, "newick")
            nwk_buf.seek(0)
            clean_tree = Phylo.read(nwk_buf, "newick")

            # 2026-10-12 修复(真 bug, P0): RTT 回归**必须**用替换数枝长，
            # 不能用 tt.tree。原因：tt.run() 返回后 TreeTime 会把 tt.tree
            # 换成 timetree（枝长 = 时间/年），于是 y 轴变成"年到根"——
            # 与采样日期必然近乎完美线性。GCVA 实测：
            #   · 用 tt.tree（timetree）: OLS R² = 0.9998, slope = 1.59e-4, p = 6.7e-193
            #     ← 这是伪相关（图上线与点完全重合），而日志/报告同时写着
            #       TreeTime 的 date2dist.r_val² = 0.0518，读者看到的图与数字自相矛盾。
            #   · 用 run() 前的替换数树: OLS R² = 0.0195, slope = 1.04e-4, p = 0.156
            #     ← 与 r_val² = 0.0518 同量级，才是真实的分子钟信号。
            rtt_tree = (Phylo.read(_io.StringIO(_rtt_nwk), "newick")
                        if _rtt_nwk else clean_tree)

            tip_names = [clade.name for clade in rtt_tree.get_terminals()]
            # Filter to tips that exist in dates (exclude synthetic nodes from rerooting)
            valid_tips = [n for n in tip_names if n in dates and dates[n] is not None]
            for tip_name in valid_tips:
                tip = rtt_tree.find_any(tip_name)
                path = rtt_tree.get_path(tip)
                dist = sum(clade.branch_length for clade in path
                           if clade.branch_length is not None)
                root_to_tip_dists.append(dist)
                sampling_dates_list.append(dates[tip_name])

            if len(root_to_tip_dists) >= 3:
                reg = linregress(sampling_dates_list, root_to_tip_dists)
                # 始终存 float (历史坑: <1e-10 时存字符串, 下游 float(p_value) 崩)
                result["p_value"] = float(reg.pvalue)
                result["rtt_n_tips"] = len(root_to_tip_dists)
                result["rtt_r_squared"] = float(reg.rvalue) ** 2
                log.emit(f"Root-to-tip regression: p={result['p_value']:.4g} "
                         f"(n={result['rtt_n_tips']} tips, "
                         f"OLS R²={result['rtt_r_squared']:.4f})")
            else:
                log.warning(f"Root-to-tip 回归跳过: 有效 tip 数 {len(root_to_tip_dists)} < 3")
        except Exception as _rtt_err:
            log.warning(f"Root-to-tip 回归计算失败 (p_value 保持 None): {_rtt_err}")
        finally:
            _sys.setrecursionlimit(_old_limit)

        # ── 9-10: 绘图 (best-effort; 核心指标此时已捕获) ──
        try:
            if clean_tree is not None and reg is not None:
                # 2026-10-12: 配色统一 OKABE（原 #00A0E9 亮蓝 / #E9522E 橙红均不在
                # 项目调色板内）；dpi 150 → 300；补 despine + y 轴网格。
                # 画布 10 in 保留 —— 排进 Fig 3 的 2 栏（4.58 in）缩放 0.46 倍，
                # 18 pt 标题正好落到 8.2 pt，字号比例本身是合适的。
                apply_sci_style(base_size=11)
                fig, ax = plt.subplots(figsize=(10, 6))
                ax.scatter(sampling_dates_list, root_to_tip_dists,
                           color=OKABE_ITO[4], s=14, alpha=0.75,
                           edgecolors="none", zorder=3)
                x_range = [min(sampling_dates_list), max(sampling_dates_list)]
                y_pred = [reg.slope * x + reg.intercept for x in x_range]
                ax.plot(x_range, y_pred, '--', color=OKABE_ITO[5], linewidth=1.8,
                        zorder=2, label="least-squares fit")
                text_lines = []
                if beta is not None:
                    text_lines.append(f"Clock rate: {beta:.3e} subs/site/year")
                # 2026-10-12: 同时给出两个 R² 并注明来源 —— TreeTime 的
                # date2dist.r_val²（管线口径，写进 metrics）与本图这条 OLS 线的
                # R² 在数值上不同（GCVA 实测 0.052 vs 0.020），只显示一个会让
                # 读者以为图与数字矛盾。P 值改为 3 位有效数字（原样输出会是
                # 6.689802503754485e-193 这种 20 位小数）。
                if r_squared is not None:
                    text_lines.append(f"TreeTime $R^2$ = {r_squared:.3f}")
                if reg is not None:
                    text_lines.append(f"OLS fit $R^2$ = {reg.rvalue ** 2:.3f}")
                if result.get("p_value") is not None:
                    text_lines.append(f"$P$ = {result['p_value']:.3g}")
                ax.text(0.02, 0.98, '\n'.join(text_lines), transform=ax.transAxes,
                        fontsize=12, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='white',
                                  edgecolor="#BBBBBB", alpha=0.9))
                ax.set_title("Root-to-tip regression", fontsize=18)
                ax.set_xlabel("Sampling date", fontsize=16)
                ax.set_ylabel("Distance to root (subs/site)", fontsize=16)
                ax.tick_params(labelsize=14)
                ax.grid(alpha=0.35, axis="y")
                despine(ax)
                ax.legend(fontsize=11, loc="lower right")
                rtt_plot = os.path.join(output_dir, "RootToTip_regression.pdf")
                save_both(fig, rtt_plot, bbox_inches='tight', dpi=300, format='pdf')
                plt.close(fig)
                result["plots"]["rtt"] = rtt_plot

            # Time-tree plot
            if clean_tree is not None:
                try:
                    fig2, ax2 = plt.subplots(figsize=(12, 8))
                    Phylo.draw(clean_tree, axes=ax2, do_show=False,
                               label_func=lambda c: c.name if c.name and not c.name.startswith('NODE_') else "",
                               branch_labels=lambda c: "")
                    ax2.set_title("Time-calibrated phylogeny", fontsize=16)
                    tree_plot = os.path.join(output_dir, "TimeTree.pdf")
                    save_both(fig2, tree_plot, bbox_inches='tight', dpi=300, format='pdf')
                    plt.close(fig2)
                    result["plots"]["tree"] = tree_plot
                except RecursionError:
                    log.warning("Time-tree plot skipped (tree too deep for Phylo.draw)")
                    # Write Newick as fallback deliverable
                    tree_nwk = os.path.join(output_dir, "TimeTree.nwk")
                    Phylo.write(tt.tree, tree_nwk, "newick")
                    result["plots"]["tree_nwk"] = tree_nwk
            log.emit(f"Plots saved: {list(result['plots'].keys())}")
        except Exception as plot_err:
            # 2026-09-15 (P2-12): 这里现在真的只影响图; p_value 已在上面的
            # 独立块里算完并写进 result, 不会再被绘图异常吞掉。
            log.warning(f"Plot generation failed (core metrics already captured): {plot_err}")

        result["success"] = True
        result["log_messages"] = log.messages
        # 2026-10-16: 收尾兜底 —— TreeTime 等第三方程序直写的 PDF 若没有同名 PNG,
        # 这里统一补一张 (幂等; 失败不抛错, 不影响 RTT 结果)。
        try:
            _made = _pdf_png_finalize(output_dir, log=log)
            if _made:
                log.emit(f"TreeTime-RTT: PDF→PNG 兜底渲染 {len(_made)} 张")
        except Exception as e:  # noqa: BLE001
            log.warning(f"PDF→PNG 兜底跳过 ({e})")
        return result

    except KeyboardInterrupt:
        raise
    except BaseException as e:
        # 2026-08-27 修复: 这里原本是 except Exception。
        # treetime.py:55-81 的 TreeTime.run() 内部把异常吞掉后 sys.exit(2),
        # 抛出的 SystemExit 是 BaseException 而非 Exception, 不 catch 会直接
        # 穿透并终止整个管线进程 (实测症状: 脚本"静默死亡", TreeTime 的
        # traceback 走 stderr 被吞, 什么线索都不剩)。
        # 虽然已给 tt.run() 传 raise_uncaught_exceptions=True 让它改抛普通异常,
        # 这里再兜一层, 防其他代码路径抛出 SystemExit。
        # 另: raise_uncaught_exceptions=True 会把未知异常包成无 message 的
        # TreeTimeUnknownError, 真正原因在 __cause__, 故展开 cause 链。
        cause = e.__cause__
        detail = f"{type(e).__name__}: {e}"
        _seen = 0
        while cause is not None and _seen < 5:
            detail += f"  <- {type(cause).__name__}: {cause}"
            cause = cause.__cause__
            _seen += 1
        log.error(f"TreeTime-RTT failed: {detail}")
        result["error"] = detail
        result["log_messages"] = log.messages
        return result


# ═══════════════════════════════════════════════════════════════════
# Module 2: TreeDater-LTT (Lineage-Through-Time + Molecular Dating)
# ═══════════════════════════════════════════════════════════════════

TREEDATER_R_SCRIPT = """
suppressMessages({
    require(treedater, quietly=TRUE)
    require(ape, quietly=TRUE)
    require(ggplot2, quietly=TRUE)
    require(jsonlite, quietly=TRUE)
})

args <- commandArgs(trailingOnly=TRUE)
tree_file <- args[1]
metadata_file <- args[2]
seqlen <- as.numeric(args[3])
output_dir <- args[4]
plot_ltt <- as.logical(args[5])
ncpu <- max(1, as.integer(args[6]))
nboot <- max(20, as.integer(args[7]))

# ── 2026-08-31 增强: temporalConstraints 开启 (outlierTips 的前提) ──
dtr <- dater(tre <- read.tree(tree_file),
             sts <- setNames(read.csv(metadata_file, header=TRUE)[,2],
                             read.csv(metadata_file, header=TRUE)[,1]),
             seqlen, clock="uncorrelated", temporalConstraints=TRUE)

out <- list()
# meanRate 可能是长度 2 向量 (失败时 0,Inf)，取第一个有限值
mr <- dtr$meanRate
if (length(mr) > 1) mr <- mr[is.finite(mr) & mr > 0][1]
out$tmrca <- if (length(dtr$TimeOfMRCA)) as.numeric(dtr$TimeOfMRCA) else NULL
out$time_to_mrca <- if (length(dtr$TimeToMRCA)) as.numeric(dtr$TimeToMRCA) else NULL
out$mean_rate <- if (length(mr) == 1 && is.finite(mr)) as.numeric(mr) else NULL
if (!is.null(dtr$omegas) && length(dtr$omegas) > 1)
    out$cov_rate <- sd(dtr$omegas) / mean(dtr$omegas)  # CoV of lineage rates
out$clock <- dtr$clock

# 2026-10-12（审稿人视角修复）：旧调用 `plot(dtr, no.mar=TRUE, cex=0.6)` 里的
# `no.mar` **不是** ape::plot.phylo 的参数（ape 里叫 `no.margin`），所以它被 `...`
# 吞掉，实际等于 `plot.phylo(dtr, cex=0.6)`。plot.phylo 自留的上边距很小，
# 1.6 倍的标题被顶出画布 —— 实测标题上半截被裁掉。
# 现改为 `par(oma=...)` 预留外层上边距 + `mtext(outer=TRUE)` 写标题：
# 标题落在预留带内，且完全不动 plot.phylo 自己画的坐标轴。
## 2026-10-16: 绘图体抽成闭包, 供 PDF / PNG 两个设备复用 (图本身一字未改)
.draw_dated <- function() {
    par(oma=c(0, 0, 3.4, 0), family="sans")
    plot(dtr, cex=0.6)
    .ntip <- if (!is.null(dtr$tip.label)) length(dtr$tip.label) else NA
    mtext("Time-calibrated phylogeny (TreeDater)", side=3, line=1.9, cex=1.5,
          font=2, outer=TRUE)
    if (!is.na(.ntip))
        mtext(sprintf("n = %d tips", .ntip), side=3, line=0.55, cex=1.0,
              outer=TRUE)
}
pdf(file.path(output_dir, "Phylogeny_dated.pdf"), width=10, height=8)
.draw_dated()
dev.off()

## 2026-10-16: 图件双格式 (PNG+PDF 同主名) —— 同图再出 PNG (300 dpi, 白底)
png(file.path(output_dir, "Phylogeny_dated.png"), width=10, height=8,
    units="in", res=300, bg="white")
.draw_dated()
dev.off()

# ── 参数化自举: TMRCA/速率 CI (ncpu 并行) + LTT 图 ──
# 注: parboot 返回对象 CI 在 meanRate_CI / timeOfMRCA_CI 字段 (bootTreedater 结构，实测)
pb <- try(parboot(dtr, nreps=nboot, ncpu=ncpu), silent=TRUE)
if (!inherits(pb, 'try-error')) {
    if (!is.null(pb$meanRate_CI))
        out$rate_ci <- as.numeric(pb$meanRate_CI)
    if (!is.null(pb$timeOfMRCA_CI))
        out$tmrca_ci <- as.numeric(pb$timeOfMRCA_CI)
    if (!is.null(pb$coef_of_variation_CI))
        out$cov_ci <- as.numeric(pb$coef_of_variation_CI)
    if (plot_ltt) {
        # 2026-10-12（审稿人视角修复）：`plot.bootTreedater(ggplot=TRUE)` 把 HPD 带
        # **硬编码成 `fill="blue", col="blue"`**（R 的 "blue" = #0000FF，不在项目
        # 调色板内，渲染出来是淡紫），并且**没有图例**、也没说明阴影代表什么。
        # 现复用它算好的数据框，重建图层：OKABE 蓝 HPD 带 + 中位线 + ML 曲线，
        # 并补图例与量纲明确的轴标签。
        .g0 <- plot(pb, ggplot=TRUE)
        .ltt_df <- .g0$data
        g <- ggplot(.ltt_df, aes(x=times)) +
             geom_ribbon(aes(ymin=lb, ymax=ub, fill="95% CI (parboot)"),
                         alpha=0.30, colour=NA) +
             geom_line(aes(y=median, colour="Median (parboot)"), linewidth=0.9) +
             geom_line(aes(y=pml, colour="Maximum likelihood"),
                       linewidth=0.8, linetype="dashed") +
             scale_fill_manual(name=NULL,
                               values=c("95% CI (parboot)"="#56B4E9")) +
             scale_colour_manual(name=NULL,
                                 values=c("Median (parboot)"="#0072B2",
                                          "Maximum likelihood"="#D55E00")) +
             theme_bw(base_size=13) +
             theme(panel.grid.minor=element_blank(),
                   plot.title=element_text(face="bold", hjust=0.5),
                   legend.position=c(0.015, 0.985),
                   legend.justification=c(0, 1),
                   legend.background=element_rect(fill="white", colour=NA),
                   legend.key=element_rect(fill="white", colour=NA)) +
             labs(title="Lineages-through-time (TreeDater)",
                  x="Time (year)", y="Number of lineages")
        ggsave(file.path(output_dir, "LTT.pdf"), plot=g, width=10, height=8, dpi=300)
        ## 2026-10-16: 图件双格式 (PNG+PDF 同主名) —— 补同主名 PNG (300 dpi, 白底)
        ggsave(file.path(output_dir, "LTT.png"), plot=g, width=10, height=8,
               dpi=300, bg="white")
    }
}

# ── 严格/松弛钟检验 (parboot 似然比; 失败不阻断) ──
rct <- try(relaxedClockTest(dtr, nreps=nboot, ncpu=ncpu), silent=TRUE)
if (!inherits(rct, 'try-error')) {
    # relaxedClockTest 返回 htest 类对象
    if (!is.null(rct$p.value)) out$relaxed_clock_p <- as.numeric(rct$p.value)
    if (!is.null(rct$statistic)) out$relaxed_clock_lr <- as.numeric(rct$statistic)
}

# ── 谱系离群检测 (枝速率 FDR; 失败不阻断) ──
ot <- try(outlierTips(dtr, alpha=0.05), silent=TRUE)
if (!inherits(ot, 'try-error') && !is.null(ot)) {
    if (nrow(ot)) {
        write.table(ot, file.path(output_dir, "outlier_tips.tsv"),
                    sep="\t", quote=FALSE, row.names=FALSE)
        out$n_outlier_tips <- nrow(ot)
    } else out$n_outlier_tips <- 0
}

cat(toJSON(out, auto_unbox=TRUE, digits=10), sep="\n")
cat("\n___TD_JSON_END___\n")
cat(capture.output(print(dtr)), sep="\n")
"""


def run_treedater_ltt(
    tree_file: str,
    metadata_file: str,
    seq_len: int,
    output_dir: str,
    plot_ltt: bool = True,
    rscript_path: str = "Rscript",
    threads: int = 40,
    n_boot: int = 100,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    运行 TreeDater 分子定年 + Lineage-Through-Time 分析。

    2026-08-31 增强 (treedater 四高价值功能):
      temporalConstraints=TRUE 开启;
      parboot 数值 CI (TMRCA/速率, ncpu 并行);
      relaxedClockTest 严格/松弛钟检验;
      outlierTips 谱系离群检测 (outlier_tips.tsv);
      CoV of rates 进 metrics。

    Parameters
    ----------
    tree_file : str
        Newick 格式系统发育树
    metadata_file : str
        CSV 格式，列: name, date
    seq_len : int
        序列长度 (bp)
    output_dir : str
        输出目录
    plot_ltt : bool
        是否绘制 LTT 图
    rscript_path : str
        Rscript 可执行文件路径
    threads : int
        parboot/relaxedClockTest 并行核数
    n_boot : int
        自举次数 (parboot + relaxedClockTest)
    log : LogCollector, optional

    Returns
    -------
    dict with keys: success, output, plots, metrics, outlier_tips, error
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "output": "", "plots": {}, "error": None}

    try:
        # Convert dates to decimal years for TreeDater (it expects numeric)
        dates_df = pd.read_csv(metadata_file)
        name_col = dates_df.columns[0]
        date_col = dates_df.columns[1]

        # 2026-09-16 (上游核对 Task#10): 旧版此处内联 `_to_decimal`, 其口径有两处
        # 与全管线唯一实现 utils/decimal_year.to_decimal_year 不一致:
        #   ① YYYY-MM-DD 用 `(tm_yday-1)/365`  —— 固定 365 天年, 闰年与跨月有漂移;
        #      统一实现用 `(mo-1 + (dy-1)/days_in_month(mo))/12` —— 当月实际天数。
        #      仅此一项, 同一条日期最多相差约 1.57 天 (2021-03-15 实测)。
        #   ② YYYY / YYYY-MM 分别退化为"年初" / "月初", 而统一实现约定取 6-15 / 15 日。
        #      二者对"只给年份"的 metadata 会差约 166 天 (0.46 年) —— 量级远超 ①,
        #      这正是本函数**不可**直接换成统一实现的原因, 见下方 _to_decimal 包装。
        #
        # 结论: 本路径改为调用 utils/decimal_year (消除 ① 的 1.57 天漂移), 但年份/年月
        # 两级**显式保留上游 VirPhyKit 原口径** (只给年 → 年初; 只给年月 → 月初),
        # 因为 VirPhyKit 的 dendropy/treetime 侧同样是"只给年即年初"。若在此处改成
        # 年中/月中, Treedater 会与同一份数据在 BEAST2 路径上算出的采样时间不一致。
        # 这是**有意差异**, 不是遗漏; 详见 UPSTREAM_CONSISTENCY_20260916.md §6。
        #
        # 2026-09-16 (老师拍板 · 接入统一入口): 原先此处内联的 `_to_decimal` 已改为
        # **委托 utils/metadata_governance.decimal_year(mode='start')** —— 口径不变
        # (只给年→年初 / 只给年月→月初 / 三段→当月实际天数), 数值**逐位一致**
        # (verify_metadata_governance.py D2 实测)。收益:
        #   · 日期解析集中一处, 增补格式 (DD-Mon-YYYY 等) 全管线同步受益;
        #   · 占位符 (字面量 'YYYY-MM-DD'/'NA') 统一被拦, 不再静默造年份;
        #   · 口径不再藏在实现里, 由 mode 参数显式表达。
        from utils.metadata_governance import decimal_year as _gov_decimal_year

        def _to_decimal(d):
            """Convert YYYY, YYYY-MM, YYYY-MM-DD to decimal year.

            委托 utils/metadata_governance.decimal_year(mode='start'):
              YYYY-MM-DD → 当月实际天数 (utils/decimal_year 统一口径);
              YYYY-MM    → yr + (mo-1)/12   (VirPhyKit 原口径, 月初);
              YYYY       → yr               (VirPhyKit 原口径, 年初)。
            解析失败 / 占位符 → None (由下方 dropna 剔除, 不静默造年份)。
            """
            try:
                return _gov_decimal_year(d, mode='start')
            except (ValueError, IndexError, TypeError):
                # 历史坑: 非数字日期串 (用户自备 dates.csv) 直接 ValueError 崩整个函数
                return None

        dates_df['_decimal_year'] = dates_df[date_col].apply(_to_decimal)
        # Drop rows without dates
        dates_df = dates_df.dropna(subset=['_decimal_year'])

        if len(dates_df) < 3:
            result["error"] = f"Only {len(dates_df)} tips have valid dates (need ≥3)"
            log.warning(result["error"])
            return result

        dated_tips = set(dates_df[name_col].values)
        n_original = 0
        try:
            from Bio import Phylo
            tree = Phylo.read(tree_file, "newick")
            all_tips = [c for c in tree.get_terminals()]
            n_original = len(all_tips)
            # Prune tree to only dated tips
            undated = [c for c in all_tips if c.name not in dated_tips]
            if undated:
                for tip in undated:
                    tree.prune(tip)
                pruned_tree = os.path.join(output_dir, "timetree_pruned.nwk")
                Phylo.write(tree, pruned_tree, "newick")
                tree_file = pruned_tree
                n_pruned = len(list(tree.get_terminals()))
                log.emit(f"Pruned tree: {n_original} → {n_pruned} tips (only dated)")
        except Exception as e:
            log.warning(f"Tree pruning failed (continuing with original): {e}")

        # Write temp dates file
        dates_temp = os.path.join(output_dir, "dates_decimal.csv")
        dates_df[[name_col, '_decimal_year']].to_csv(dates_temp, index=False, header=[name_col, 'date'])
        log.emit(f"Converted {len(dates_df)} dates to decimal years")
        with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as f:
            f.write(TREEDATER_R_SCRIPT)
            r_script_path = f.name

        cmd = [
            rscript_path, r_script_path,
            tree_file, dates_temp,
            str(seq_len), output_dir,
            "TRUE" if plot_ltt else "FALSE",
            str(threads), str(n_boot)
        ]
        log.emit(f"Running: {' '.join(cmd)}")

        proc = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=1800)
        os.unlink(r_script_path)

        if proc.returncode != 0:
            log.error(f"TreeDater R script failed: {proc.stderr}")
            result["error"] = proc.stderr
            return result

        result["output"] = proc.stdout.strip()
        # ── 解析 JSON 指标块 (___TD_JSON_END___ 之前) ──
        result["metrics"] = {}
        if "___TD_JSON_END___" in proc.stdout:
            try:
                json_blob = proc.stdout.split("___TD_JSON_END___")[0].strip()
                # 取最后一个完整 JSON 行 (前面可能有 R 的启动信息)
                for line in reversed(json_blob.splitlines()):
                    line = line.strip()
                    if line.startswith("{") and line.endswith("}"):
                        result["metrics"] = json.loads(line)
                        break
            except Exception as _e:
                log.warning(f"TreeDater metrics JSON 解析失败 (非阻断): {_e}")
        # 2026-10-11 修 (P0): treedater 指标原先**只存在内存**, 从不落盘 ——
        # 于是任何「checkpoint 跳过 clock 后重跑 report」都会丢掉 td_mean_rate /
        # td_relaxed_clock_p / td_tmrca (headline_metrics 的**首选**来源),
        # 汇总表 Beta/R²/P/TMRCA 全空 (PSTVd / GCVA 实测双双中招)。
        # 现把 metrics 持久化到 <output_dir>/metrics.json, 供 skip 分支回填。
        if result["metrics"]:
            try:
                _mj = os.path.join(output_dir, "metrics.json")
                with open(_mj, "w") as _mf:
                    json.dump(result["metrics"], _mf, indent=2, default=str)
                result["metrics_json"] = _mj
            except OSError as _e:
                log.warning(f"TreeDater metrics 落盘失败 (非阻断): {_e}")
        _otsv = os.path.join(output_dir, "outlier_tips.tsv")
        if os.path.exists(_otsv):
            result["outlier_tips"] = _otsv
            n_ot = result["metrics"].get("n_outlier_tips")
            if n_ot is not None:
                log.emit(f"TreeDater outlier tips: {n_ot} (q<0.05) → outlier_tips.tsv")
        result["plots"]["phylogeny"] = os.path.join(output_dir, "Phylogeny_dated.pdf")
        if plot_ltt:
            ltt_path = os.path.join(output_dir, "LTT.pdf")
            if os.path.exists(ltt_path):
                result["plots"]["ltt"] = ltt_path

        result["success"] = True
        log.emit("TreeDater completed successfully")
        return result

    except subprocess.TimeoutExpired:
        result["error"] = "TreeDater 超时 (>10分钟)"
        log.error(result["error"])
        return result
    except Exception as e:
        err_str = str(e)
        if 'Missing sample time' in err_str:
            hint = ("TreeDater: 部分序列缺少采样时间。"
                    "请检查 dates.csv 是否包含所有 tip 的日期，"
                    "或使用 TreeTime 先推断缺失日期。")
        elif 'non-numeric' in err_str:
            hint = ("TreeDater: 日期格式必须为数值。"
                    "请使用 YYYY (如 2021) 或小数年 (如 2021.5) 格式。")
        else:
            hint = f"TreeDater R 脚本失败: {err_str[:200]}"
        result["error"] = hint
        log.error(hint)
        return result


# ═══════════════════════════════════════════════════════════════════
# Module 3: RSPP-Viz (Root State Posterior Probability)
# ═══════════════════════════════════════════════════════════════════

def run_rspp_viz(
    tree_file: str,
    output_dir: str,
    chart_type: str = "bar",
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    从 BEAST 注释树解析根状态后验概率并绘图。

    Parameters
    ----------
    tree_file : str
        BEAST 注释树 (.tre / .tree)，含 set.prob={} 注释
    output_dir : str
        输出目录
    chart_type : str
        "bar" 或 "pie"
    log : LogCollector, optional

    Returns
    -------
    dict with keys: success, root_states, plot, error
    """
    import re
    import random
    import seaborn as sns

    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "root_states": {}, "plot": None, "error": None}

    try:
        with open(tree_file, 'r', encoding='utf-8') as f:
            content = f.read()

        # Parse .set={} and .set.prob={} (VirPhyKit readTreeFile logic)
        # Handles both "Region.set={...}" and "set={...}"
        set_patterns = re.findall(r'(?:\w+\.)?set=\{([^\}]+)\}', content)
        prob_patterns = re.findall(r'(?:\w+\.)?set\.prob=\{([^\}]+)\}', content)
        if not set_patterns or not prob_patterns:
            result["error"] = "No set={} or set.prob={} annotations found"
            log.error(result["error"])
            return result

        # Take last occurrence (root node in BEAST tree)
        last_set = set_patterns[-1]
        last_prob = prob_patterns[-1]

        # Parallel lists: regions and their probabilities
        regions_list = [f.strip().strip('"') for f in last_set.split(',')]
        prob_list = [float(v.strip()) for v in last_prob.split(',')]

        states = {r: p for r, p in zip(regions_list, prob_list)}
        result["root_states"] = states

        n = len(states)
        labels = list(states.keys())
        values = list(states.values())
        h_start = random.random()
        colors = sns.husl_palette(n_colors=n, h=h_start, s=random.uniform(0.4, 0.7), l=random.uniform(0.5, 0.8))
        hex_colors = [f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}" for r, g, b in colors]

        # Plot
        matplotlib.rcParams['font.family'] = 'Arial'
        matplotlib.rcParams['pdf.fonttype'] = 42

        labels = list(states.keys())
        values = list(states.values())
        sorted_pairs = sorted(zip(values, labels, hex_colors), reverse=True)
        values_s, labels_s, colors_s = zip(*sorted_pairs)

        if chart_type == "pie":
            fig, ax = plt.subplots(figsize=(8, 8))
            wedges, texts = ax.pie(values_s, labels=labels_s, colors=colors_s,
                                   autopct='%1.1f%%', startangle=90)
            ax.set_title("Root State Posterior Probability", fontsize=16, fontweight='bold')
        else:
            fig, ax = plt.subplots(figsize=(max(8, n * 1.2), 6))
            bars = ax.bar(range(n), values_s, color=colors_s, edgecolor='white', linewidth=1.5)
            ax.set_xticks(range(n))
            ax.set_xticklabels(labels_s, rotation=45, ha='right', fontsize=11)
            ax.set_ylabel("Posterior Probability", fontsize=13)
            ax.set_title("Root State Posterior Probability", fontsize=16, fontweight='bold')
            ax.set_ylim(0, max(values_s) * 1.15)
            for bar, val in zip(bars, values_s):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                        f'{val:.3f}', ha='center', va='bottom', fontsize=10)

        plot_path = os.path.join(output_dir, f"RSPP_{chart_type}.pdf")
        save_both(fig, plot_path, bbox_inches='tight', dpi=150, format='pdf')
        plt.close(fig)
        result["plot"] = plot_path
        result["success"] = True
        log.emit(f"RSPP-Viz plot saved: {plot_path}")
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"RSPP-Viz failed: {e}")
        return result


# ═══════════════════════════════════════════════════════════════════
# Module 4: GeoSubsampler
# ═══════════════════════════════════════════════════════════════════

def run_geo_subsampler(
    fasta_file: str,
    output_dir: str,
    num_seqs: int = 50,
    equal_sampling: bool = False,
    replicates: int = 1,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    地理均衡子抽样。

    Parameters
    ----------
    fasta_file : str
        FASTA 文件，序列名需含地点信息
    output_dir : str
        输出目录
    num_seqs : int
        抽取序列数（非等量模式）
    equal_sampling : bool
        是否每个地区等量抽样
    replicates : int
        bootstrap 重复次数（等量模式）
    log : LogCollector, optional

    Returns
    -------
    dict with keys: success, output_files, error
    """
    import random
    import re
    from collections import OrderedDict

    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "output_files": [], "error": None}

    try:
        # Parse FASTA
        dict_fas = OrderedDict()
        with open(fasta_file) as f:
            line = f.readline()
            while line:
                while line and not line.startswith('>'):
                    line = f.readline()
                if not line: break
                fas_name = line.strip()
                fas_seq = ""
                line = f.readline()
                while line and not line.startswith('>'):
                    fas_seq += re.sub(r"\s", "", line)
                    line = f.readline()
                dict_fas[fas_name] = fas_seq

        if equal_sampling:
            site_dict = OrderedDict()
            for fas_name, fas_seq in dict_fas.items():
                match = re.search(r"(?:^|>|_)([A-Za-z0-9]+)(?:_|$)", fas_name)
                if not match:
                    log.warning(f"Cannot parse site from: {fas_name}")
                    continue
                site = match.group(1)
                site_dict.setdefault(site, []).append((fas_name, fas_seq))

            min_count = min(len(seqs) for seqs in site_dict.values())
            if min_count == 0:
                result["error"] = "Some sites have zero sequences"
                return result

            for i in range(replicates):
                sampled = []
                for site, seqs in site_dict.items():
                    sampled.extend(random.sample(seqs, min_count))
                out_name = f"extract_rep{i+1}.fas" if replicates > 1 else "extract.fas"
                out_path = os.path.join(output_dir, out_name)
                with open(out_path, 'w') as fo:
                    fo.write('\n'.join([f"{n}\n{s}" for n, s in sampled]))
                result["output_files"].append(out_path)
        else:
            items = list(dict_fas.items())
            sampled = random.sample(items, min(num_seqs, len(items)))
            out_path = os.path.join(output_dir, "extract.fas")
            with open(out_path, 'w') as fo:
                fo.write('\n'.join([f"{n}\n{s}" for n, s in sampled]))
            result["output_files"].append(out_path)

        result["success"] = True
        log.emit(f"GeoSubsampler: {len(result['output_files'])} file(s) produced")
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"GeoSubsampler failed: {e}")
        return result


# ═══════════════════════════════════════════════════════════════════
# Module 5: BSP-Viz (Bayesian Skyline Plot)
# ═══════════════════════════════════════════════════════════════════

def _natural_key(col):
    """popSize2 < popSize10 (按列名尾部的整数排序)。"""
    m = re.search(r'(\d+)$', str(col))
    return (int(m.group(1)) if m else 0, str(col))


def _parse_beast_log(log_file, log):
    """读 BEAST log → (states, DataFrame)。跳过 # 注释行, state 列缺失时为 None。"""
    import numpy as np
    rows = []
    header = None
    with open(log_file, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if header is None:
                header = parts
                continue
            rows.append(parts)
    if header is None or not rows:
        return None, None
    import pandas as pd
    df = pd.DataFrame(rows, columns=header)
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    states = df['state'].astype('int64').tolist() if 'state' in df.columns else None
    return states, df


def _parse_beast_trees(trees_file, log):
    """读 BEAST .trees → [(state|None, Bio.Phylo tree), ...]。"""
    from Bio import Phylo
    out = []
    pat = re.compile(r'tree\s+STATE_(\d+)', re.IGNORECASE)
    with open(trees_file, encoding='utf-8', errors='replace') as f:
        for line in f:
            m = pat.search(line)
            if 'tree' not in line.lower():
                continue
            # 取首个 '(' 到末尾 ';', 剥掉 Newick 注释段 [&...]
            start = line.find('(')
            if start < 0:
                continue
            nwk = line[start:line.rfind(';') + 1]
            if not nwk.endswith(';'):
                continue
            nwk = re.sub(r'\[&[^\]]*\]', '', nwk)
            nwk = re.sub(r'\[[^\]]*\]', '', nwk)
            try:
                tr = Phylo.read(io.StringIO(nwk), 'newick')
            except Exception as e:
                log.warning(f"[BSP] 树解析失败 (state="
                            f"{m.group(1) if m else '?'}): {e}")
                continue
            out.append((int(m.group(1)) if m else None, tr))
    return out


def _coalescent_ages(tree):
    """树 → (升序聚合事件年龄, 根高)。年龄 = 距今 (root 最老, tips=0)。"""
    depths = tree.depths()
    terms = tree.get_terminals()
    if not terms:
        return [], 0.0
    max_depth = max(depths.get(t, 0.0) for t in terms)
    internals = [c for c in tree.get_nonterminals()]
    ages = sorted(max_depth - depths.get(c, 0.0) for c in internals)
    return ages, max_depth


def run_bsp_viz(
    log_file: str,
    output_dir: str,
    burnin_pct: float = 10.0,
    rscript_path: str = "Rscript",          # 兼容旧签名; 本实现为纯 Python, 不再用 R
    log: Optional[LogCollector] = None,
    trees_file: Optional[str] = None,
) -> Dict:
    """
    Bayesian Skyline Plot —— Tracer 式正确重建 (2026-10-14 重写, 修复 P1)。

    口径 (与 Tracer / Drummond et al. 2005 一致), 对每个后验样本 k:
      · 从 .trees 的第 k 棵树取聚合事件年龄 (升序, 距今);
      · 用该样本 log 行的 groupSize 把聚合区间自今向古划成 N 组;
      · 组 j 的上边界年龄 = 覆盖完前 j 组时的聚合年龄, 组的 Neτ = popSize_j。
    跨样本: 逐组取中位时间 + Neτ 中位数 + 95% HPD (最短区间, utils.ess.calculate_95hpd)。
    图: 0→根高的阶梯曲线 + 95% HPD 带 (log-y)。

    历史坑 (旧实现已废): 按 unique(rootHeight) 分箱、箱内取"5 组跨组"分位数 ——
    x 轴只盖住根高后验区间 (GCVA 实测 28–46 ybp), HPD 是样本内组间离散度,
    既非任何组的后验、也不覆盖近期群体历史 (reports/UPSTREAM_CONSISTENCY_20261014.md §3.1)。

    依赖: log 须同时含 popSize 与 groupSize 列 (beast1_bridge 2026-10-14 起
    fileLog 两者都记; 旧 XML 的输出缺 groupSize → 拒出图并给出重跑指引)。

    Parameters
    ----------
    log_file : str — BEAST .log (含 state 列更佳)
    trees_file : str — BEAST .trees (缺省时探测 log 同名 .trees)
    output_dir : str
    burnin_pct : float — log 与 trees 同比例去 burn-in
    rscript_path : 兼容旧签名, 未使用
    log : LogCollector

    Returns
    -------
    dict with keys: success, plot, summary_tsv, n_samples, n_groups,
        root_height_range, error
    """
    import numpy as np
    try:
        from utils.ess import calculate_95hpd
    except ImportError:                                    # 脚本方式直跑
        from ess import calculate_95hpd

    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "plot": None, "summary_tsv": None,
              "n_samples": 0, "n_groups": 0, "root_height_range": None,
              "error": None}
    output_file = os.path.join(output_dir, "BSP.pdf")

    def _fail(msg):
        result["error"] = msg
        log.error(f"BSP-Viz: {msg}")
        # 退化守卫: 本次没算出可用图 → 不留上一次的旧图冒充新结果
        if os.path.exists(output_file):
            try:
                os.remove(output_file)
                log.warning(f"BSP-Viz: 删除陈旧图 {output_file}")
            except OSError as _e:
                log.warning(f"BSP-Viz: 陈旧图删除失败 {output_file}: {_e}")
        return result

    if not os.path.exists(log_file):
        return _fail(f"log 不存在: {log_file}")

    # ── 1. log: popSize / groupSize 列 ──
    log_states, df = _parse_beast_log(log_file, log)
    if df is None:
        return _fail(f"log 无可用数据行: {log_file}")
    pop_cols = sorted([c for c in df.columns
                       if re.search(r'popSize|pop\.size', c, re.I)],
                      key=_natural_key)
    grp_cols = sorted([c for c in df.columns
                       if re.search(r'groupSize|group\.size', c, re.I)],
                      key=_natural_key)
    if not pop_cols:
        return _fail("log 中找不到 popSize 列 (非 skyline log?)")
    if not grp_cols:
        return _fail("log 中找不到 groupSize 列 —— 组边界不可恢复, 无法做正确的 "
                     "skyline 重建。请在 BEAST1 XML 的 fileLog 里加 "
                     "<parameter idref=\"skyline.groupSize\"/> 后重跑 "
                     "(beast1_bridge 2026-10-14 起已自动带上), 或改用 Tracer "
                     "的 Bayesian Skyline Reconstruction 导出 TSV。")
    n_groups = len(pop_cols)
    if len(grp_cols) != n_groups:
        return _fail(f"popSize({n_groups}) 与 groupSize({len(grp_cols)}) 列数不一致")
    n_burn = int(len(df) * burnin_pct / 100.0)
    df = df.iloc[n_burn:].reset_index(drop=True)
    if log_states is not None:
        log_states = log_states[n_burn:]        # 与 df 同步去 burn-in (行号才对得上)
    log.emit(f"[BSP] log: {len(df)} 样本 (burn-in {burnin_pct:.0f}% = {n_burn} 行), "
             f"{n_groups} 组")

    # ── 2. trees: 与 log 按样本对齐 ──
    if not trees_file or not os.path.exists(trees_file):
        _cand = os.path.join(os.path.dirname(log_file),
                             os.path.splitext(os.path.basename(log_file))[0]
                             + '.trees')
        trees_file = _cand if os.path.exists(_cand) else None
    if not trees_file:
        return _fail("找不到 .trees 后验树文件 (组时间边界必须从后验树取)")
    samples = _parse_beast_trees(trees_file, log)
    if not samples:
        return _fail(f".trees 无可解析样本: {trees_file}")
    t_burn = int(len(samples) * burnin_pct / 100.0)
    samples = samples[t_burn:]
    tree_states = [s for s, _t in samples]
    trees = [t for _s, t in samples]
    if tree_states and tree_states[0] is None:
        tree_states = [None] * len(trees)
    if log_states is not None and any(s is not None for s in tree_states):
        log_map = {s: i for i, s in enumerate(log_states)}
        pairs = [(log_map[s], i) for i, s in enumerate(tree_states)
                 if s in log_map]
        if len(pairs) < max(3, len(df) // 2):
            log.warning(f"[BSP] state 对齐仅命中 {len(pairs)} 对, 退回按位置对齐")
            pairs = [(i, i) for i in range(min(len(df), len(trees)))]
    else:
        pairs = [(i, i) for i in range(min(len(df), len(trees)))]
    if len(pairs) < 3:
        return _fail(f"log×trees 有效配对不足 ({len(pairs)} < 3)")
    log.emit(f"[BSP] trees: {len(trees)} 棵 (配对 {len(pairs)})")

    # ── 3. 逐样本重建组边界 ──
    grp_heights = [[] for _ in range(n_groups)]      # 每组的边界年龄 (跨样本)
    grp_pops = [[] for _ in range(n_groups)]         # 每组的 Neτ (跨样本)
    skipped = 0
    for li, ti in pairs:
        row = df.iloc[li]
        gsizes = [int(round(float(row[c]))) for c in grp_cols]
        pops = [float(row[c]) for c in pop_cols]
        coal, _root_h = _coalescent_ages(trees[ti])
        if sum(gsizes) != len(coal) or any(g <= 0 for g in gsizes):
            skipped += 1
            continue
        cum = 0
        for j, g in enumerate(gsizes):
            cum += g
            grp_heights[j].append(coal[cum - 1])     # 组的上边界年龄
            grp_pops[j].append(pops[j])
    if skipped:
        log.warning(f"[BSP] {skipped} 个样本 groupSize 与树聚合区间数不匹配, 已跳过")
    if len(grp_pops[0]) < 3:
        return _fail(f"有效样本不足 ({len(grp_pops[0])} < 3), 无法重建 skyline")

    # ── 4. 逐组统计 ──
    summary_rows = []
    med_h, med_n, lo_n, hi_n = [], [], [], []
    for j in range(n_groups):
        h = np.asarray(grp_heights[j], dtype=float)
        p = np.asarray(grp_pops[j], dtype=float)
        h_med = float(np.median(h))
        p_lo, p_hi = calculate_95hpd(p)
        p_med = float(np.median(p))
        med_h.append(h_med)
        med_n.append(p_med)
        lo_n.append(float(p_lo))
        hi_n.append(float(p_hi))
        summary_rows.append({
            "group": j + 1, "time_median": h_med,
            "Ne_tau_median": p_med,
            "Ne_tau_hpd_lower": float(p_lo), "Ne_tau_hpd_upper": float(p_hi),
            "n_samples": len(p)})
    if not all(np.isfinite(med_n)) or len(set(np.round(med_n, 12))) < 2:
        return _fail("退化的 skyline (非有限或零方差中位数) → 不产图")
    result["n_samples"] = len(grp_pops[0])
    result["n_groups"] = n_groups
    result["root_height_range"] = [float(min(med_h)), float(max(med_h))]

    # 汇总 TSV 落盘 (VirPhyKit BSP 模块期望的 Time/Median/Lower/Upper 同构产物)
    tsv_path = os.path.join(output_dir, "BSP_summary.tsv")
    try:
        pd.DataFrame(summary_rows).to_csv(tsv_path, sep="\t", index=False)
        result["summary_tsv"] = tsv_path
        log.emit(f"[BSP] 汇总表 → {tsv_path}")
    except OSError as _e:
        log.warning(f"[BSP] 汇总表落盘失败 (非阻断): {_e}")

    # ── 5. 绘图: 0→根高阶梯曲线 + 95% HPD 带 ──
    # OKABE 蓝系 (与 LTT/RTT 面板一致); dpi 300; 图例与量纲齐全。
    apply_sci_style(base_size=13)
    xs = [0.0] + med_h                    # n+1 个边界: 组 j 覆盖 (h_{j-1}, h_j]
    ys = med_n + [med_n[-1]]              # step='post' 需与 xs 等长
    los = lo_n + [lo_n[-1]]
    his = hi_n + [hi_n[-1]]
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.fill_between(xs, los, his, step='post', alpha=0.30,
                    color="#56B4E9", label="95% HPD")
    ax.step(xs, ys, where='post', linewidth=2.0, color="#0072B2",
            label="Posterior median")
    ax.set_xlim(0, max(med_h) * 1.03)
    ax.set_yscale('log')
    ax.set_xlabel("Time (years before present)", fontsize=14)
    ax.set_ylabel("Effective population size (" + r"$N_e\tau$" + ")",
                  fontsize=14)
    ax.set_title("Bayesian skyline plot", fontsize=16, fontweight='bold')
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=11, loc="best")
    despine(ax)
    save_both(fig, output_file, bbox_inches='tight', dpi=300, format='pdf')
    plt.close(fig)

    result["plot"] = output_file
    result["success"] = True
    log.emit(f"[BSP] {n_groups} 组 × {result['n_samples']} 样本, "
             f"根高中位 {med_h[-1]:.2f} ybp → {output_file}")
    return result


# ═══════════════════════════════════════════════════════════════════
# Module 6: MJRM Generator (Markov Jumps & Rewards Matrix)
# ═══════════════════════════════════════════════════════════════════

def generate_mjrm(traits: str, output_dir: str, input_xml: Optional[str] = None,
                  output_name: str = "beast_with_mjrm.xml") -> Dict:
    """
    为 BEAST XML 生成 Markov Jumps 和 Rewards Matrix。

    Parameters
    ----------
    traits : str
        逗号分隔的离散性状，如 "China,USA,Japan"
    output_dir : str
        输出目录
    input_xml : str, optional
        BEAUti 1.x XML 模板路径。不提供则输出纯文本矩阵
    output_name : str
        输出文件名

    Returns
    -------
    dict with keys: success, output_file, error
    """
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "output_file": None, "error": None}

    # 自研 MJRM 配置生成 (等效 VirPhyKit MakovMJump.generate_config)
    try:
        from utils.mjrm_config import generate_config
        success, msg = generate_config(traits, output_dir, output_name, input_xml)
        if success:
            result["output_file"] = os.path.join(output_dir, output_name)
            result["success"] = True
        else:
            result["error"] = msg
        return result
    except Exception as e:
        result["error"] = str(e)
        return result


# ═══════════════════════════════════════════════════════════════════
# Module 7: Geography — Mantel Test (遗传-地理距离关联)
# ═══════════════════════════════════════════════════════════════════

def run_mantel_test(
    fasta_file: str,
    metadata_csv: str,
    output_dir: str,
    n_permutations: int = 9999,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    Mantel 检验：地理距离 vs 遗传距离的相关性。

    Parameters
    ----------
    fasta_file : str
        FASTA 多序列比对 (序列名 = 样本名)
    metadata_csv : str
        CSV 含 name, location 列
    output_dir : str
    n_permutations : int
        置换检验次数
    log : LogCollector, optional

    Returns
    -------
    dict: {success, mantel_r, p_value, geo_dist_matrix, gen_dist_matrix, plot}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "mantel_r": None, "p_value": None,
              "geo_dist_matrix": None, "gen_dist_matrix": None, "plot": None, "error": None}

    try:
        import numpy as np
        from scipy.spatial.distance import squareform, pdist
        from geopy.distance import geodesic
        from geopy.geocoders import Nominatim
    except ImportError:
        result["error"] = "Missing scipy or geopy. Run: pip install scipy geopy"
        log.error(result["error"])
        return result

    # ── 1. Parse sample locations ──
    sample_locs = {}
    try:
        with open(metadata_csv, encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = row.get('name', '').strip()
                loc = row.get('location', '').strip()
                if name and loc:
                    sample_locs[name] = loc
    except Exception as e:
        result["error"] = f"Failed to read metadata: {e}"
        return result

    # ── 2. Parse sequences ──
    seqs = {}
    try:
        for rec in SeqIO.parse(fasta_file, "fasta"):
            seqs[rec.id] = str(rec.seq).upper()
    except Exception as e:
        result["error"] = f"Failed to read fasta: {e}"
        return result

    # Intersection: samples with both sequence and location
    common = sorted(set(seqs.keys()) & set(sample_locs.keys()))
    if len(common) < 5:
        result["error"] = f"Only {len(common)} samples have both sequence and location (need ≥5)"
        log.error(result["error"])
        return result
    log.emit(f"Mantel test: {len(common)} samples with sequence + location")

    # ── 3. Compute genetic distance matrix (p-distance) ──
    seq_list = [seqs[name] for name in common]
    n = len(common)
    gen_dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            # 2026-09-15 修复: 统一走 p_distance_pair (跳过 gap, N 按 IUPAC 简并)
            d = p_distance_pair(seq_list[i], seq_list[j])
            gen_dist[i, j] = d
            gen_dist[j, i] = d
    gen_dist_condensed = squareform(gen_dist)

    # ── 4. Compute geographic distance matrix ──
    # Real coordinate lookup for Chinese locations (lat, lon)
    COORD_LOOKUP = {
        'China, Ningxia, Yinchuan': (38.47, 106.27),
        'China, Ningxia, Zhongning': (37.49, 105.68),
        'China, Neimenggu': (40.82, 111.75),
        'China, Beijing, Beijing': (39.90, 116.40),
        'China, Beijing': (39.90, 116.40),
        'China, Beijing, Unknown': (39.90, 116.40),
        'China, Guangdong, Guangzhou': (23.13, 113.26),
        'China, Qinghai, Xining': (36.62, 101.78),
        'China, Gansu, Wuwei': (37.93, 102.64),
        'China, Liaoning, Shenyang': (41.80, 123.43),
        'China, ningxia, Yinchuan': (38.47, 106.27),
        # ── 裸省份名 (PSTVd metadata 格式, 2026-09-02 补) ──
        'Ningxia': (38.47, 106.27),          # 宁夏 (银川)
        'Gansu': (36.06, 103.83),            # 甘肃 (兰州)
        'Qinghai': (36.62, 101.78),          # 青海 (西宁)
        'Neimenggu': (40.82, 111.75),        # 内蒙古 (呼和浩特)
        'Inner Mongolia': (40.82, 111.75),
        'Beijing': (39.90, 116.40),
        'Guangdong': (23.13, 113.26),
        'Jiangsu': (32.06, 118.78),          # 江苏 (南京)
        'Shandong': (36.67, 117.02),
        'Hebei': (38.04, 114.51),
        'Xinjiang': (43.79, 87.63),
        'Yunnan': (25.04, 102.71),
        'Sichuan': (30.66, 104.07),
        'Hunan': (28.23, 112.94),
        'Hubei': (30.59, 114.31),
        'Henan': (34.75, 113.63),
        'Shaanxi': (34.27, 108.95),
        'Shanxi': (37.87, 112.55),
        'Liaoning': (41.80, 123.43),
        'Fujian': (26.08, 119.30),
        'Zhejiang': (30.27, 120.16),
        'Anhui': (31.86, 117.28),
        'Heilongjiang': (45.80, 126.53),
        'Jilin': (43.88, 125.32),
        'Guangxi': (22.82, 108.32),
        'Guizhou': (26.65, 106.63),
        'Hainan': (20.04, 110.32),
        'Chongqing': (29.56, 106.55),
        'Tianjin': (39.08, 117.20),
        'Shanghai': (31.23, 121.47),
        'Tibet': (29.65, 91.13),
    }
    import math
    def haversine_km(lat1, lon1, lat2, lon2):
        R = 6371
        dlat = math.radians(lat2 - lat1)
        dlon = math.radians(lon2 - lon1)
        a = math.sin(dlat/2)**2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon/2)**2
        return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))

    geo_dist = np.zeros((n, n))
    coords = {}
    skipped_locs = set()
    # 分层地名解析 (2026-09-02, geo_resolver): 旧 COORD_LOOKUP 作 L0 遗留精确层,
    # 其后走 省/市中英文名+拼音变体 → 段拆分 → 模糊匹配 → 缓存 → (可选)在线
    try:
        from utils.geo_resolver import get_resolver
        _gres = get_resolver()
    except ImportError:
        import importlib.util as _ilu
        _spec = _ilu.spec_from_file_location(
            'geo_resolver', os.path.join(os.path.dirname(__file__), 'geo_resolver.py'))
        _mod = _ilu.module_from_spec(_spec)
        _spec.loader.exec_module(_mod)
        _gres = _mod.get_resolver()
    _src_stats = {}
    for name in common:
        loc_str = sample_locs[name]
        coord = None
        if loc_str in COORD_LOOKUP:
            coord = COORD_LOOKUP[loc_str]
            _src_stats['legacy'] = _src_stats.get('legacy', 0) + 1
        else:
            coord, src = _gres.resolve(loc_str)
            if coord:
                _src_stats[src] = _src_stats.get(src, 0) + 1
        if coord:
            coords[name] = coord
        else:
            # 解析失败: 不伪造坐标 (假坐标会制造虚假地理距离), 排除并列明
            skipped_locs.add(loc_str)
    _gres.flush()
    if _src_stats:
        log.emit(f"[Mantel] 坐标解析来源: {_src_stats}")
    if skipped_locs:
        log.emit(f"[Mantel] WARNING: {len(skipped_locs)} unrecognized location(s), "
                 f"samples excluded: {sorted(skipped_locs)}")
        log.emit("[Mantel] Add these to utils/geo_resolver.py (PROVINCE/CITY_LOOKUP) or cache/geo_cache.json to include them.")
        common = [c for c in common if c in coords]
        n = len(common)
        if n < 3:
            result["error"] = (f"Mantel test aborted: only {n} samples with known "
                               f"coordinates (need >=3). Unknown: {sorted(skipped_locs)}")
            log.emit(f"[Mantel] {result['error']}")
            return result
        # 过滤样本后重算遗传距离子矩阵 (保序, 简单可靠)
        seq_list = [seqs[name] for name in common]
        gen_dist = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                # 2026-09-15 修复: 统一走 p_distance_pair (跳过 gap, N 按 IUPAC 简并)
                gen_dist[i, j] = gen_dist[j, i] = p_distance_pair(seq_list[i], seq_list[j])
        gen_dist_condensed = squareform(gen_dist)
        # geo_dist 也必须用新 n 重建 (2026-09-02 崩溃修复:
        # 旧代码用旧 n 的 geo_dist 直接算 → 维度不匹配 ValueError)
        geo_dist = np.zeros((n, n))
        log.emit(f"[Mantel] {n} samples retained after coordinate filtering")

    for i in range(n):
        for j in range(i + 1, n):
            lat1, lon1 = coords[common[i]]
            lat2, lon2 = coords[common[j]]
            d = haversine_km(lat1, lon1, lat2, lon2)
            geo_dist[i, j] = d
            geo_dist[j, i] = d
    geo_dist_condensed = squareform(geo_dist)

    # ── 5. Mantel test ──
    r_obs = np.corrcoef(gen_dist_condensed, geo_dist_condensed)[0, 1]

    # Permutation
    perm_vals = []
    rng = np.random.RandomState(42)
    for _ in range(n_permutations):
        permuted = rng.permutation(gen_dist_condensed)
        perm_vals.append(np.corrcoef(permuted, geo_dist_condensed)[0, 1])
    perm_vals = np.array(perm_vals)

    p_value = (np.sum(np.abs(perm_vals) >= np.abs(r_obs)) + 1) / (n_permutations + 1)

    result["mantel_r"] = float(r_obs)
    result["p_value"] = float(p_value)
    result["n_samples"] = n
    log.emit(f"Mantel r = {r_obs:.4f}, p = {p_value:.4f} (n={n}, {n_permutations} perms)")

    # ── 6. Plot ──
    # 2026-10-12: 配色统一 OKABE（原 #00A0E9 / #E9522E 不在调色板内）；热图
    # YlOrRd → 项目连续色 viridis；dpi 150 → 300。
    # 画布 14 in → 9.32 in：T 面板在 Fig 3 里占 2 栏（4.58 in），源 14 in 缩放
    # 0.33 倍会把 13 pt 标题压到 4.3 pt。改成 9.32 in 后缩放 0.49，标题 ≈6.4 pt。
    apply_sci_style(base_size=11)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.32, 4.6))

    # Scatter
    ax1.scatter(geo_dist_condensed, gen_dist_condensed, alpha=0.5, s=12,
                color=OKABE_ITO[4], edgecolors='none',
                label=f"{len(geo_dist_condensed):,} sample pairs")
    z = np.polyfit(geo_dist_condensed, gen_dist_condensed, 1)
    x_line = np.linspace(0, geo_dist_condensed.max(), 100)
    ax1.plot(x_line, z[0]*x_line + z[1], '--', color=OKABE_ITO[5], linewidth=1.6,
             label="least-squares fit")
    ax1.set_xlabel("Geographic distance (km)", fontsize=13)
    ax1.set_ylabel("Genetic distance (p-distance)", fontsize=13)
    # 2026-10-12（审稿人视角）：Mantel 的 P 值在"配对数极大"时极易显著，
    # 必须同时给出效应量（r）与配对数，否则读者会把 r≈0 误读成"存在地理隔离"。
    # ⚠ 配对数**只写在图例里**：标题里再写一遍会与右图的
    #   "Genetic distance matrix  (n = ...)" 标题在 4.58 in 半栏宽下相撞
    #   （实测 "(7,140 pairs)" 与 "Genetic" 直接叠字），且属信息重复。
    ax1.set_title(f"Mantel test:  $r$ = {r_obs:.3f},  $P$ = {p_value:.4f}",
                  fontsize=13, fontweight='bold')
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=9, loc="upper left", frameon=False)
    despine(ax1)

    # Genetic distance heatmap
    # 2026-10-12（审稿人视角）：旧版在 127×127 的矩阵上等间隔抽 28 个样本名做刻度，
    # 抽出来的又都是被截断到 8 字符的 accession 前缀（CRR11261…）—— 既糊又无信息量，
    # 而且同一 accession 的多个样本会重复出现，读者会误以为图错。
    # 现直接去掉刻度标签，把维度写进标题（矩阵本身的可读性不依赖逐样本标签）。
    im = ax2.imshow(gen_dist, cmap=SEQ_CMAP, aspect='auto',
                    interpolation='nearest')
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_title(f"Genetic distance matrix  (n = {n})",
                  fontsize=13, fontweight='bold')
    plt.colorbar(im, ax=ax2, shrink=0.8, label='p-distance')

    plt.tight_layout()
    plot_path = os.path.join(output_dir, "mantel_test.pdf")
    save_both(fig, plot_path, bbox_inches='tight', dpi=300, format='pdf')
    plt.close(fig)
    result["plot"] = plot_path

    # Save matrices
    geo_mat_path = os.path.join(output_dir, "geo_distance.csv")
    gen_mat_path = os.path.join(output_dir, "genetic_distance.csv")
    np.savetxt(geo_mat_path, geo_dist, delimiter=',', header=','.join(common), comments='', fmt='%.1f')
    np.savetxt(gen_mat_path, gen_dist, delimiter=',', header=','.join(common), comments='', fmt='%.6f')
    result["geo_dist_matrix"] = geo_mat_path
    result["gen_dist_matrix"] = gen_mat_path

    result["success"] = True
    return result


# ═══════════════════════════════════════════════════════════════════
# Module 8: Geography — Tree Geo-Coloring (树上地理着色)
# ═══════════════════════════════════════════════════════════════════

def run_tree_geography(
    tree_file: str,
    metadata_csv: str,
    fasta_file: str,
    output_dir: str,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    在系统发育树上按地理来源给叶子着色。

    Parameters
    ----------
    tree_file : str
        Newick 树文件
    metadata_csv : str
        CSV 含 name, location
    fasta_file : str
        FASTA 多序列比对 (用于提取叶子名称)
    output_dir : str
    log : LogCollector, optional

    Returns
    -------
    dict: {success, plot, n_regions, region_counts}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "plot": None, "n_regions": 0,
              "region_counts": {}, "error": None}

    try:
        from Bio import Phylo
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
    except ImportError as e:
        result["error"] = f"Import error: {e}"
        return result

    # ── 1. Load tree ──
    tree = Phylo.read(tree_file, "newick")
    tips = {clade.name for clade in tree.get_terminals()}
    log.emit(f"Tree has {len(tips)} tips")

    # ── 2. Load location map ──
    sample_locs = {}
    try:
        with open(metadata_csv, encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = row.get('name', '').strip()
                loc = row.get('location', '').strip()
                if name and loc:
                    # Simplify: keep province-level detail
                    parts = [p.strip() for p in loc.replace(':', ',').split(',')]
                    if len(parts) >= 2:
                        simple = parts[1]  # e.g. 'Ningxia'
                    else:
                        simple = parts[0]
                    sample_locs[name] = simple
    except Exception as e:
        result["error"] = f"Failed to read metadata: {e}"
        return result

    # ── 3. Color mapping ──
    regions = sorted(set(sample_locs.values()))
    if len(regions) < 2:
        result["error"] = f"Only {len(regions)} unique region(s), need ≥2 for coloring"
        log.error(result["error"])
        return result

    n_regions = len(regions)
    cmap = plt.get_cmap('tab20', n_regions) if n_regions <= 20 else plt.get_cmap('gist_rainbow', n_regions)
    region_colors = {r: cmap(i) for i, r in enumerate(regions)}
    result["n_regions"] = n_regions
    result["region_counts"] = {r: sum(1 for v in sample_locs.values() if v == r) for r in regions}

    log.emit(f"{n_regions} regions: {result['region_counts']}")

    # ── 4. Draw tree with colored tips ──
    fig, ax = plt.subplots(figsize=(14, max(8, len(tips) * 0.2)))

    def label_func(clade):
        name = clade.name
        if name and name in sample_locs:
            region = sample_locs[name]
            color = region_colors[region]
            # Color the label
            return name
        return name if name else ""

    # Use Phylo.draw with custom coloring
    Phylo.draw(tree, axes=ax, do_show=False,
               label_func=lambda c: "",  # Don't show labels in default draw
               branch_labels=lambda c: "",
               show_confidence=False)

    # Manually add colored tip labels
    from matplotlib import transforms
    ax.set_xlim(ax.get_xlim())
    y_positions = {}
    # 历史坑: 等距估算 y + 枝长当 x → 标签全部错位, 且 except: pass 静默吞异常。
    # 改用 Bio.Phylo 官方布局坐标 (DepthFirstTraverser 同款算法: y=叶序号, x=到根累积深度)。
    from Bio.Phylo.Newick import Clade as _Clade
    _depths = tree.depths()
    _tips = tree.get_terminals()
    for i, clade in enumerate(_tips):
        if not (clade.name and clade.name in sample_locs):
            continue
        region = sample_locs[clade.name]
        color = region_colors[region]
        y = float(i)
        x = float(_depths.get(clade, 0.0))
        ax.annotate(clade.name, xy=(x, y), fontsize=8,
                   color=color, fontweight='bold',
                   xytext=(5, 0), textcoords='offset points')

    # Legend
    legend_patches = [Patch(color=region_colors[r], label=f"{r} (n={result['region_counts'][r]})")
                      for r in regions]
    ax.legend(handles=legend_patches, loc='upper right', fontsize=9,
              title="Geographic Origin", title_fontsize=10)
    ax.set_title("Phylogeny colored by geographic origin", fontsize=16, fontweight='bold')

    plot_path = os.path.join(output_dir, "tree_geography.pdf")
    save_both(fig, plot_path, bbox_inches='tight', dpi=150, format='pdf')
    plt.close(fig)
    result["plot"] = plot_path
    result["success"] = True
    log.emit(f"Tree geo-coloring: {plot_path}")
    return result


# ═══════════════════════════════════════════════════════════════════
# Module 9: Geography — VirSpaceTime (时空分布可视化)
# ═══════════════════════════════════════════════════════════════════

VIRSPACETIME_R_SCRIPT = '''
library(ggplot2)
library(scales)
library(ggsci)
library(tidyr)

args <- commandArgs(trailingOnly=TRUE)
input_csv <- args[1]
output_file <- args[2]

# 2026-10-12（审稿人视角修复，四处）：
#  ① `check.names=FALSE` —— 默认 TRUE 会用 make.names() 把表头 "China, Beijing"
#     改写成 "China..Beijing"（空格/逗号都变点），图例里于是出现双点标签；
#  ② 显示名统一取**最后一段**（最具体地名），与 Python 侧
#     `sample_plot._canonical_label` 同口径 —— 否则 V 的 "China.Beijing.Beijing"
#     和 A/B 的 "Beijing" 对不上；
#  ③ 配色改为"按地名排序后依次取 OKABE_ITO"，与 `sample_plot.loc_color_map`
#     完全一致 —— 否则同一地点在 A/B 与 V 里颜色不同（多面板图的硬伤）；
#  ④ 图例原为单行 7 项、base_size=15 → 右侧被裁掉；现压到 2 行、缩小字号，
#     并给气泡图补"大小 = 样本数"的图例。
data <- read.csv(input_csv, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

loc_cols <- names(data)[-1]

short_label <- function(x) {
    parts <- trimws(strsplit(x, ",")[[1]])
    parts <- parts[nzchar(parts)]
    parts[length(parts)]
}

df <- pivot_longer(data, cols=all_of(loc_cols), names_to="Location", values_to="Count")
df$Location <- vapply(df$Location, short_label, character(1))
df <- df[df$Count > 0,]

okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
           "#D55E00", "#CC79A7", "#000000", "#999999")
locs <- sort(unique(df$Location))
loc_map <- setNames(rep(okabe, length.out=length(locs)), locs)

# 整数年份刻度，两幅子图共用（旧版上图整数、下图小数，同一张图两种口径）
yr <- sort(unique(df$Year))
step <- max(1, floor(length(yr) / 6))
brk <- yr[seq(1, length(yr), by=step)]

# Plot 1: stacked bars by location
p1 <- ggplot(df, aes(x=Year, y=Count, fill=Location)) +
    geom_bar(stat="identity", position="stack", width=0.8) +
    scale_fill_manual(values=loc_map) +
    scale_x_continuous(breaks=brk) +
    theme_bw(base_size=13) +
    labs(x="Year", y="Number of isolates",
         title="Temporal distribution by location") +
    theme(legend.position="bottom", legend.title=element_blank(),
          legend.text=element_text(size=8),
          legend.key.size=unit(0.28, "cm"),
          panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold", hjust=0.5)) +
    guides(fill=guide_legend(nrow=2, byrow=TRUE))

# Plot 2: dot plot (location vs year, size = count)
p2 <- ggplot(df, aes(x=Year, y=Location)) +
    geom_point(aes(color=Location, size=Count), alpha=0.85) +
    scale_color_manual(values=loc_map, guide="none") +
    scale_x_continuous(breaks=brk) +
    scale_size_continuous(name="isolates", range=c(1.5, 9),
                          breaks=pretty(range(df$Count), n=4)) +
    theme_bw(base_size=13) +
    labs(x="Year", y="", title="Spatiotemporal sample distribution") +
    theme(panel.grid.minor=element_blank(),
          legend.position="right",
          legend.title=element_text(size=8),
          legend.text=element_text(size=7),
          legend.key.size=unit(0.28, "cm"),
          plot.title=element_text(face="bold", hjust=0.5))

# Combine
library(patchwork)
combined <- p1 / p2 + plot_layout(heights=c(2, 1))

# 画布 = Fig 3 的 2 栏宽度量级（源 9.32 in 排进 4.58 in → 0.49 倍，
# 13 pt 落到 ≈6.4 pt，可读）。原 12×9 in + 13 pt 排进去只有 5 pt。
ggsave(output_file, plot=combined, width=9.32, height=7.4, dpi=300)
# 2026-10-16: 图件双格式 (PNG+PDF 同主名) —— output_file 由调用方传 .pdf,
# 这里按同主名补一张 PNG (300 dpi, 白底); PDF 那条参数逐字节不变。
ggsave(sub("\\\\.pdf$", ".png", output_file, ignore.case=TRUE), plot=combined,
       width=9.32, height=7.4, dpi=300, bg="white")
cat("VirSpaceTime plot saved to", output_file, "\\n")
'''


def run_virspacetime(
    metadata_csv: str,
    output_dir: str,
    rscript_path: str = "Rscript",
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    时空分布可视化：按年份 × 地点绘制样本分布。

    Parameters
    ----------
    metadata_csv : str
        CSV 含 name, date, location 列
    output_dir : str
    rscript_path : str
    log : LogCollector, optional

    Returns
    -------
    dict: {success, plot, yearly_table, error}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "plot": None, "yearly_table": None, "error": None}

    try:
        # ── Preprocess: extract year from date, count by location ──
        records = []
        with open(metadata_csv, encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                date_str = row.get('date', '').strip()
                loc = row.get('location', '').strip()
                if date_str and loc:
                    year = date_str[:4]
                    if year.isdigit():
                        records.append({'year': int(year), 'location': loc})

        if not records:
            result["error"] = "No date+location data found"
            return result

        df = pd.DataFrame(records)
        yearly_table = df.pivot_table(index='year', columns='location',
                                      aggfunc='size', fill_value=0)
        yearly_table.index.name = 'Year'

        # Save table
        table_path = os.path.join(output_dir, "spatiotemporal_counts.csv")
        yearly_table.to_csv(table_path)
        result["yearly_table"] = table_path

        # ── Run R script for visualization ──
        output_plot = os.path.join(output_dir, "VirSpaceTime.pdf")
        r_failed = True
        r_script = None
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as f:
                f.write(VIRSPACETIME_R_SCRIPT)
                r_script = f.name
            cmd = [rscript_path, r_script, table_path, output_plot]
            log.emit(f"Running VirSpaceTime (R): {' '.join(cmd)}")
            proc = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=120)
            if proc.returncode == 0 and os.path.exists(output_plot):
                r_failed = False
                result["plot"] = output_plot
            elif proc.returncode != 0:
                log.warning(f"VirSpaceTime R 失败 (rc={proc.returncode}): "
                            f"{(proc.stderr or proc.stdout)[-200:]}")
        except (FileNotFoundError, subprocess.TimeoutExpired) as e:
            # 历史坑: 元组里 Exception 冗余 + 裸 pass 吞掉一切(含 unlink 失败)
            log.warning(f"VirSpaceTime R 异常: {e}")
        finally:
            if r_script and os.path.exists(r_script):
                try:
                    os.unlink(r_script)
                except OSError:
                    pass

        if r_failed:
            log.emit("VirSpaceTime R unavailable, using Python fallback")
            run_virspacetime_python(yearly_table, output_dir, result, log)

        result["success"] = True
        log.emit(f"VirSpaceTime: {result['plot']}")
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"VirSpaceTime failed: {e}")
        return result


def run_virspacetime_python(yearly_table, output_dir, result, log):
    """Python fallback for VirSpaceTime plot"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9),
                                     gridspec_kw={'height_ratios': [2, 1]})

    locations = yearly_table.columns.tolist()
    years = yearly_table.index.tolist()
    n_loc = len(locations)

    cmap = plt.get_cmap('tab20', n_loc) if n_loc <= 20 else plt.get_cmap('gist_rainbow', n_loc)
    colors = [cmap(i) for i in range(n_loc)]

    # Stacked bar
    bottom = None
    for i, loc in enumerate(locations):
        vals = yearly_table[loc].values
        ax1.bar(years, vals, bottom=bottom, color=colors[i], label=loc, width=0.8)
        bottom = vals if bottom is None else bottom + vals
    ax1.set_xlabel("Year", fontsize=13)
    ax1.set_ylabel("Number of isolates", fontsize=13)
    ax1.set_title("Temporal distribution by location", fontsize=14, fontweight='bold')
    ax1.legend(loc='upper left', fontsize=9, ncol=2)

    # Dot plot
    for i, loc in enumerate(locations):
        sizes = yearly_table[loc].values * 20 + 10
        # 历史坑: c=[colors[i]] 是长度 1 的列表, years>1 点时 matplotlib ValueError → fallback 整个失败
        ax2.scatter(years, [i]*len(years), s=sizes, c=colors[i],
                   alpha=0.7, edgecolors='white', linewidth=0.5)
    ax2.set_yticks(range(n_loc))
    ax2.set_yticklabels(locations, fontsize=9)
    ax2.set_xlabel("Year", fontsize=13)
    ax2.set_title("Spatiotemporal sample distribution", fontsize=14, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plot_path = os.path.join(output_dir, "VirSpaceTime.pdf")
    save_both(fig, plot_path, bbox_inches='tight', dpi=150, format='pdf')
    plt.close(fig)
    result["plot"] = plot_path


# ═══════════════════════════════════════════════════════════════════
# Module 10: Host Differentiation (宿主分化分析)
# ═══════════════════════════════════════════════════════════════════

def run_tree_host(
    tree_file: str,
    metadata_csv: str,
    fasta_file: str,
    output_dir: str,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    在系统发育树上按宿主物种着色。

    Parameters
    ----------
    tree_file : str — Newick 树
    metadata_csv : str — CSV 含 name, host 列
    fasta_file : str — FASTA (用于叶子名验证)
    output_dir : str
    log : LogCollector, optional

    Returns
    -------
    dict: {success, plot, n_hosts, host_counts}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "plot": None, "n_hosts": 0,
              "host_counts": {}, "error": None}

    try:
        from Bio import Phylo

        tree = Phylo.read(tree_file, "newick")
        tips = {clade.name for clade in tree.get_terminals()}
        log.emit(f"Tree has {len(tips)} tips")

        # Load host map
        sample_hosts = {}
        with open(metadata_csv, encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = row.get('name', '').strip()
                host = row.get('host', '').strip()
                if name and host:
                    # Remove genus prefix for cleaner labels
                    simple_host = host.replace('Lycium ', 'L. ')
                    sample_hosts[name] = simple_host

        hosts = sorted(set(sample_hosts.values()))
        if len(hosts) < 2:
            result["error"] = f"Only {len(hosts)} host(s), need >=2"
            log.error(result["error"])
            return result

        n_hosts = len(hosts)
        cmap = plt.get_cmap('Set2', n_hosts) if n_hosts <= 8 else plt.get_cmap('tab20', n_hosts)
        host_colors = {h: cmap(i) for i, h in enumerate(hosts)}
        result["n_hosts"] = n_hosts
        result["host_counts"] = {h: sum(1 for v in sample_hosts.values() if v == h) for h in hosts}
        log.emit(f"{n_hosts} hosts: {result['host_counts']}")

        # Draw tree
        fig, ax = plt.subplots(figsize=(14, max(8, len(tips) * 0.2)))
        Phylo.draw(tree, axes=ax, do_show=False,
                   label_func=lambda c: "", branch_labels=lambda c: "",
                   show_confidence=False)

        # Color tips by host (同 tree_geography 修复: 官方布局坐标, 不再静默吞异常)
        _depths_h = tree.depths()
        for i, clade in enumerate(tree.get_terminals()):
            if not (clade.name and clade.name in sample_hosts):
                continue
            host = sample_hosts[clade.name]
            color = host_colors[host]
            y = float(i)
            x = float(_depths_h.get(clade, 0.0))
            ax.annotate(clade.name, xy=(x, y), fontsize=7,
                       color=color, fontweight='bold',
                       xytext=(5, 0), textcoords='offset points')

        from matplotlib.patches import Patch
        legend_patches = [Patch(color=host_colors[h],
                                label=f"{h} (n={result['host_counts'][h]})")
                         for h in hosts]
        ax.legend(handles=legend_patches, loc='upper right', fontsize=9,
                  title="Host Species", title_fontsize=10)
        ax.set_title("Phylogeny colored by host species", fontsize=16, fontweight='bold')

        plot_path = os.path.join(output_dir, "tree_host.pdf")
        save_both(fig, plot_path, bbox_inches='tight', dpi=150, format='pdf')
        plt.close(fig)
        result["plot"] = plot_path
        result["success"] = True
        log.emit(f"Tree host-coloring: {plot_path}")
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"Tree host failed: {e}")
        return result


def run_host_differentiation(
    fasta_file: str,
    metadata_csv: str,
    output_dir: str,
    n_permutations: int = 9999,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    宿主分化定量分析：组间 vs 组内遗传距离 + Permutation 检验。

    Parameters
    ----------
    fasta_file : str — FASTA 多序列比对
    metadata_csv : str — CSV 含 name, host 列
    output_dir : str
    n_permutations : int
    log : LogCollector, optional

    Returns
    -------
    dict: {success, between_mean, within_mean, p_value, plot, host_groups}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "between_mean": None, "within_mean": None,
              "p_value": None, "plot": None, "host_groups": {}, "error": None}

    try:
        import numpy as np

        # Load host map
        sample_hosts = {}
        with open(metadata_csv, encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = row.get('name', '').strip()
                host = row.get('host', '').strip()
                if name and host:
                    sample_hosts[name] = host

        # Load sequences
        seqs = {}
        for rec in SeqIO.parse(fasta_file, "fasta"):
            sid = rec.id
            if sid in sample_hosts:
                seqs[sid] = str(rec.seq).upper()

        # Group by host
        host_groups = {}
        for sid, host in sample_hosts.items():
            if sid in seqs:
                host_groups.setdefault(host, []).append(sid)

        # Filter hosts with >=2 samples
        host_groups = {h: ids for h, ids in host_groups.items() if len(ids) >= 2}
        if len(host_groups) < 2:
            result["error"] = f"Need >=2 host groups with >=2 samples each, got {len(host_groups)}"
            log.error(result["error"])
            return result

        log.emit(f"Host groups: { {h: len(v) for h, v in host_groups.items()} }")
        result["host_groups"] = {h: len(v) for h, v in host_groups.items()}

        # Compute all pairwise p-distances
        all_samples = []
        for ids in host_groups.values():
            all_samples.extend(ids)
        n = len(all_samples)
        dist_matrix = np.zeros((n, n))
        for i in range(n):
            si = seqs[all_samples[i]]
            for j in range(i + 1, n):
                # 2026-09-15 修复: 旧版 `diffs = sum(a != b ...)` 在 MSA 上把 gap/N
                # 当真实差异 → 宿主 within/between 距离与置换 p 值系统性偏高。
                d = p_distance_pair(si, seqs[all_samples[j]])
                dist_matrix[i, j] = d
                dist_matrix[j, i] = d

        # Within vs Between
        within_dists = []
        between_dists = []
        host_of = {}
        idx = 0
        for host, ids in host_groups.items():
            for sid in ids:
                host_of[idx] = host
                idx += 1

        for i in range(n):
            for j in range(i + 1, n):
                if host_of[i] == host_of[j]:
                    within_dists.append(dist_matrix[i, j])
                else:
                    between_dists.append(dist_matrix[i, j])

        within_mean = np.mean(within_dists) if within_dists else 0
        between_mean = np.mean(between_dists) if between_dists else 0
        obs_diff = between_mean - within_mean

        result["within_mean"] = float(within_mean)
        result["between_mean"] = float(between_mean)

        # Permutation test
        rng = np.random.RandomState(42)
        perm_diffs = []
        host_labels = np.array([host_of[i] for i in range(n)])
        for _ in range(n_permutations):
            shuffled = rng.permutation(host_labels)
            w_dists = []
            b_dists = []
            for i in range(n):
                for j in range(i + 1, n):
                    if shuffled[i] == shuffled[j]:
                        w_dists.append(dist_matrix[i, j])
                    else:
                        b_dists.append(dist_matrix[i, j])
            perm_diff = np.mean(b_dists) - np.mean(w_dists)
            perm_diffs.append(perm_diff)
        perm_diffs = np.array(perm_diffs)

        p_value = (np.sum(perm_diffs >= obs_diff) + 1) / (n_permutations + 1)
        result["p_value"] = float(p_value)

        effect_size = obs_diff / np.std(np.concatenate([within_dists, between_dists])) if within_dists and between_dists else 0
        log.emit(f"Within: {within_mean:.5f}, Between: {between_mean:.5f}, diff={obs_diff:.5f}, p={p_value:.4f}, effect={effect_size:.3f}")

        # Plot
        # 2026-10-12: 配色统一 OKABE（原 #2ecc71/#e74c3c/#3498db 均不在项目调色板内）；
        # dpi 150 → 300；补 y 轴网格。画布 14 in 保留 —— 排进 Fig 3 的通栏
        # （9.32 in）缩放 0.67 倍，13 pt 标题落到 8.7 pt，比例本身合适。
        apply_sci_style(base_size=10)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Boxplot
        box_data = [within_dists, between_dists]
        bp = ax1.boxplot(box_data, patch_artist=True, widths=0.5)
        bp['boxes'][0].set_facecolor(OKABE_ITO[4])      # blue   — within host
        bp['boxes'][1].set_facecolor(OKABE_ITO[5])      # vermillion — between hosts
        for b in bp['boxes']:
            b.set_edgecolor("#333333")
            b.set_linewidth(0.8)
        for el in ('whiskers', 'caps', 'medians'):
            for art in bp[el]:
                art.set_color("#333333")
                art.set_linewidth(0.8)
        ax1.set_xticks([1, 2])
        # 2026-10-12（审稿人视角）：箱线图必须给出每组样本量 —— 否则读者无法判断
        # "between − within = 0.007" 是建立在多少对比较上的。
        ax1.set_xticklabels([f'Within host\n({len(within_dists):,} pairs)',
                             f'Between hosts\n({len(between_dists):,} pairs)'],
                            fontsize=12)
        ax1.set_ylabel('p-distance', fontsize=13)
        ax1.set_title(f'Host differentiation\nbetween − within = {obs_diff:.5f}, '
                      f'$P$ = {p_value:.4f}',
                      fontsize=13, fontweight='bold')
        ax1.grid(alpha=0.35, axis='y')
        despine(ax1)

        # Permutation histogram
        ax2.hist(perm_diffs, bins=50, color=OKABE_ITO[4], alpha=0.75,
                 edgecolor='white', linewidth=0.3)
        ax2.axvline(obs_diff, color=OKABE_ITO[5], linewidth=2.0, linestyle='--',
                    label=f'Observed difference = {obs_diff:.5f}')
        ax2.set_xlabel('Between − within mean distance', fontsize=13)
        ax2.set_ylabel('Frequency', fontsize=13)
        # 2026-10-12: 把 P 值同时写进右图标题（右图才是 P 的来源），
        # 并用千分位标出置换次数。
        ax2.set_title(f'Permutation test  ({n_permutations:,} permutations)\n'
                      f'$P$ = {p_value:.4f}',
                      fontsize=13, fontweight='bold')
        ax2.legend(fontsize=10, loc="upper left")
        ax2.grid(alpha=0.35, axis='y')
        despine(ax2)

        plt.tight_layout()
        plot_path = os.path.join(output_dir, "host_differentiation.pdf")
        save_both(fig, plot_path, bbox_inches='tight', dpi=300, format='pdf')
        plt.close(fig)
        result["plot"] = plot_path
        result["success"] = True
        log.emit(f"Host differentiation: {plot_path}")
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"Host differentiation failed: {e}")
        return result


# ═══════════════════════════════════════════════════════════════════
# Module 11: NCBI SeqHarvester + SeqIDRenamer + SeqGrouper (桥接 VirPhyKit 源码)
# ═══════════════════════════════════════════════════════════════════
# 以下函数通过专用 bridge 文件直接复用 VirPhyKit 源码逻辑
# 而非重新实现


def seqid_renamer(fasta_file, rename_map, output_dir):
    """序列 ID 重命名 (自研, 等效 VirPhyKit SeqIDRenamer)。"""
    import tempfile
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False)
    for old, new in rename_map.items(): tmp.write(f"{old}\t{new}\n")
    tmp.close()
    try:
        from utils.seqid_renamer import rename_sequences
        return rename_sequences(fasta_file, tmp.name, output_dir)
    finally: os.unlink(tmp.name)


# ═══════════════════════════════════════════════════════════════════
# Module RRT: Region Randomization Test (桥接 VirPhyKit 源码)
# ═══════════════════════════════════════════════════════════════════

def run_rrt(
    original_tree: str,
    randomized_trees: List[str],
    output_dir: str,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    桥接 VirPhyKit RRT: 区域随机化检验。

    从 BEAST MCC 树解析 set={} 和 set.prob={} 注释，
    比较真实数据 vs 区域随机化数据的根状态后验概率。

    复用 VirPhyKit src/RRT/function_rrt.py 的 generate_table + plot_graph_from_csv 逻辑。
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.backends.backend_pdf

    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "passed": None, "csv_table": None, "plot": None, "error": None}

    try:
        # Parse original tree annotations (VirPhyKit browse_original_data_file logic)
        def parse_tree_annotations(filepath):
            with open(filepath, encoding='utf-8', errors='replace') as f:
                content = f.read()
            # 2026-09-15 (审查 P2-13): 旧版正则 \b\w+\.set=\{...\} 强制要求前缀
            # (如 "Region.set={...}"), 裸 "set={...}" 直接漏掉 → RRT 误报"无注释"。
            # 与同文件 _parse_tree_annotations (rspp 路径) 统一为可选前缀。
            set_matches = re.findall(r'(?:\w+\.)?set=\{\s*([^}]*)\s*\}', content)
            prob_matches = re.findall(r'(?:\w+\.)?set\.prob=\{\s*([^}]*)\s*\}', content)
            last_set = set_matches[-1] if set_matches else None
            last_prob = prob_matches[-1] if prob_matches else None
            if last_set:
                last_set = ','.join([c.strip().strip('"') for c in last_set.split(',')])
            if last_prob:
                prob_vals = [float(p.strip()) for p in last_prob.split(',')]
                if all(abs(prob_vals[i]-prob_vals[i+1])<0.01 for i in range(len(prob_vals)-1)):
                    return {"set": last_set, "set_prob": last_prob, "skip": True}
            return {"set": last_set, "set_prob": last_prob, "skip": False}

        orig = parse_tree_annotations(original_tree)
        if not orig["set"]:
            result["error"] = "No set={} annotation found in original tree"; return result
        if orig["skip"]:
            # 上游 :38-42 的 heuristic：相邻后验差 < 0.01 → "too small and no detection
            # is needed"。这里**保持与上游一致**（passed=True），但额外给出
            # `skipped=True` 供调用方区分"无需检验"与"检验通过" ——
            # `_review20260916/AG-C_geo.md:43` 指出：若把两者混为一谈，
            # 全平后验会被当成"地理信号显著"的假阳性。
            result["passed"] = True
            result["skipped"] = True
            result["note"] = ("original-tree posterior is near uniform (adjacent difference < 0.01) → upstream convention calls this "
                              "'test not needed', not 'test passed'")
            result["error"] = "Posterior probability too small, RRT not needed"
            return result

        # Parse randomized trees
        rand_data = []
        for rt in randomized_trees:
            d = parse_tree_annotations(rt)
            if d["set"]:
                rand_data.append(d)
        log.emit(f"RRT: 1 original + {len(rand_data)} randomized trees")

        # generate_table (VirPhyKit logic)
        all_data = []
        for i, d in enumerate(rand_data, 1):
            row = {"Replicates": f"Random{i}"}
            countries = d["set"].split(',')
            probs = d["set_prob"].split(',')
            for c, p in zip(countries, probs):
                row[c.strip()] = float(p.strip())
            all_data.append(row)

        df = pd.DataFrame(all_data)
        real_row = {"Replicates": "Real"}
        countries = orig["set"].split(',')
        probs = [float(p.strip()) for p in orig["set_prob"].split(',')]
        for c, p in zip(countries, probs):
            real_row[c.strip()] = p

        target_country = countries[probs.index(max(probs))].strip()
        max_prob = max(probs)

        # ── Min/Max 行：只在**随机行**上取值（口径 = 上游 `df.iloc[:-2, 1:]`）──
        rand_numeric = df.select_dtypes(include='number')
        min_vals = rand_numeric.min() if not rand_numeric.empty else pd.Series(dtype=float)
        max_vals = rand_numeric.max() if not rand_numeric.empty else pd.Series(dtype=float)

        csv_path = os.path.join(output_dir, "rrt_results.csv")
        # VirPhyKit `generate_table`（上游 :79-89）的 CSV 契约：
        #   行序 = Random1..N → 空行 → Real → Min → Max；
        #   行标签**一律写进 `Replicates` 列**；
        #   上游 `plot_graph_from_csv`(:103-106) 是**按行标签**找 Real/Min/Max 的，
        #   标签放错列 → 直接 `ValueError("Required rows (Real, Min, Max) not found")`。
        # 2026-09-16 实测修复：此前 Min/Max 用 `{'node': 'Min'}` —— 键名写错，生成的 CSV
        # Replicates 列里只有 Real，Min/Max 落到多出来的 `node` 列 → 与上游
        # plot_graph_from_csv 的互操作**实际是断的**（原注释声称"同款契约"，未达成）。
        min_row = {'Replicates': 'Min'}
        max_row = {'Replicates': 'Max'}
        for col in rand_numeric.columns:
            min_row[col] = rand_numeric[col].min()
            max_row[col] = rand_numeric[col].max()
        df_out = pd.concat(
            [df,                                    # Random1..N
             pd.DataFrame([{}]),                    # 空行分隔（上游 :80-81）
             pd.DataFrame([real_row]),              # Real
             pd.DataFrame([min_row]),               # Min
             pd.DataFrame([max_row])],              # Max
            ignore_index=True)
        # utf-8-sig：上游 :89 就是 utf-8-sig；缺 BOM 时 Excel 打开中文区划名会乱码
        df_out.to_csv(csv_path, index=False, encoding='utf-8-sig')

        # Check pass/fail (VirPhyKit logic)
        rand_target = [r.get(target_country, 0.0) for r in all_data]
        max_rand = max(rand_target) if rand_target else 0.0
        result["passed"] = max_prob > max_rand

        # Plot
        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['font.family'] = 'Arial'
        fig, ax = plt.subplots(figsize=(10, 6))
        # x 轴用**合并后**的列（上游 :110 `df.columns[1:]` 取的是含 Real 的合并表）：
        # 只在 Real 里出现的区划也要画出来。缺值的用 NaN（matplotlib 会断线），
        # 用 0 会画出一条并不存在的"跌到 0"。
        cols = [c for c in df_out.columns if c != 'Replicates']
        ax.plot(cols, [real_row.get(c, float('nan')) for c in cols], 'o-',
                label='Real', color='#5BD4D8')
        if not min_vals.empty:
            ax.plot(cols, [min_vals.get(c, float('nan')) for c in cols], 'o--',
                    label='Min', color='#FFDBC6')
            ax.plot(cols, [max_vals.get(c, float('nan')) for c in cols], 'o--',
                    label='Max', color='#FE692E')
        ax.set_title(f"Region Randomization Test ({'PASSED' if result['passed'] else 'FAILED'})")
        ax.set_xlabel("Region"); ax.set_ylabel("Posterior Probability")
        ax.legend(); ax.grid(True, alpha=0.3)
        plt.xticks(rotation=45, ha='right')

        plot_path = os.path.join(output_dir, "rrt_plot.pdf")
        save_both(fig, plot_path, bbox_inches='tight', dpi=150, format='pdf')
        plt.close(fig)

        result["csv_table"] = csv_path
        result["plot"] = plot_path
        result["success"] = True
        log.emit(f"RRT: {'PASSED' if result['passed'] else 'FAILED'} → {plot_path}")
        return result

    except Exception as e:
        result["error"] = str(e); return result


# ═══════════════════════════════════════════════════════════════════
# Module TempMig: Temporal Migration Tracker（自研: tempmig_full MOT 复刻）
# ═══════════════════════════════════════════════════════════════════

def run_tempmig(
    mcc_tree: str,
    output_dir: str,
    trait: str = "type",
    dates_csv: str = None,
    perl_path: str = "perl",
    rscript_path: str = "Rscript",
    python_path: str = "python3",
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    从 BEAST MCC 树重建时空迁移矩阵（自研 tempmig_full，等效 VirPhyKit TempMig/MOT）。

    输出: migration_over_time.csv (迁移矩阵×年份) + migration_over_time.pdf (时序迁移图)。

    Parameters
    ----------
    mcc_tree : str — BEAST MCC 树 (含 trait 注释, 叶子名需与 dates_csv 匹配)
    output_dir : str
    trait : str — trait 名称 (默认 type)
    dates_csv : str — 采样时间表 (必需; 首列 tip 名, 次列日期)
    perl_path, rscript_path, python_path : str — 保留参数，仅兼容旧调用（不再使用）
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "migration_matrix": None, "plot": None, "error": None}

    if not dates_csv or not os.path.exists(dates_csv):
        result["error"] = ("自研 TempMig (MOT 复刻) 需要 dates_csv 采样时间表: "
                           "首列 tip 名, 次列日期")
        log.error(result["error"])
        return result

    try:
        from utils.tempmig_full import run_tempmig_full
        r = run_tempmig_full(mcc_tree, dates_csv, output_dir, trait=trait)
        if r.get('success'):
            result['success'] = True
            result['migration_matrix'] = r.get('matrix')
            result['plot'] = r.get('plot')
            result['total_migrations'] = r.get('total_migrations')
            result['top_routes'] = r.get('top_routes')
            log.emit(f"TempMig (自研 MOT): 迁移矩阵 → {r.get('matrix')}")
        else:
            result['error'] = r.get('error') or 'tempmig_full 失败'
            log.error(result['error'])
        return result
    except Exception as e:
        result['error'] = str(e)
        log.error(f"TempMig failed: {e}")
        return result


# ═══════════════════════════════════════════════════════════════════
# Module: Gene-Level Analysis (逐基因时钟性比较)
# ═══════════════════════════════════════════════════════════════════

def slice_genes_from_alignment(
    gbk_path: str,
    aln_path: str,
    out_dir: str,
    log: Optional[LogCollector] = None,
) -> List[Dict]:
    """从 GenBank CDS 坐标切分全基因组 MSA，产出各基因比对。"""
    log = log or LogCollector()
    os.makedirs(out_dir, exist_ok=True)
    results = []
    try:
        cds_list = []
        ref_id = None
        for rec in SeqIO.parse(gbk_path, 'genbank'):
            ref_id = rec.id
            log.emit(f"Genome: {len(rec.seq)} bp, {rec.id}")
            for _fi, feat in enumerate(rec.features):
                if feat.type == 'CDS':
                    # 2026-09-15 修复: 旧版直接取 /gene, 缺该限定符时为 '?' →
                    # 下游拼出 '?.mafft.fasta' (Windows 实测 OSError [Errno 22]),
                    # 且含 / \ : 等字符也会写文件失败。现做多级回退 + 文件名清洗。
                    gene = str(feat.qualifiers.get('gene', [''])[0]).strip()
                    if not gene or gene in ('?', '.', '-'):
                        gene = str(feat.qualifiers.get('locus_tag', [''])[0]).strip()
                    if not gene or gene in ('?', '.', '-'):
                        gene = str(feat.qualifiers.get('protein_id', [''])[0]).strip()
                    if not gene or gene in ('?', '.', '-'):
                        gene = f"CDS{_fi}"
                    gene = re.sub(r'[^A-Za-z0-9_.\-]', '_', gene) or f"CDS{_fi}"
                    product = feat.qualifiers.get('product', ['?'])[0]
                    start = int(feat.location.start)
                    end = int(feat.location.end)
                    cds_list.append({'gene': gene, 'start': start, 'end': end,
                                   'len': end - start, 'product': product})
        if not cds_list:
            log.error("No CDS features found")
            return results
        aln_records = list(SeqIO.parse(aln_path, 'fasta'))
        # 2026-09-15 修复: 旧版参考未命中时 ref_idx 静默留 0 → 拿比对第一条当参考,
        # CDS 坐标整体错位且无任何提示。另 `ref_id in r.id` 是子串匹配, 会误命中
        # (如 NC_003636 命中 NC_0036361)。现改为 精确 → 前缀 两级, 命中失败即报错。
        ref_idx = None
        for i, r in enumerate(aln_records):
            if r.id == ref_id:
                ref_idx = i
                break
        if ref_idx is None:
            for i, r in enumerate(aln_records):
                if r.id.startswith(ref_id) or ref_id.startswith(r.id):
                    ref_idx = i
                    break
        if ref_idx is None:
            log.error(f"参考序列 {ref_id} 不在比对中 (比对共 {len(aln_records)} 条); "
                      f"无法建立 CDS→比对列 坐标映射, 跳过基因切分")
            return results
        ref_aln = str(aln_records[ref_idx].seq).upper()
        aln_to_genome = {}
        gpos = 0
        for apos, base in enumerate(ref_aln):
            if base != '-':
                aln_to_genome[apos] = gpos; gpos += 1
        for cds in cds_list:
            gene = cds['gene']
            aln_start = aln_end = None
            for apos, gv in aln_to_genome.items():
                if cds['start'] <= gv < cds['end']:
                    if aln_start is None: aln_start = apos
                    aln_end = apos
            if aln_start is None: continue
            gene_recs = [SeqRecord(seq=Seq(str(r.seq)[aln_start:aln_end+1]), id=r.id, description='')
                        for r in aln_records]
            out_fa = os.path.join(out_dir, f"{gene}.mafft.fasta")
            SeqIO.write(gene_recs, out_fa, 'fasta')
            results.append({'gene': gene, 'genome_start': cds['start'], 'genome_end': cds['end'],
                          'aln_length': aln_end - aln_start + 1, 'fasta': out_fa, 'product': cds['product']})
        return results
    except Exception as e:
        log.error(f"Gene slicing failed: {e}")
        return results


def run_gene_rtt_batch(
    gene_alignments: List[Dict], metadata_csv: str, out_dir: str,
    iqtree_bin: str = "iqtree2", threads: int = 8,
    log: Optional[LogCollector] = None,
) -> List[Dict]:
    """对每个基因运行 IQ-TREE + TreeTime-RTT，返回汇总表。

    2026-08-31 加速: 基因间并行 (ThreadPoolExecutor)。
    小比对 (几百 bp) 单基因吃不满全部线程, 基因间并行收益远大于单基因加线程;
    每基因 IQ-TREE 线程 = max(1, threads // n_genes)。
    """
    log = log or LogCollector()
    n_genes = max(1, len(gene_alignments))
    per_gene_thr = max(1, threads // n_genes)

    def _run_one(ga):
        gene = ga['gene']
        gene_out = os.path.join(out_dir, gene)
        os.makedirs(gene_out, exist_ok=True)
        tree_prefix = os.path.join(gene_out, gene)
        tree_file = f"{tree_prefix}.treefile"
        if not os.path.exists(tree_file):
            # 2026-09-15 (P1-12): shell=True 拼接的路径必须加引号 ——
            # out-dir 含空格 (如 Windows `D:\桌面\...`) 会被拆成多个参数;
            # 且 gene 名来自 GenBank /gene 限定符, 含 ; & 等即构成命令注入面。
            # 2026-09-16 说明: 这里**刻意不加** -B/-alrt —— 基因树的唯一消费者是
            # TreeTime-RTT (只要拓扑+枝长, 不用支持值), 加自举纯属浪费算力;
            # 论文里报告支持值的是**主树** (phylo_pipeline.run_iqtree, 已补
            # `-B 1000 -alrt 1000`)。此处的 `-m MFP --quiet` 不是遗漏。
            import shlex as _sh
            cmd = (f'{_sh.quote(str(iqtree_bin))} -s {_sh.quote(str(ga["fasta"]))} '
                   f'-pre {_sh.quote(str(tree_prefix))} '
                   f'-nt {int(per_gene_thr)} -m MFP --quiet')
            try:
                r = subprocess.run(cmd, shell=True, capture_output=True, text=True, encoding='utf-8', errors='replace',
                                   timeout=600, cwd=gene_out)
            except Exception:
                r = None
            if r is None or r.returncode != 0:
                # 历史坑: 失败条目缺 r_squared/beta 键 → 汇总循环 r["r_squared"] KeyError
                # 且异常上传 run_stage_genes 无包裹 → 一个基因失败崩整个 genes stage
                return {'gene': gene, 'product': ga.get('product', ''),
                        'length': ga.get('aln_length'), 'success': False,
                        'error': 'IQ-TREE failed', 'beta': None,
                        'r_squared': None, 'p_value': None}
        rtt = run_treetime_rtt(fasta_file=ga['fasta'], tree_file=tree_file,
                              dates_file=metadata_csv, output_dir=gene_out,
                              log=None)  # 并行时不逐基因刷日志
        return {'gene': gene, 'product': ga.get('product', ''),
                'length': ga.get('aln_length'), 'success': rtt['success'],
                'beta': rtt.get('beta'), 'r_squared': rtt.get('r_squared'),
                'p_value': rtt.get('p_value')}

    log.emit(f"  基因级并行: {n_genes} 基因 × {per_gene_thr} 线程/基因")
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=n_genes) as ex:
        results = list(ex.map(_run_one, gene_alignments))
    for ga, r in zip(gene_alignments, results):
        log.emit(f"\n--- Gene: {ga['gene']} ({ga.get('product','')}) ---")
        log.emit(f"  R²={r.get('r_squared')}, β={r.get('beta')}")
    log.emit(f"\n{'Gene':6s}  {'R²':>8s}  {'β':>12s}" + "\n" + "-"*30)
    for r in results:
        r2 = f'{r["r_squared"]:.4f}' if r['r_squared'] else 'N/A'
        beta = f'{r["beta"]:.3e}' if r['beta'] else 'N/A'
        log.emit(f'{r["gene"]:6s}  {r2:>8s}  {beta:>12s}')
    return results


# ═══════════════════════════════════════════════════════════════════
# Utility
# ═══════════════════════════════════════════════════════════════════

def check_dependencies() -> Dict[str, bool]:
    """检查 virome_phylo_pipeline 运行时依赖 (无需 VirPhyKit)。

    注：R 相关的 TreeDater/MOT 等功能需要系统 Rscript + 相应 R 包。
    """
    deps = {
        "treetime": HAS_TREETIME,
        "biopython": True,
        "matplotlib": True,
        "pandas": True,
        "scipy": True,
        "geopy": False,
        "rscript": False,
    }
    try:
        import geopy
        deps["geopy"] = True
    except ImportError:
        pass
    try:
        subprocess.run(["Rscript", "--version"], capture_output=True, timeout=5, check=True)
        deps["rscript"] = True
    except Exception:
        pass
    return deps
