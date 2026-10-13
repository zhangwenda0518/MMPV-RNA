#!/usr/bin/env python3
"""
treetime_bridge.py — TreeTime CLI 完整时间树桥 (timetree.nexus)
================================================================

phymapr TreeTime 路线的第一环 (2026-10-12 新增, 参考 ADHS-Taylor/phymap-workflow)。

与 utils/virphy_bridge.py::run_treetime_rtt 的分工:
  · run_treetime_rtt —— TreeTime **Python API**, 硬编码 Jukes-Cantor/infer_gtr=False,
    只产出无注解 Newick (timetree_inferred.nwk), 服务于 clock 阶段的 RTT 回归口径;
  · run_treetime_full (本模块) —— 调 **`python -m treetime` CLI**, 打开 clock-filter
    (IQD 离群 tip 剔除)、coalescent(skyline) 先验、--confidence 置信区间与固定 RNG 种子,
    产出带 `[&date=...]` 注解的 **timetree.nexus** —— 这是 phymapr
    `generate_phylo_transmission()` 的硬性输入 (treeio::read.beast 要求 date 注解),
    全仓库此前没有任何路径能产出该格式。

日期口径说明: 本模块只做"有没有可解析日期"的过滤 (委托
utils/metadata_governance.decimal_year 作有效性判定, 与 TreeDater 同源),
**不改写日期值** —— 树的采样时间标度由 TreeTime CLI 自己 parse_dates
(YYYY / YYYY-MM / YYYY-MM-DD / 小数年通吃), 与 run_treetime_rtt 同口径;
仅在无日期 tip 存在时按有日期 tip 集剪树 (照抄 TreeDater 的 prune 模式)。
"""

import os
import re
import subprocess
import sys
from typing import Dict, Optional

from utils.virphy_bridge import LogCollector

# TreeTime CLI 写出的时间树文件名, 按命中顺序回退 (phymap-workflow 同款清单)
_TIMETREE_CANDIDATES = (
    "timetree.nexus", "timetree.nex", "timetree.nwk", "annotated_tree.nexus",
)


def _dated_tip_set(dates_csv: str, log: LogCollector):
    """读 dates.csv → ({有可解析日期的 tip 名集合}, 全部 tip 名集合)。

    有效性判定委托 metadata_governance.decimal_year (与 TreeDater 同一实现):
    占位符/纯文本日期返回 None → 该 tip 视为无日期。
    """
    import pandas as pd

    df = pd.read_csv(dates_csv)
    if df.shape[1] < 2:
        raise ValueError(f"dates 文件至少需要 name/date 两列: {dates_csv}")
    name_col, date_col = df.columns[0], df.columns[1]

    from utils.metadata_governance import decimal_year as _gov_decimal_year

    def _is_valid(d):
        try:
            return _gov_decimal_year(d, mode="start") is not None
        except (ValueError, IndexError, TypeError):
            return False

    df["_ok"] = df[date_col].apply(_is_valid)
    dated = set(df.loc[df["_ok"], name_col].astype(str))
    all_tips = set(df[name_col].astype(str))
    log.emit(f"dates: {len(dated)}/{len(all_tips)} tips 有可解析日期")
    return dated, all_tips


def _prune_to_dated(tree_file: str, dated: set, output_dir: str, log: LogCollector):
    """把树剪到只剩有日期 tip (照抄 run_treedater_ltt 的 Bio.Phylo prune 模式)。

    返回 (剪后树路径, n_before, n_after); 剪枝失败返回 (原路径, 0, 0) 并只告警 ——
    与 TreeDater 一致: 树剪不动时降级用原树, 让上游工具自己报错。
    """
    try:
        from Bio import Phylo

        tree = Phylo.read(tree_file, "newick")
        tips = list(tree.get_terminals())
        n_before = len(tips)
        undated = [c for c in tips if c.name not in dated]
        if not undated:
            return tree_file, n_before, n_before
        for tip in undated:
            tree.prune(tip)
        pruned = os.path.join(output_dir, "timetree_input_dated.nwk")
        Phylo.write(tree, pruned, "newick")
        n_after = len(list(tree.get_terminals()))
        log.emit(f"树已按有日期 tip 剪枝: {n_before} → {n_after}")
        return pruned, n_before, n_after
    except Exception as e:  # noqa: BLE001
        log.warning(f"树剪枝失败 (降级用原树, TreeTime 可能因缺日期报错): {e}")
        return tree_file, 0, 0


def run_treetime_full(
    tree_file: str,
    aln_file: str,
    dates_csv: str,
    output_dir: str,
    clock_filter: float = 4.0,
    coalescent: str = "skyline",
    clock_rate: Optional[float] = None,
    rng_seed: int = 42,
    timeout: int = 3600,
    log: Optional[LogCollector] = None,
) -> Dict:
    """跑 `python -m treetime` CLI 产出带 date 注解的 timetree.nexus。

    Parameters
    ----------
    tree_file : str        ML 树 (IQ-TREE treefile, newick)
    aln_file : str         MAFFT 全长比对 (等长, gap 允许)
    dates_csv : str        管线时间通道表 (列 name,date[,location], 原始日期串)
    output_dir : str       产物目录 (TreeTime CLI 直接写入)
    clock_filter : float   IQD 离群 tip 剔除阈值 (0=关; TreeTime CLI 默认 4.0)
    coalescent : str       coalescent 先验 (const / skyline / optimal / 0.0)
    clock_rate : float     手动钟速率 (慢进化病原用; None=让 TreeTime 回归估计)
    rng_seed : int         固定种子 (管线惯例: 可复现)
    timeout : int          子进程超时秒数

    Returns
    -------
    dict: {success, timetree_nexus, metrics{clock_rate,r_squared,n_tips_total,
          n_tips_pruned,...}, output_dir, error}
          永不抛异常, 失败一律 {"success": False, "error": ...}。
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result: Dict = {"success": False, "timetree_nexus": None,
                    "output_dir": output_dir, "metrics": {}, "error": None}

    try:
        import treetime  # noqa: F401  探测依赖 (CLI 由同一解释器拉起)
    except ImportError:
        result["error"] = ("treetime 未安装 (pip install phylo-treetime)。"
                           "注: clock 阶段的 run_treetime_rtt 同样依赖它, "
                           "若 clock 正常则本探测不应失败")
        log.error(result["error"])
        return result

    for label, path in (("树", tree_file), ("比对", aln_file), ("日期", dates_csv)):
        if not path or not os.path.exists(path):
            result["error"] = f"输入缺失: {label} {path}"
            log.error(result["error"])
            return result

    try:
        dated, _all = _dated_tip_set(dates_csv, log)
    except Exception as e:  # noqa: BLE001
        result["error"] = f"日期表读取失败: {e}"
        log.error(result["error"])
        return result
    if len(dated) < 3:
        result["error"] = f"仅 {len(dated)} 个 tip 有可解析日期 (需 ≥3)"
        log.error(result["error"])
        return result

    # 无日期 tip → 剪树 + 写只含 (name,date) 两列的过滤日期表
    tree_use, n_total, n_after = _prune_to_dated(tree_file, dated, output_dir, log)
    n_pruned = max(0, n_total - n_after) if n_total else 0

    import pandas as pd
    df = pd.read_csv(dates_csv)
    name_col, date_col = df.columns[0], df.columns[1]
    df_f = df[df[name_col].astype(str).isin(dated)][[name_col, date_col]]
    dates_tsv = os.path.join(output_dir, "dates_filtered.tsv")
    df_f.to_csv(dates_tsv, sep="\t", index=False, header=["name", "date"])

    cmd = [
        sys.executable, "-m", "treetime",
        "--tree", tree_use,
        "--aln", aln_file,
        "--dates", dates_tsv,
        "--outdir", output_dir,
        "--clock-filter", str(clock_filter),
        "--coalescent", str(coalescent),
        "--confidence",
        "--rng-seed", str(rng_seed),
        "--verbose", "1",
    ]
    if clock_rate:
        cmd += ["--clock-rate", str(clock_rate)]

    log.emit(f"Running: {' '.join(cmd)}")
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              encoding="utf-8", errors="replace", timeout=timeout)
    except subprocess.TimeoutExpired:
        result["error"] = f"TreeTime CLI 超时 (>{timeout}s)"
        log.error(result["error"])
        return result
    except Exception as e:  # noqa: BLE001
        result["error"] = f"TreeTime CLI 启动失败: {e}"
        log.error(result["error"])
        return result

    if proc.stdout:
        (log.emit if log.emit else print)(proc.stdout.strip()[-1500:])

    timetree = None
    for name in _TIMETREE_CANDIDATES:
        cand = os.path.join(output_dir, name)
        if os.path.exists(cand) and os.path.getsize(cand) > 0:
            timetree = cand
            break

    if proc.returncode != 0 or timetree is None:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-15:]
        result["error"] = ("TreeTime CLI 失败"
                           + (f" (exit {proc.returncode})" if proc.returncode else " (无 timetree 产物)")
                           + ": " + " | ".join(tail))
        log.error(result["error"])
        return result

    m = result["metrics"]
    m["clock_filter_iqd"] = clock_filter
    m["coalescent"] = coalescent
    m["rng_seed"] = rng_seed
    m["n_tips_total"] = n_total or len(dated)
    m["n_tips_used"] = len(dated)
    if n_pruned:
        m["n_tips_pruned_undated"] = n_pruned
    if clock_rate:
        m["clock_rate_manual"] = clock_rate

    # molecular_clock 报告 (TreeTime CLI 产物; 键名随版本浮动, 宽松解析)
    for mc_name in ("molecular_clock.txt", "molecular_clock.tsv", "molecular_clock.csv"):
        mc = os.path.join(output_dir, mc_name)
        if os.path.exists(mc):
            m["molecular_clock_report"] = mc
            try:
                with open(mc, encoding="utf-8", errors="replace") as f:
                    txt = f.read()
                mr = re.search(r"rate[^\d\n]*([\d.eE+\-]+)", txt, re.I)
                if mr:
                    m["clock_rate"] = float(mr.group(1))
            except Exception:  # noqa: BLE001
                pass
            break
    # stdout 兜底抓速率/R² (TreeTime 打印 'Inferred association rate: x (R^2 = y)' 类行)
    blob = f"{proc.stdout or ''}\n{proc.stderr or ''}"
    if "clock_rate" not in m:
        mr = re.search(r"[Rr]ate[^\d\n]*?([\d.]+e?[-+]?\d*)\s*\(?.*?R\^?2?\s*=\s*([\d.]+)", blob)
        if mr:
            try:
                m["clock_rate"] = float(mr.group(1))
                m["r_squared"] = float(mr.group(2))
            except ValueError:
                pass

    result["success"] = True
    result["timetree_nexus"] = timetree
    log.emit(f"TreeTime 时间树完成: {timetree}")
    return result
