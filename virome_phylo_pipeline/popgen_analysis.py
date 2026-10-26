#!/usr/bin/env python3
"""popgen_analysis.py — 病毒群体遗传学分析模块 v2 (纯 pypopart, 无 R 依赖)

整合群体遗传学框架, 工具链统一 pypopart:
  1. 遗传多样性: π, θ-W, S, Hd (pypopart)
  2. 中性检验: Tajima's D + Fu's Fs (pypopart)
  3. 遗传分化: Fst (pypopart 单倍型频率法, Nei) + 位点法补充
  4. 保守性图谱: 逐位点 Shannon 熵

用法:
  from popgen_analysis import PopGenAnalyzer
  pg = PopGenAnalyzer(fasta, metadata_csv)
  d = pg.diversity()          # π/S/θ/Hd/Tajima's D/Fu's Fs
  f = pg.fst('location')      # Fst 矩阵 + 总 Fst
  pg.amova('location')        # AMOVA
  pg.entropy_profile()        # Shannon 熵谱

CLI:
  python -m popgen_analysis --fasta X --metadata Y --group location
"""
import argparse
import json
import math
import sys
import os
import numpy as np
import csv
from collections import Counter, defaultdict
from pathlib import Path

# pypopart editable 安装后从 site-packages 的 .pth 找到; 显式路径仅兼容旧布局
# 历史坑: 写死 /tmp/pypopart/src (已被 /tmp 清理), 导致 ModuleNotFoundError
# 2026-09-15 (审查 P2-22): 兜底路径改为**相对仓库根目录**推导, 不再硬编码
# /home/zhangwenda/... (换机器即静默失效)。
from utils.dataset_config import env as get_env
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PYPOPART_PATH = (get_env('pypopart_path')
                 or os.path.join(_REPO_ROOT, 'biosoft', 'pypopart', 'src'))
if os.path.isdir(PYPOPART_PATH) and PYPOPART_PATH not in sys.path:
    sys.path.insert(0, PYPOPART_PATH)

from pypopart.io import load_alignment
from pypopart.core.alignment import Alignment
from pypopart.core.graph import HaplotypeNetwork
from pypopart.stats import popgen as pg_stats

from utils.seq_ids import base_sample_id


# 参考/外群识别模式 (小写子串匹配)。
# 2026-09-15 (审查 P2-15): 比对里常把参考基因组 (或外群) 一并放进来做坐标锚定。
# 参考序列参与 π/S/θ-W/Hd/Shannon 熵会把多样性系统性抬高 (它和所有样本都不同),
# 而旧版没有任何开关能把它排除, 也没有任何提示。
REFERENCE_ID_PATTERNS = (
    'reference', 'ref_', '_ref', 'outgroup', '外群', '参考', 'nc_', 'ref|',
)


def looks_like_reference(seq_id):
    """按命名约定判断是否像参考序列/外群。"""
    s = str(seq_id).lower()
    return any(p in s for p in REFERENCE_ID_PATTERNS)


# ═══════════════════════════════════════════════════════════════════
# 核心统计 (模块级, 供 全局 diversity() 与 按群体 diversity_by_population 共用;
# 2026-10-12 新增: 配套 popgen_stats.json / popgen_plots.py 可视化)
# ═══════════════════════════════════════════════════════════════════

def _fu_fs(n, k):
    """Fu (1997) Fs — 经典 Ewens 口径。

    θ̂ 由观测单倍型数 k 拟合 Ewens 期望 E[K]=Σ_{j=0}^{n-1} θ/(θ+j);
    Fs = ln Σ_{i≤k} P(i | θ̂, n), P 由 Ewens 递推
    P(k|j) = (j-1)/(θ+j-1)·P(k|j-1) + θ/(θ+j-1)·P(k-1|j-1)。
    负值 = 单倍型数超出中性期望 = 群体扩张信号。
    注: Fu 氏 Fs 属**单倍型数检验族** (Fu 1997), Arlequin 口径下 θ̂ 由观测
    单倍型数 k 经 MLE 拟合 Ewens 期望 E[K] 得到, Fs = ln P(K≤k | θ̂, n)
    —— 故 Fs 只依赖 (n, k), 这是该检验族的固有性质而非缺陷。
    pypopart 自带的 calculate_fu_fs 是 ±ln(k/E[k]) 的"simplified approximation"
    (其源码自认, 且 if/else 两分支恒等), 本实现用完整 Ewens 递推, 与
    Arlequin/DnaSP 经典定义同族; 中性预期在 ln(0.5)≈-0.69 附近小幅波动,
    扩张群体显著更负。引用时注明口径。全唯一单倍型 (k≥n, θ̂→∞) 时 Fs
    未定义, 返回 None。
    """
    if n < 2 or k < 1 or k >= n:
        return None
    lo, hi = 1e-8, max(1e3, float(n) * 10.0)
    for _ in range(200):  # θ̂: 对 E[K] 单调, 几何二分
        mid = (lo * hi) ** 0.5
        e = sum(mid / (mid + i) for i in range(n))
        if e < k:
            lo = mid
        else:
            hi = mid
    theta = (lo * hi) ** 0.5
    P = [0.0] * (n + 1)
    P[1] = 1.0
    for j in range(2, n + 1):
        new = [0.0] * (n + 1)
        repeat = (j - 1) / (theta + j - 1)
        new_allele = theta / (theta + j - 1)
        for kk in range(1, j + 1):
            new[kk] = repeat * P[kk] + (new_allele * P[kk - 1] if kk > 1 else 0.0)
        P = new
    cum = sum(P[1:k + 1])
    return round(math.log(cum), 4) if cum > 0 else None


_ALLELE_IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def _ac_matrix(seqs):
    """ACGT 等位基因计数矩阵 (位点 × 4, int64) —— scikit-allel 的输入格式。

    2026-10-13 (修 10-10 P0-3/P0-5): 非 ACGT (gap/N/简并) 不计入任何等位基因,
    与 _diversity_core 的 pairwise-deletion 口径一致 (只认大写 ACGT)。
    """
    L = len(seqs[0]) if seqs else 0
    ac = np.zeros((L, 4), dtype=np.int64)
    if seqs:
        arr = np.frombuffer(''.join(seqs).encode('ascii'), dtype='S1'
                            ).reshape(len(seqs), L)
        for ch, k in _ALLELE_IDX.items():
            ac[:, k] = (arr == ch.encode('ascii')).sum(axis=0)
    return ac


def _diversity_core(seqs_used):
    """一组序列的 π/θ-W/S/Hd/k/Tajima's D/Fu's Fs (仅 ACGT, pairwise-deletion)。

    数学与原 diversity() 逐字一致 (2026-09-15 审查 P2-15 口径), 只是抽出成
    模块级函数供 全局/按群体 两处调用; Fu's Fs 由恒 None 改为 _fu_fs 实算。
    """
    seqs = [s for _, s in seqs_used]
    n = len(seqs)
    L = len(seqs[0]) if seqs else 0

    # π: pairwise deletion (仅 ACGT 位点参与比较)
    total_pi = 0.0
    comparisons = 0
    for i in range(n):
        for j in range(i + 1, n):
            valid = [(a, b) for a, b in zip(seqs[i], seqs[j])
                     if a in 'ACGT' and b in 'ACGT']
            if valid:
                total_pi += sum(1 for a, b in valid if a != b) / len(valid)
                comparisons += 1
    pi = total_pi / comparisons if comparisons else 0.0

    # S: 仅 ACGT 列的多态位点数
    S = 0
    for pos in range(L):
        col = [s[pos] for s in seqs if s[pos] in 'ACGT']
        if len(col) >= 2 and len(set(col)) > 1:
            S += 1

    # θ-W = S / a1,  a1 = sum(1/i)
    a1 = sum(1.0 / i for i in range(1, n)) if n > 1 else 0.0
    theta_w = S / a1 if a1 > 0 else 0.0

    # Tajima's D (2026-10-13 修 10-10 P0-3): 交由 scikit-allel 计算, 不再手写公式。
    # 旧实现分子混用「每位点 pi」与「绝对 theta_w = S/a1」, 量纲不一致 →
    # 符号/数值可翻转 (10-10 审查实测旧式 -3.31 vs 正确 +1.89)。
    # allel.tajima_d 内部 = 绝对 theta_pi(逐位点成对差异之和) − 绝对 theta_w,
    # 与 S 同量纲 (Nei 1987 标准式); 输入计数矩阵只含 ACGT 列,
    # pairwise-deletion 语义与上方 pi/S 完全一致。
    tajima_d = None
    if S > 0 and n >= 2:
        try:
            from allel.stats.diversity import tajima_d as _allel_tajima_d
        except ImportError as _ie:
            raise SystemExit("[popgen_analysis] Tajima's D 需要 scikit-allel "
                             "(pixi.toml 已登记; pip install scikit-allel): " + str(_ie))
        _d = float(_allel_tajima_d(_ac_matrix(seqs), min_sites=1))
        tajima_d = None if math.isnan(_d) else _d

    # 单倍型
    haps = {}
    for sid, s in seqs_used:
        haps.setdefault(s, []).append(sid)
    nh = len(haps)
    freqs = list(map(len, haps.values()))
    Hd = (1 - sum(f * (f - 1) for f in freqs) / (n * (n - 1))) if n > 1 else 0.0

    return {
        'n': n, 'L': L, 'haplotypes': nh, 'Hd': round(Hd, 6),
        'k': round(pi * L, 4), 'pi': round(pi, 6),
        'theta_w': round(theta_w, 6), 'S': S,
        'tajima_D': round(tajima_d, 4) if tajima_d is not None else None,
        'fu_fs': _fu_fs(n, nh), 'method': 'manual (ACGT pairwise-deletion)',
    }


def _norm_seq_key(seq_id):
    """序列名 → 元数据键 (单一口径)。

    2026-09-15 (审查 P2-16) 修复: 旧版 _read_meta 用 `split('_')[0]`,
    而 _seq_id 对含 '/' 的 id 用 `split('/')[1]` —— 同一个 id 两条路径
    得到的键不同, 元数据静默匹配不上 (分组全落 Other / Fst 分组为空)。
    2026-09-15 (二次修复): 进一步与 data_collector / build_haplo_outputs
    统一到 utils.seq_ids.base_sample_id —— 旧规则对 `CRR1126135.OR489165.1`
    这类带 accession 版本号的 id 返回整串, 与只写 `CRR1126135` 的元数据
    对不上。
    """
    return base_sample_id(seq_id)


class PopGenAnalyzer:
    """群体遗传学分析器 (输入: MAFFT 比对 fasta + 元数据 CSV; 工具链: pypopart)"""

    def __init__(self, fasta_file, metadata_csv=None, name_col='name',
                 exclude_ids=None, auto_exclude_reference=False):
        """
        Parameters
        ----------
        exclude_ids : iterable of str, optional
            要排除的序列 id (精确匹配)。用于剔除参考基因组/外群。
        auto_exclude_reference : bool, default False
            按命名约定 (REFERENCE_ID_PATTERNS) 自动识别并排除参考/外群。
            默认关闭以保持既有数字口径不变; 但会**告警**提示检测到了什么。
        """
        self.fasta_path = fasta_file
        self.aln = load_alignment(str(fasta_file))
        self.sequences = self._read_fasta(fasta_file)
        self.n = len(self.sequences)
        self.L = len(self.sequences[0][1]) if self.sequences else 0
        self.meta = self._read_meta(metadata_csv, name_col) if metadata_csv else {}
        self._network = None

        # ── 参考/外群排除 (P2-15) ──────────────────────────────────────
        self.exclude_ids = set(str(x) for x in (exclude_ids or ()))
        self.detected_reference_ids = [sid for sid, _ in self.sequences
                                       if looks_like_reference(sid)]
        if auto_exclude_reference:
            self.exclude_ids |= set(self.detected_reference_ids)
        self.excluded_ids = [sid for sid, _ in self.sequences if sid in self.exclude_ids]
        self.sequences_used = [(sid, s) for sid, s in self.sequences
                               if sid not in self.exclude_ids]
        self.n_used = len(self.sequences_used)

    def _read_fasta(self, path):
        seqs = []
        cur_id, cur_seq = None, []
        # 2026-09-15: 补 encoding (Linux LANG=C 下 locale 为 ASCII)
        for line in open(path, encoding='utf-8', errors='replace'):
            line = line.strip()
            if line.startswith('>'):
                if cur_id is not None:
                    seqs.append((cur_id, ''.join(cur_seq)))
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.upper())
        if cur_id is not None:
            seqs.append((cur_id, ''.join(cur_seq)))
        return seqs

    def _read_meta(self, path, name_col):
        meta = {}
        with open(path, encoding='utf-8-sig') as f:
            for row in csv.DictReader(f):
                key = _norm_seq_key(row.get(name_col, ''))
                meta[key] = row
        return meta

    def _seq_id(self, full_id):
        # 与 _read_meta 共用同一口径 (P2-16)
        return _norm_seq_key(full_id)

    def _populations(self, group_col):
        """序列 → 分组映射 (标准化省名/宿主; 通用回退见下)

        2026-09-15 (P2-15): 遍历 `sequences_used` (已排除参考/外群), 与
        `fst` / `amova` 里 `zip(self.sequences_used, mask)` 严格同序同长。
        2026-10-12 (真实数据实测修复): 旧版对不在省份/宿主清单里的取值
        一律落 'Other' —— 换任何非 GCVA 数据集 (如国家级 location 的
        PSTVd 公共数据) 分组全塌成单组, Fst 直接 "group size < min_n"。
        现补通用回退: 清单未命中 → 逗号串取第二段 (一级行政区, 如
        "China, Ningxia, Yinchuan" → Ningxia), 无逗号用原值 (如国家名)。
        清单命中的数据集行为不变。
        """
        provs = ['Ningxia', 'Beijing', 'Guangdong', 'Gansu', 'Neimenggu', 'Qinghai']
        pops = {}
        for sid, _ in self.sequences_used:
            mm = self.meta.get(self._seq_id(sid))
            g = mm.get(group_col, '') if mm else ''
            ng = 'Other'
            for p in provs:
                if p in g:
                    ng = p
                    break
            if 'barbarum' in g:
                ng = 'barbarum'
            elif 'ruthenicum' in g:
                ng = 'ruthenicum'
            elif 'chinense' in g:
                ng = 'chinense'
            elif ng == 'Other' and g.strip():
                segs = [x.strip() for x in g.split(',') if x.strip()]
                ng = segs[1] if len(segs) > 1 else segs[0]
            pops[sid] = ng
        return pops

    def _build_network(self):
        if self._network is None:
            from pypopart.algorithms.mjn import MedianJoiningNetwork
            net = MedianJoiningNetwork(distance_method='hamming', epsilon=0).build_network(self.aln)
            self._network = net
        return self._network

    # ════════════════════════════════════════════════════
    # 1. 遗传多样性 + 中性检验
    # ════════════════════════════════════════════════════
    def diversity(self):
        """π, θ-W, S, Hd, k, Tajima's D, Fu's Fs (修正 N 处理; 全局, 不分组)

        修复: pypopart 的 calculate_tajimas_d 只删 '-' 不删 'N',
        导致 N 被当作多态碱基 (S/π 高估)。这里手动实现正确的
        pairwise-deletion 统计 (仅 ACGT)。数学在模块级 `_diversity_core`
        (2026-10-12 抽出, 与按群体版本共用同一实现)。

        2026-09-15 (审查 P2-15): 统计只针对 `self.sequences_used`
        (即已排除参考/外群的集合, 见 __init__ 的 exclude_ids /
        auto_exclude_reference)。若检测到疑似参考序列但仍被计入, 会告警。
        """
        seqs_used = self.sequences_used
        n = len(seqs_used)

        # 参考/外群未排除时的显式告警 (不静默改数字, 但必须让人看见)
        _unexcluded = [sid for sid in self.detected_reference_ids
                       if sid not in self.exclude_ids]
        if _unexcluded:
            print(f"[popgen_analysis] 警告: 比对中疑似参考/外群序列未排除 "
                  f"({len(_unexcluded)} 条: {_unexcluded[:3]}), π/S/θ/Hd 会被系统性抬高。"
                  f" 需要排除请用 --exclude-reference 或 --exclude-ids", file=sys.stderr)
        if self.excluded_ids:
            print(f"[popgen_analysis] 已排除 {len(self.excluded_ids)} 条序列: "
                  f"{self.excluded_ids[:5]}", file=sys.stderr)

        d = _diversity_core(seqs_used)
        # 口径透明化 (P2-15): 让下游能看出这个 n 是否已剔除参考/外群
        d['n_excluded'] = len(self.excluded_ids)
        d['excluded_ids'] = list(self.excluded_ids)
        d['detected_reference_ids'] = list(self.detected_reference_ids)
        return d

    def diversity_by_population(self, group_col, min_n=2):
        """按群体拆分的多样性 + 中性检验 (2026-10-12 新增, 配套可视化)。

        复用 `_populations` 分组映射 (与 fst/amova 严格同序同长), 每个群体
        n≥min_n (默认 2, Hd 需要 ≥2) 单独跑 `_diversity_core` —— 口径与
        全局 diversity() 完全一致。返回按 n 降序的列表。
        """
        pops = self._populations(group_col)
        buckets = {}
        for sid, s in self.sequences_used:
            buckets.setdefault(pops.get(sid, 'Other'), []).append((sid, s))
        out = []
        for pop, members in sorted(buckets.items(), key=lambda kv: -len(kv[1])):
            if len(members) < max(2, min_n):
                continue
            d = _diversity_core(members)
            d['pop'] = pop
            out.append(d)
        return out

    # ════════════════════════════════════════════════════
    # 2. 遗传分化 Fst
    # ════════════════════════════════════════════════════
    def fst(self, group_col, min_n=5, n_perm=999, seed=42):
        """Fst: pypopart 单倍型频率法 (Nei) + 位点法补充 (permutation)"""
        pops = self._populations(group_col)
        tab = Counter(pops.values())
        keep = [g for g, c in tab.items() if c >= min_n]
        if len(keep) < 2:
            return {'fst': None, 'p': None, 'groups': dict(tab), 'error': 'group size < min_n'}

        # 方法1: Nei 单倍型频率法 Fst (DnaSP/Arlequin 同款标准公式)
        # 注: pypopart 的 calculate_fst_matrix(network) 因 construct_network 从不传
        # population_map → 单倍型全 Unassigned, 无法按群体算 Fst (工具缺陷);
        # 这里用带群体映射的 identify_haplotypes_from_alignment 自行按标准公式计算。
        hap_fst = {}
        try:
            from pypopart.core.haplotype import identify_haplotypes_from_alignment
            haps = identify_haplotypes_from_alignment(self.aln, population_map=pops)
            pop_sizes = Counter(g for g in pops.values() if g in keep)
            N = sum(pop_sizes.values())
            if N >= 2 and len(pop_sizes) >= 2:
                # Ht: 全样本单倍型多样度
                total_hap = Counter()
                for h in haps:
                    for pop, cnt in h.get_frequency_by_population().items():
                        if pop in keep:
                            total_hap[h.id] += cnt
                Ht = 1 - sum((c / N) ** 2 for c in total_hap.values())
                # Hs: 群体内多样度加权平均
                Hs = 0.0
                for pop, size in pop_sizes.items():
                    within = Counter()
                    for h in haps:
                        within[h.id] += h.get_frequency_by_population().get(pop, 0)
                    hs_g = 1 - sum((c / size) ** 2 for c in within.values()) if size > 0 else 0.0
                    Hs += (size / N) * hs_g
                if Ht > 0:
                    hap_fst['Nei_overall'] = round((Ht - Hs) / Ht, 4)
            # 网络构建 (pypopart 已 patch: 单倍型序列保留完整比对等长;
            # 网络仅用于单倍型关系/可视化, 不参与 Fst 数值)
            self._build_network()
        except Exception as _hap_err:
            # 2026-09-15 (审查 P2-14) 修复: 旧版裸 `except Exception: hap_fst = {}`
            # 静默吞异常 —— 单倍型法 Fst 失败时输出里什么都没有, 看起来就像
            # "没有单倍型法结果"而不是"算崩了"。
            hap_fst = {}
            print(f"[popgen_analysis] 警告: 单倍型法 Fst 计算失败 ({type(_hap_err).__name__}: "
                  f"{_hap_err}); 已回退为仅位点法 Fst", file=sys.stderr)

        # 方法2: 位点法 Fst + permutation
        # 2026-10-13 (修 10-10 P0-5) 重写, 公式实现全部交由 scikit-allel:
        # 旧版 fst_one 列内剔 gap 后群体标签向量未同步剔除 → 含 gap 的位点
        # 标签错位 (10-10 实测 Fst=1 vs 正确 1/3)。现从「群体×位点×等位基因」
        # 计数矩阵出发 (gap 序列在计数中天然缺席, 标签不可能错位), 估计量用
        # Hudson (1992)/Bhatia (2013): 全局值 = 各群体两两 Hudson Fst 的
        # 比值和聚合。选 Hudson 而非 Weir&Cockerham 是因为 wc 实现只支持
        # 二倍体 (allel 硬性限制), 病毒数据是单倍体; Hudson 直接吃等位基因
        # 计数, 任意倍性皆可。旧 Nei 1−Hs/Ht 逐位点均值口径废弃
        # (单倍型法 'fst_haplotype' 保留)。
        try:
            from allel.stats.fst import hudson_fst as _hudson_fst
        except ImportError as _ie:
            raise SystemExit("[popgen_analysis] 位点法 Fst 需要 scikit-allel "
                             "(pixi.toml 已登记; pip install scikit-allel): " + str(_ie))

        rng = np.random.RandomState(seed)
        mask = np.isin(list(pops.values()), keep)
        _sub_seqs = [s for (sid, s), m in zip(self.sequences_used, mask) if m]
        _sub_labs = np.array([pops[sid] for (sid, s), m in zip(self.sequences_used, mask) if m])
        # 每样本的「位点 × 等位基因」计数 (int8 堆叠, n×L×4);
        # 按组求和即得组计数矩阵 —— gap 序列在该样本天然全零, 标签不可能错位。
        _stack = np.stack([_ac_matrix([s]) for s in _sub_seqs]
                          ).astype(np.int8) if _sub_seqs else None
        _ac_all = _ac_matrix(_sub_seqs)
        n_sites_used = int(((_ac_all > 0).sum(axis=1) >= 2).sum())
        if n_sites_used == 0 or _stack is None:
            return {'fst': None, 'p': None, 'error': 'no polymorphic sites'}

        def _fst_from_labels(labs):
            acs = {}
            for g in keep:
                idx = np.flatnonzero(labs == g)
                if len(idx):
                    acs[g] = _stack[idx].sum(axis=0)
            if len(acs) < 2:
                return np.nan
            gs = sorted(acs)
            nums, dens = [], []
            for i in range(len(gs)):
                for j in range(i + 1, len(gs)):
                    num, den = _hudson_fst(acs[gs[i]], acs[gs[j]])
                    ok = np.isfinite(den) & (den > 0)
                    if ok.any():
                        nums.append(float(np.sum(num[ok])))
                        dens.append(float(np.sum(den[ok])))
            return float(np.sum(nums) / np.sum(dens)) if dens else np.nan

        obs_fst = _fst_from_labels(_sub_labs)
        if not np.isfinite(obs_fst):
            return {'fst': None, 'p': None, 'error': 'no polymorphic sites'}
        perms = []
        for _ in range(n_perm):
            _f = _fst_from_labels(rng.permutation(_sub_labs))
            if np.isfinite(_f):
                perms.append(_f)
        pval = (sum(1 for x in perms if x >= obs_fst) + 1) / (len(perms) + 1)

        return {
            'fst_haplotype': hap_fst,      # Nei 单倍型频率法
            'fst_site': obs_fst, 'p': pval,  # 位点法 (Hudson 两两均值) + permutation
            'groups': {g: int(tab[g]) for g in keep},
            'sites': n_sites_used, 'perm': len(perms),
        }

    # ════════════════════════════════════════════════════
    # 3. AMOVA (群体遗传方差分层)
    # ════════════════════════════════════════════════════
    def fst_pairwise(self, group_col, min_n=5):
        """两两群体间位点法 Fst (2026-10-12 新增, 配套 Fst 热图)。

        2026-10-13 (修 10-10 P0-5 同源问题) 改用 scikit-allel 的 Hudson (1992)
        估计量 (Bhatia et al. 2013), 公式不手写: 旧实现 pooled Ht 分母用整组
        大小, 当位 gap 序列未从分母剔除 → 有 gap 的位点被系统性稀释; 汇总
        口径也由逐位点均值改为比值和 (Σnum/Σden), 与 fst() 的 WC 口径同风格。
        只算观测值 (无置换) —— 显著性仍看全局 fst() 的置换 p, 热图读数值即可。
        返回 {'pops': [...], 'matrix': [[...]]} (对角线 0, n<min_n 的组不入)。
        """
        pops = self._populations(group_col)
        tab = Counter(pops.values())
        keep = sorted(g for g, c in tab.items() if c >= min_n)
        if len(keep) < 2:
            return {'pops': keep, 'matrix': [], 'error': 'group size < min_n'}

        try:
            from allel.stats.fst import hudson_fst as _hudson_fst
        except ImportError as _ie:
            raise SystemExit("[popgen_analysis] Fst 热图需要 scikit-allel "
                             "(pixi.toml 已登记; pip install scikit-allel): " + str(_ie))

        by_group = {g: [] for g in keep}
        for sid, s in self.sequences_used:
            g = pops.get(sid)
            if g in by_group:
                by_group[g].append(s)
        # 每组的「位点 × 等位基因」计数矩阵只算一次 (gap 序列天然不计入)
        acs = {g: _ac_matrix(ss) for g, ss in by_group.items() if ss}

        def _pair_fst(ac1, ac2):
            # 2026-10-13 修 (用户报告 E 面板大面积空白/数值虚高): hudson_fst 是
            # **双等位**估计量 (要求 (n,2) 计数), 直接喂 4 列 ACGT 矩阵会算出
            # 错误数值且多数配对 NaN → None (实测 Ningxia 5 对全消失, 余下数值
            # 虚高 ~17 倍)。现按「两群体合并计数恰为 2 个 ACGT 等位」筛选双等位
            # 位点, 多等位位点保守跳过 (标准做法), 再喂 Hudson。
            pooled = ac1 + ac2
            mask = (pooled > 0).sum(axis=1) == 2
            if not mask.any():
                return None
            idx = np.argsort(-pooled[mask], axis=1)[:, :2]   # 每位点两等位 (按 pooled 计数序, 群间一致)
            b1 = np.take_along_axis(ac1[mask], idx, axis=1)
            b2 = np.take_along_axis(ac2[mask], idx, axis=1)
            num, den = _hudson_fst(b1, b2)
            # 2026-10-14 修 (用户报告 E 面板 10 对 NaN): 某群体在位点上全 gap 时
            # 计数行全零 → Hudson 内部 0/0 = NaN。过滤: 两群在该位点都必须有计数。
            tot1 = b1.sum(axis=1)
            tot2 = b2.sum(axis=1)
            ok = (np.isfinite(num) & np.isfinite(den) & (den > 0)
                  & (tot1 > 0) & (tot2 > 0))
            if not ok.any():
                return None
            return round(float(np.sum(num[ok]) / np.sum(den[ok])), 4)

        mat = [[0.0] * len(keep) for _ in keep]
        for i, g1 in enumerate(keep):
            for j, g2 in enumerate(keep):
                if j <= i or g1 not in acs or g2 not in acs:
                    continue
                v = _pair_fst(acs[g1], acs[g2])
                mat[i][j] = v
                mat[j][i] = v
        return {'pops': keep, 'matrix': mat}

    # ════════════════════════════════════════════════════
    # 2b. 私有单倍型 / 单倍型×群体关联 / 相似性聚类
    #     (2026-10-12 新增, 借 EasyHap research.py/stats.py 思想;
    #      配套 popgen_stats.json 可视化与摘要)
    # ════════════════════════════════════════════════════

    def _hap_freq_table(self):
        """单倍型频率表: [(hap_id, seq, count, [member_sid...]), ...] 按频次降序。"""
        buckets = {}
        for sid, s in self.sequences_used:
            buckets.setdefault(s, []).append(sid)
        table = []
        for rank, (seq, members) in enumerate(
                sorted(buckets.items(), key=lambda kv: -len(kv[1])), 1):
            table.append((f"H{rank}", seq, len(members), members))
        return table

    def private_haplotypes(self, group_col):
        """私有单倍型: 仅在一个群体出现的单倍型 (借 EasyHap
        write_population_haplotype_statistics 的 PrivateHaplotypes 表思想)。

        私有单倍型是地域性传播/本地循环的最直观证据 (如 "Ningxia 私有单倍型
        41 个 vs Neimenggu 0 个"), 比 Fst 更可读。返回:
        {'n_private', 'n_total_haps', 'by_group': {g: {n_haps,count,freq_max}},
         'rows': [{group,hap_id,count,freq}...]} (rows 按频次降序)。
        """
        pops = self._populations(group_col)
        groups_in_hap = {}
        for hap_id, seq, cnt, members in self._hap_freq_table():
            gs = {pops.get(m, 'Other') for m in members}
            if len(gs) == 1:
                g = next(iter(gs))
                groups_in_hap.setdefault(g, []).append((hap_id, cnt))
        n_by_g = Counter(pops.values())
        by_group = {}
        rows = []
        for g, haps in groups_in_hap.items():
            by_group[g] = {'n_haps': len(haps), 'count': sum(c for _, c in haps),
                           'freq_max': round(max(c for _, c in haps) /
                                             max(1, n_by_g.get(g, 1)), 4)}
            for hap_id, cnt in sorted(haps, key=lambda x: -x[1]):
                rows.append({'group': g, 'hap_id': hap_id, 'count': cnt,
                             'freq': round(cnt / max(1, n_by_g.get(g, 1)), 4)})
        rows.sort(key=lambda r: -r['count'])
        return {'n_private': len(rows),
                'n_total_haps': len(set(s for _, s in self.sequences_used)),
                'by_group': by_group, 'rows': rows}

    def haplotype_enrichment(self, group_col, min_hap_count=3):
        """主要单倍型 × 群体富集: 每单倍型对每群体做 2×2 Fisher 精确检验,
        全部检验统一 BH 校正 (借 EasyHap write_trait_association +
        stats.fisher_exact_2x2/BH 思想; 用 scipy.stats.fisher_exact 实现)。

        回答 "哪个单倍型在哪个群体富集" —— Fst/AMOVA 是全基因组汇总,
        本检验给单倍型级指认。返回 {'n_tests', 'rows':
        [{hap_id,count,group,in_group,group_n,p,q,star}...]} 按 q 升序。
        """
        from scipy.stats import fisher_exact

        pops = self._populations(group_col)
        N = len(self.sequences_used)
        n_by_g = Counter(pops.values())
        tests = []
        for hap_id, seq, cnt, members in self._hap_freq_table():
            if cnt < max(2, min_hap_count):
                continue
            mem_g = Counter(pops.get(m, 'Other') for m in members)
            for g in sorted(set(pops.values())):
                a = mem_g.get(g, 0)               # 组内 & 有该单倍型
                b = cnt - a                        # 组外 & 有
                c = n_by_g.get(g, 0) - a          # 组内 & 无
                d = N - a - b - c                 # 组外 & 无
                _, p = fisher_exact([[a, b], [c, d]])
                tests.append({'hap_id': hap_id, 'count': cnt, 'group': g,
                              'in_group': a, 'group_n': n_by_g.get(g, 0), 'p': float(p)})
        # BH 校正
        tests.sort(key=lambda t: t['p'])
        m = len(tests)
        q_prev = 1.0
        for i in range(m - 1, -1, -1):
            q = min(q_prev, tests[i]['p'] * m / (i + 1))
            tests[i]['q'] = round(q, 6)
            q_prev = q
        for t in tests:
            t['star'] = ('***' if t['q'] < 0.001 else '**' if t['q'] < 0.01
                         else '*' if t['q'] < 0.05 else '')
        return {'n_tests': m, 'rows': tests}

    def similarity_clusters(self, group_col, max_diff_ratio: float = 0.03):
        """Hamming 距离阈值的相似性连通分量 (借 EasyHap
        stats.connected_component_clusters 思想; numpy 分块向量化)。

        描述性 strain 群: 簇内两两差异比例 ≤ max_diff_ratio (默认 3%)。
        返回 {'n_clusters', 'threshold', 'assignment': {sid: cluster_id},
        'composition': [{cluster, n, group_sizes:{g:c}}...]} (簇按大小降序)。
        """
        seqs = [s for _, s in self.sequences_used]
        sids = [sid for sid, _ in self.sequences_used]
        n = len(seqs)
        pops = self._populations(group_col)
        arr = np.frombuffer(''.join(seqs).encode('ascii'), dtype=np.uint8).reshape(n, -1)
        parent = list(range(n))

        def _find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        chunk = max(1, min(64, int(4e7 // max(1, arr.shape[1]))))
        for i in range(0, n, chunk):
            bi = arr[i:i + chunk]
            for j in range(i, n, chunk):
                bj = arr[j:j + chunk]
                diff = (bi[:, None, :] != bj[None, :, :]).mean(axis=2)
                for a in range(bi.shape[0]):
                    for b in range(bj.shape[0]):
                        if j + b > i + a and diff[a, b] <= max_diff_ratio:
                            ra, rb = _find(i + a), _find(j + b)
                            if ra != rb:
                                parent[ra] = rb
        roots = {}
        for x in range(n):
            r = _find(x)
            roots.setdefault(r, []).append(x)
        comp = sorted(roots.values(), key=len, reverse=True)
        assignment, composition = {}, []
        for cid, members in enumerate(comp, 1):
            cname = f"C{cid}"
            gs = Counter(pops.get(sids[x], 'Other') for x in members)
            for x in members:
                assignment[sids[x]] = cname
            composition.append({'cluster': cname, 'n': len(members),
                                'group_sizes': dict(gs.most_common())})
        return {'n_clusters': len(comp), 'threshold': max_diff_ratio,
                'assignment': assignment, 'composition': composition}

    # ════════════════════════════════════════════════════
    # 3. AMOVA (群体遗传方差分层)
    # ════════════════════════════════════════════════════
    def amova(self, group_col, min_n=5):
        pops = self._populations(group_col)
        tab = Counter(pops.values())
        keep = [g for g, c in tab.items() if c >= min_n]
        if len(keep) < 2:
            return {'error': 'group size < min_n', 'groups': dict(tab)}
        network = self._build_network()
        # 设置 population (group_col 分组), 供 AMOVA 方差分层
        for hap_id in network.nodes:
            hap = network.get_haplotype(hap_id)
            if hap is None:
                continue
            for sid in list(hap.sample_ids):
                pop = pops.get(sid, pops.get(_norm_seq_key(sid), 'Unknown'))
                hap.add_sample(sid, pop)
        try:
            r = pg_stats.calculate_amova(network, self.aln)
            return {
                'within': getattr(r, 'variance_within_pops', None),
                'among': getattr(r, 'variance_among_pops', None),
                'phi_st': getattr(r, 'phi_st', None),
                'percent_within': getattr(r, 'percent_within_pops', None),
                'percent_among': getattr(r, 'percent_among_pops', None),
                'groups': dict(tab),
            }
        except Exception as e:
            return {'error': str(e)}

    # ════════════════════════════════════════════════════
    # 4. 保守性图谱 (Shannon 熵)
    # ════════════════════════════════════════════════════
    def sliding_tracks(self, window=100, step=50, group_col=None,
                       min_n=2, reference=None):
        """滑窗 π / Fst 轨道 (2026-10-22 新增, 供单倍型状态热图轨道层)。

        与 hastat view -t pi/-t fst 的滑窗口径对齐 (论文 Fig.2 构图: π 按群体
        各一条线, Fst 按群体两两各一条线): 窗口按**比对列**滑动, 每个窗口对
        每个群体算一次组内 π (pairwise deletion, 仅 ACGT), 对每个群体对算一次
        Hudson (1992) Fst (scikit-allel; 双等位列筛选, 多等位保守跳过)。
        窗口数过多时按 step 倍增降采样, 保证轨道点数 ≤ 400 (与热图 max_sites
        同量级)。

        返回 {'pi': {group: [(pos, val), ...]},
               'fst': {pair_label: [(pos, val), ...]}} —— pair_label 形如
        "A|B" (群体名按字典序, '|' 分隔, 避免群体名本身含下划线时分不清)。
        无元数据/无群体时 'pi' 退化为 {'all': [...]} (全体序列一个 π 轨道),
        'fst' 为空 dict。pos 为窗口中心的**参考坐标** (1-based, 经比对列→参考
        坐标映射); 缺参考映射时退化为比对列坐标 (调用方需知情: 无 genes 时
        热图横轴也是列坐标)。
        """
        seqs = [s for _, s in self.sequences_used]
        n = len(seqs)
        out = {'pi': {}, 'fst': {}}
        if n < 2 or self.L < window:
            out['pi'] = {'all': []}
            return out
        # 比对列 → 参考坐标 (有参考可比时; 失败则退化为列坐标)
        col2ref = None
        if reference is not None:
            try:
                from utils.seq_match import resolve_reference
                _rn, ref_seq = resolve_reference(self.sequences_used, reference,
                                                 who="popgen.sliding_tracks")
                col2ref = {}
                seen = 0
                for i, c in enumerate(ref_seq):
                    if c != '-':
                        seen += 1
                        col2ref[i] = seen
            except Exception as e:  # noqa: BLE001
                print(f"[popgen] 滑窗轨道参考映射失败, 退化为列坐标: {e}",
                      file=sys.stderr)
        # 群体划分 (π 按群体 / Fst 按群体对); 无元数据时 π 只有 'all' 一个轨道
        pops = {}
        if group_col:
            try:
                pops = self._populations(group_col)
            except Exception as e:  # noqa: BLE001
                print(f"[popgen] 滑窗轨道分组失败, 仅出全体 π: {e}", file=sys.stderr)
        tab = Counter(pops.values())
        keep = sorted(g for g, c in tab.items() if c >= max(1, min_n))
        by_group = {g: [] for g in keep}
        if pops:
            for sid, s in self.sequences_used:
                g = pops.get(sid)
                if g in by_group:
                    by_group[g].append(s)
        multi = len(by_group) >= 2
        # 单群体也要出该群体的 π 轨道 (论文 B 面板下方只有一条 π 线的情形)
        pi_keys = keep if keep else ['all']
        for g in keep:
            out['pi'][g] = []
        if not keep:
            out['pi']['all'] = []
        pair_keys = []
        if multi:
            for a_i in range(len(keep)):
                for b_i in range(a_i + 1, len(keep)):
                    pair_keys.append((keep[a_i], keep[b_i]))
                    out['fst'][f"{keep[a_i]}|{keep[b_i]}"] = []

        def _pos_of(col):
            if col2ref is not None and col in col2ref:
                return col2ref[col]
            return col + 1

        def _pi_of(win_seqs):
            m = len(win_seqs)
            if m < 2:
                return None
            npairs = ndiff = 0
            for i in range(m - 1):
                for j in range(i + 1, m):
                    npairs += 1
                    ndiff += sum(1 for a, b in zip(win_seqs[i], win_seqs[j])
                                 if a != b and a in 'ACGT' and b in 'ACGT')
            return (ndiff / npairs / window) if npairs else None

        # 窗口数上限 400: step 倍增直到满足 (窗口不重叠也不稀疏到没信息)
        eff_step = max(1, step)
        while (self.L - window) // eff_step + 1 > 400 and eff_step < self.L:
            eff_step *= 2
        try:
            from allel.stats.fst import hudson_fst as _hudson_fst
        except ImportError:
            _hudson_fst = None
            if pair_keys:
                print("[popgen] 滑窗 Fst 轨道需要 scikit-allel, 仅出 π",
                      file=sys.stderr)
                pair_keys = []
        for start in range(0, self.L - window + 1, eff_step):
            end = start + window
            pos = _pos_of((start + end) // 2)
            win = [s[start:end] for s in seqs]
            for g in pi_keys:
                if g == 'all':
                    v = _pi_of(win)
                else:
                    v = _pi_of([s[start:end] for s in by_group.get(g, [])])
                if v is not None:
                    out['pi'][g].append((pos, round(v, 6)))
            for g1, g2 in pair_keys:
                if _hudson_fst is None:
                    continue
                ac1 = _ac_matrix([s[start:end] for s in by_group[g1]])
                ac2 = _ac_matrix([s[start:end] for s in by_group[g2]])
                pooled = ac1 + ac2
                mask = (pooled > 0).sum(axis=1) == 2
                if not mask.any():
                    continue
                idx = np.argsort(-pooled[mask], axis=1)[:, :2]
                b1 = np.take_along_axis(ac1[mask], idx, axis=1)
                b2 = np.take_along_axis(ac2[mask], idx, axis=1)
                num, den = _hudson_fst(b1, b2)
                ok = (np.isfinite(num) & np.isfinite(den) & (den > 0)
                      & (b1.sum(axis=1) > 0) & (b2.sum(axis=1) > 0))
                if ok.any():
                    v = float(np.sum(num[ok]) / np.sum(den[ok]))
                    out['fst'][f"{g1}|{g2}"].append(
                        (pos, round(max(0.0, min(1.0, v)), 4)))
        # 空的轨道键删掉 (绘图侧按空列表跳过, 但干净点更好)
        out['pi'] = {k: v for k, v in out['pi'].items() if v}
        out['fst'] = {k: v for k, v in out['fst'].items() if v}
        return out

    def hap_group_map(self, group_col, min_n=1):
        """单倍型 → 主导群体 映射 (2026-10-22, 供热图行色条)。

        按该单倍型各成员序列的群体计数取多数; 并列时取字典序第一个 (确定性)。
        只收 ≥ min_n 条的群体; 单倍型无元数据/无匹配群体 → 不出现在映射里
        (绘图侧该行无色条, 与论文里未着色行一致)。
        """
        pops = self._populations(group_col)
        out = {}
        for hid, _seq, _cnt, members in self._hap_freq_table():
            gs = Counter(pops.get(m, '') for m in members if pops.get(m))
            if not gs:
                continue
            g, c = sorted(gs.items(), key=lambda kv: (-kv[1], kv[0]))[0]
            if c >= max(1, min_n):
                out[hid] = g
        return out

    def entropy_profile(self, out_tsv=None):
        # 2026-09-15 (P2-15): 与 diversity 同口径, 用已排除参考/外群的序列集
        seqs = [s for _, s in self.sequences_used]
        ent = []
        for j in range(self.L):
            col = [s[j] for s in seqs if s[j] != '-']
            if len(col) < 2:
                ent.append(0.0)
                continue
            cnt = Counter(col)
            m = len(col)
            H = -sum((f / m) * math.log2(f / m) for f in cnt.values())
            ent.append(H)
        ent = np.array(ent)
        if out_tsv:
            with open(out_tsv, 'w') as f:
                f.write('position\tentropy\n')
                for i, e in enumerate(ent):
                    f.write(f'{i+1}\t{e:.4f}\n')
        return {
            'sites': self.L, 'variable': int((ent > 0).sum()),
            'high_entropy': int((ent > 0.5).sum()), 'entropy': ent,
        }


    # ════════════════════════════════════════════════════
    # 5. 基因变异映射 (变异热点 → 基因区)
    # ════════════════════════════════════════════════════
    def gene_mapping(self, genes_csv=None, genes=None, reference=None):
        """变异位点 → 基因区聚合 (每基因: 变异位点数/高频变异/π/每kb密度)

        genes: list of (gene_name, start, end) 参考坐标 (1-based)
        genes_csv: 或提供 CSV (gene,start,end)
        reference: 坐标基准序列名; 缺省取第一条并告警 (2026-09-15 新增)
        """
        if genes is None:
            genes = []
            if genes_csv and os.path.exists(genes_csv):
                with open(genes_csv) as f:
                    for row in csv.DictReader(f):
                        genes.append((row['gene'], int(row['start']), int(row['end'])))
        if not genes:
            return {'error': 'no gene coordinates'}

        seqs = [s for _, s in self.sequences]
        # 比对列 → 参考坐标映射
        # 2026-09-15 修复: 旧版无条件 `ref_seq = seqs[0]` —— MAFFT 输出顺序不保证
        # 参考在首位, 基准一旦漂移, 全部基因坐标整体错位且无任何提示。
        # 现支持显式 reference (三级匹配), 缺省时明确告警。
        from utils.seq_match import resolve_reference
        ref_name, ref_seq = resolve_reference(
            self.sequences, reference, who="popgen.gene_mapping")
        ref_pos = {}
        seen = 0
        for i, c in enumerate(ref_seq):
            if c != '-':
                seen += 1
                ref_pos[seen] = i
        max_ref = max(ref_pos)

        # 每列变异计数 (共享变异 ≥2)
        site_var = []
        for i in range(self.L):
            col = [s[i] for s in seqs if s[i] != '-']
            if len(col) < 2:
                continue
            cnt = Counter(col)
            nvar = len(col) - cnt.most_common(1)[0][1]
            if nvar >= 2:
                site_var.append((i, nvar))

        results = []
        for g, start, end in genes:
            if start not in ref_pos or end not in ref_pos:
                continue
            cs, ce = ref_pos[start], ref_pos[end] + 1
            var_in = [i for i, n in site_var if cs <= i < ce]
            hi_in = [i for i, n in site_var if cs <= i < ce and n >= 5]
            # π per site
            total_pairs = len(seqs) * (len(seqs) - 1) / 2
            pi_sum = 0.0
            for i in range(cs, ce):
                col = [s[i] for s in seqs if s[i] != '-']
                if len(col) < 2:
                    continue
                cnt = Counter(col)
                m = len(col)
                pi_sum += (m * m - sum(f * f for f in cnt.values())) / 2
            ncol = ce - cs
            pi = (pi_sum / total_pairs / ncol) if ncol else 0.0
            length = end - start + 1
            results.append({
                'gene': g, 'start': start, 'end': end, 'length': length,
                'var_sites': len(var_in), 'high_freq': len(hi_in),
                'pi': round(pi, 6), 'per_kb': round(len(var_in) / length * 1000, 2),
            })
        return {'genes': results,
                'total_var': len(site_var), 'total_high': sum(1 for _, n in site_var if n >= 5)}


def main():
    ap = argparse.ArgumentParser(description='病毒群体遗传学分析 (pypopart)')
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--metadata', default=None)
    ap.add_argument('--group', default='location')
    ap.add_argument('--min-n', type=int, default=5)
    ap.add_argument('--perm', type=int, default=999)
    # 2026-09-15 (审查 P2-15): 参考/外群排除开关
    ap.add_argument('--exclude-reference', action='store_true',
                    help='按命名约定自动排除参考基因组/外群序列 (默认关闭)')
    ap.add_argument('--exclude-ids', default=None,
                    help='额外排除的序列 id, 逗号分隔 (精确匹配)')
    ap.add_argument('--fst-site', type=float, default=None,
                    help='[2026-10-12] 注入全局位点法 Fst 值 (配合 --skip-global-fst, '
                         '来源=历史报告/论文口径, json 标注 source=manual)')
    ap.add_argument('--fst-p', type=float, default=None,
                    help='注入全局 Fst 的置换 p 值 (配合 --fst-site)')
    ap.add_argument('--skip-global-fst', action='store_true',
                    help='[2026-10-12] 跳过全局 Fst (含 MJN 网络) —— 近全唯一单倍型数据集'
                         '(病毒全长常见)的 MJN 准中介枚举是指数爆炸 (实测 PSTVd 348GB / '
                         'GCVA 32GB 帽被杀, 均 _quasi_medians MemoryError); pairwise '
                         'Fst 矩阵与按群体统计不受影响。全局值可沿用历史报告。')
    ap.add_argument('--genes-csv', default=None,
                    help='[2026-10-12] 基因坐标 CSV (gene,start,end), 供单倍型状态热图画基因结构带')
    ap.add_argument('--window', type=int, default=100,
                    help='[2026-10-22] 单倍型热图 π/Fst 滑窗轨道窗口 (bp, 比对列); 0 或负值关闭轨道')
    ap.add_argument('--step', type=int, default=50,
                    help='[2026-10-22] 滑窗步长 (bp); 窗口数超 400 时自动倍增降采样')
    ap.add_argument('--no-tracks', action='store_true',
                    help='[2026-10-22] 关闭 π/Fst 滑窗轨道 (只出基因带+热图)')
    ap.add_argument('--stats-json', default=None,
                    help='[2026-10-12] 额外输出结构化统计 JSON (全局+按群体+Fst 矩阵), '
                         '供 utils/popgen_plots.py 渲染多面板统计图')
    args = ap.parse_args()

    _excl = [x.strip() for x in (args.exclude_ids or '').split(',') if x.strip()]
    pg = PopGenAnalyzer(args.fasta, args.metadata,
                        exclude_ids=_excl,
                        auto_exclude_reference=args.exclude_reference)
    d = pg.diversity()
    print(f"序列: {d['n']}, 位点: {d['L']}, 单倍型: {d['haplotypes']} (Hd={d['Hd']:.4f})")
    if d.get('n_excluded'):
        print(f"  (已排除 {d['n_excluded']} 条: {d.get('excluded_ids')})")
    print(f"π = {d['pi']:.6f} ({d['pi']*100:.3f}%)  k = {d['k']:.2f}  S = {d['S']}  θ-W = {d['theta_w']:.6f}")
    # 2026-09-15 修复: S==0 (无多态位点) 时 tajima_D 为 None, 原 f-string 直接
    # `{d['tajima_D']:.4f}` → TypeError 直接崩掉 CLI。
    _td = d.get('tajima_D')
    _td_s = f"{_td:.4f}" if _td is not None else "NA (S=0, 无多态位点)"
    print(f"Tajima's D = {_td_s}  Fu's Fs = {d.get('fu_fs')}")

    f = None
    if args.metadata and args.skip_global_fst:
        print("[popgen_analysis] --skip-global-fst: 跳过全局 Fst (MJN 爆炸防护), "
              "全局值请沿用历史报告或fst_pairwise")
    elif args.metadata:
        f = pg.fst(args.group, min_n=args.min_n, n_perm=args.perm)
        if f.get('fst_site') is not None:
            print(f"Fst(位点法) = {f['fst_site']:.4f}, p = {f['p']:.3f} (组: {f['groups']})")
            if f.get('fst_haplotype'):
                print(f"Fst(单倍型法, Nei) = {f['fst_haplotype']}")
        else:
            print(f"Fst: {f.get('error')}")

    # 结构化统计输出 (2026-10-12): 全局 + 按群体 (含 Fu's Fs) + pairwise Fst 矩阵
    if args.stats_json:
        out = {'group_col': args.group, 'global': d, 'populations': [], 'fst': None,
               'fst_pairwise': None}
        try:
            out['populations'] = pg.diversity_by_population(args.group, min_n=2)
        except Exception as e:  # noqa: BLE001
            print(f"[popgen_analysis] 按群体统计失败: {e}", file=sys.stderr)
        if args.metadata:
            if f is None and args.fst_site is not None:
                f = {'fst_site': args.fst_site, 'p': args.fst_p,
                     'fst_haplotype': {},
                     'source': 'manual injection (historical report)'}
            out['fst'] = f
            try:
                out['fst_pairwise'] = pg.fst_pairwise(args.group, min_n=args.min_n)
            except Exception as e:  # noqa: BLE001
                print(f"[popgen_analysis] pairwise Fst 失败: {e}", file=sys.stderr)
        # 2026-10-12 借 EasyHap 补充: 私有单倍型 / 单倍型×群体富集 (Fisher+BH) /
        # 相似性聚类。宿主列若在元数据里则额外做一轮 host 富集。
        if args.metadata:
            try:
                out['private_haplotypes'] = pg.private_haplotypes(args.group)
            except Exception as e:  # noqa: BLE001
                print(f"[popgen_analysis] 私有单倍型失败: {e}", file=sys.stderr)
            try:
                out['enrichment'] = pg.haplotype_enrichment(args.group, min_hap_count=3)
                _has_host = any('host' in (r or {}) and (r or {}).get('host')
                                for r in pg.meta.values())
                if _has_host:
                    out['enrichment_host'] = pg.haplotype_enrichment('host', min_hap_count=3)
            except Exception as e:  # noqa: BLE001
                print(f"[popgen_analysis] 单倍型富集失败: {e}", file=sys.stderr)
            try:
                out['clusters'] = pg.similarity_clusters(args.group)
            except Exception as e:  # noqa: BLE001
                print(f"[popgen_analysis] 相似性聚类失败: {e}", file=sys.stderr)
            # 配套 CSV (私有单倍型 / 富集检验 / 聚类组成)
            try:
                _base = os.path.splitext(args.stats_json)[0]
                pv = out.get('private_haplotypes') or {}
                if pv.get('rows'):
                    _df_head = 'group\thap_id\tcount\tfreq'
                    _lines = [_df_head] + [f"{r['group']}\t{r['hap_id']}\t{r['count']}\t{r['freq']}"
                                           for r in pv['rows']]
                    open(_base + '.private_haplotypes.tsv', 'w', encoding='utf-8').write('\n'.join(_lines) + '\n')
                en = out.get('enrichment') or {}
                # 无可检验单倍型时也写空表 (带表头+注释行) —— 与 "没跑过" 可区分,
                # report 声明不再常驻 not-found
                if not en.get('rows'):
                    open(_base + '.haplotype_enrichment.tsv', 'w', encoding='utf-8').write(
                        'hap_id\tcount\tgroup\tin_group\tgroup_n\tp\tq\tstar\n'
                        '# n_tests=0 (no haplotype reached min_hap_count)\n')
                if en.get('rows'):
                    _head = 'hap_id\tcount\tgroup\tin_group\tgroup_n\tp\tq\tstar'
                    _lines = [_head] + [f"{r['hap_id']}\t{r['count']}\t{r['group']}\t{r['in_group']}\t{r['group_n']}\t{r['p']:.4g}\t{r['q']:.4g}\t{r['star']}"
                                        for r in en['rows']]
                    open(_base + '.haplotype_enrichment.tsv', 'w', encoding='utf-8').write(
                        'hap_id\tcount\tgroup\tin_group\tgroup_n\tp\tq\tstar\n'
                        '# n_tests=0 (no haplotype reached min_hap_count)\n')
                # 2026-10-22 修 (PSTVd 冒烟实测): 聚类组成 CSV 在 en.get('rows')
                # 分支内引用未定义的 `cl` → 一旦有可检验单倍型, 整个 CSV 导出块
                # NameError 中止 (私有单倍型/富集 CSV 也连带丢)。改为读
                # out['clusters'] 自身, 与前面 JSON 同源; en 分支只写富集表。
                cl = out.get('clusters') or {}
                if cl.get('composition'):
                    _head = 'cluster\tn\tgroup_sizes'
                    _lines = [_head] + [f"{c['cluster']}\t{c['n']}\t{json.dumps(c['group_sizes'], ensure_ascii=False)}"
                                        for c in cl['composition']]
                    open(_base + '.similarity_clusters.tsv', 'w', encoding='utf-8').write('\n'.join(_lines) + '\n')
            except Exception as e:  # noqa: BLE001
                print(f"[popgen_analysis] 统计 CSV 导出失败: {e}", file=sys.stderr)
        with open(args.stats_json, 'w', encoding='utf-8') as fh:
            json.dump(out, fh, ensure_ascii=False, indent=2)
        print(f"[popgen_analysis] 结构化统计 → {args.stats_json}")
        # 单倍型 FASTA + 等位状态热图 (2026-10-12, 借 EasyHap
        # plot_haplotype_heatmap 思想; 位点多时纯栅格, 基因结构带可选 --genes-csv;
        # 2026-10-22 重构为 hastat 论文 Fig.2 构图: 自上而下
        # Fst 轨道 (按群体两两) → π 轨道 (按群体) → 基因结构带 (下缘变异位点
        # 倒三角) → 单倍型矩阵 (左侧行色条按主导群体着色), 全共享基因组横轴)
        try:
            _base = os.path.splitext(args.stats_json)[0]
            _table = pg._hap_freq_table()
            with open(_base + '.haplotypes.fasta', 'w', encoding='utf-8') as hf:
                for hid, seq, cnt, _m in _table:
                    hf.write(f'>{hid} count={cnt}\n{seq}\n')
            _genes = None
            if getattr(args, 'genes_csv', None) and os.path.exists(args.genes_csv):
                with open(args.genes_csv, encoding='utf-8-sig') as gf:
                    _genes = [(r['gene'], int(r['start']), int(r['end']))
                              for r in csv.DictReader(gf)]
            # 滑窗轨道 (--window/--step 控制; 默认 100/50): π 按群体、Fst 按
            # 群体两两; 无元数据时退化为单一 'all' π 轨道、无 Fst。
            _tracks = None
            _hap_groups = None
            if args.metadata:
                try:
                    _hap_groups = pg.hap_group_map(args.group, min_n=1)
                except Exception as e:  # noqa: BLE001
                    print(f"[popgen_analysis] 单倍型主导群体映射失败 (行无色条): {e}",
                          file=sys.stderr)
            if not getattr(args, 'no_tracks', False) and args.window > 0:
                try:
                    _tracks = pg.sliding_tracks(
                        window=args.window, step=args.step,
                        group_col=args.group if args.metadata else None,
                        min_n=args.min_n)
                except Exception as e:  # noqa: BLE001
                    print(f"[popgen_analysis] 滑窗轨道失败 (热图继续): {e}",
                          file=sys.stderr)
            from utils.popgen_plots import render_hap_state_heatmap
            _hr = render_hap_state_heatmap(
                [(h, s, c) for h, s, c, _m in _table],
                _base + '.hap_state_heatmap.pdf', _base + '.hap_state_heatmap.png',
                title='Haplotype allele states', genes=_genes, tracks=_tracks,
                hap_groups=_hap_groups)
            _tag = []
            if _genes:
                _tag.append(f"{len(_genes)} genes")
            if _tracks:
                _tag.append(f"π×{len(_tracks.get('pi') or {})}/"
                            f"Fst×{len(_tracks.get('fst') or {})} tracks")
            if _hap_groups:
                _tag.append(f"{len(set(_hap_groups.values()))} group colors")
            print(f"[popgen_analysis] 单倍型状态热图 → {_base + '.hap_state_heatmap.png'} "
                  f"({_hr.get('n_haps', 0)} x {_hr.get('n_sites', 0)}"
                  + (f", {', '.join(_tag)}" if _tag else '') + ")")
        except Exception as e:  # noqa: BLE001
            print(f"[popgen_analysis] 单倍型状态热图失败: {e}", file=sys.stderr)


if __name__ == '__main__':
    main()
