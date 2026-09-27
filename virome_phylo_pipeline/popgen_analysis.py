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
        """序列 → 分组映射 (标准化省名/宿主)

        2026-09-15 (P2-15): 遍历 `sequences_used` (已排除参考/外群), 与
        `fst` / `amova` 里 `zip(self.sequences_used, mask)` 严格同序同长。
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
        """π, θ-W, S, Hd, k, Tajima's D, Fu's Fs (修正 N 处理)

        修复: pypopart 的 calculate_tajimas_d 只删 '-' 不删 'N',
        导致 N 被当作多态碱基 (S/π 高估)。这里手动实现正确的
        pairwise-deletion 统计 (仅 ACGT)。

        2026-09-15 (审查 P2-15): 统计只针对 `self.sequences_used`
        (即已排除参考/外群的集合, 见 __init__ 的 exclude_ids /
        auto_exclude_reference)。若检测到疑似参考序列但仍被计入, 会告警。
        """
        import math
        seqs_used = self.sequences_used
        n = len(seqs_used)
        seqs = [s for _, s in seqs_used]
        L = len(seqs[0]) if seqs else 0

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
        a1 = sum(1.0 / i for i in range(1, n))
        theta_w = S / a1 if a1 > 0 else 0.0

        # Tajima's D (Nei 1987, 用 ACGT 统计)
        tajima_d = None
        if S > 0 and n >= 2:
            # D = (pi - theta_w) / sqrt(V)
            a2 = sum(1.0 / (i**2) for i in range(1, n))
            b1 = (n + 1) / (3 * (n - 1))
            b2 = 2 * (n**2 + n + 3) / (9 * n * (n - 1))
            c1 = b1 - 1 / a1
            c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1**2)
            e1 = c1 / a1
            e2 = c2 / (a1**2 + a2)
            var = e1 * S + e2 * S * (S - 1)
            tajima_d = (pi - theta_w) / math.sqrt(var) if var > 0 else 0.0

        # 单倍型
        haps = {}
        for sid, s in self.sequences:
            haps.setdefault(s, []).append(sid)
        nh = len(haps)
        freqs = list(map(len, haps.values()))
        Hd = (1 - sum(f * (f - 1) for f in freqs) / (n * (n - 1))) if n > 1 else 0.0

        return {
            'n': n, 'L': L, 'haplotypes': nh, 'Hd': round(Hd, 6),
            'k': round(pi * L, 4), 'pi': round(pi, 6),
            'theta_w': round(theta_w, 6), 'S': S,
            'tajima_D': round(tajima_d, 4) if tajima_d is not None else None,
            'fu_fs': None, 'method': 'manual (ACGT pairwise-deletion)',
            # 口径透明化 (P2-15): 让下游能看出这个 n 是否已剔除参考/外群
            'n_excluded': len(self.excluded_ids),
            'excluded_ids': list(self.excluded_ids),
            'detected_reference_ids': list(self.detected_reference_ids),
        }

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

        # 方法2: 位点法 Fst + permutation (保留, 与文献对齐)
        rng = np.random.RandomState(seed)
        mask = np.isin(list(pops.values()), keep)
        seqs = [s for (sid, s), m in zip(self.sequences_used, mask) if m]
        gvec = np.array([pops[sid] for (sid, s), m in zip(self.sequences_used, mask) if m])

        poly_cols = []
        for j in range(self.L):
            col = [s[j] for s in seqs if s[j] != '-']
            if len(col) < 2 or len(set(col)) < 2:
                continue
            poly_cols.append(col)

        def fst_one(col, gv):
            cnt = Counter(col)
            m = len(col)
            Ht = 1 - sum((f / m) ** 2 for f in cnt.values())
            if Ht <= 0:
                return None
            Hs = 0.0
            for g in set(gv):
                gi = [i for i in range(m) if gv[i] == g]
                colg = [col[i] for i in gi]
                if len(colg) < 2:
                    continue
                cntg = Counter(colg)
                Hs += (len(colg) / m) * (1 - sum((f / len(colg)) ** 2 for f in cntg.values()))
            return 1 - Hs / Ht

        obs = [x for x in (fst_one(c, gvec) for c in poly_cols) if x is not None]
        if not obs:
            return {'fst': None, 'p': None, 'error': 'no polymorphic sites'}
        obs_fst = float(np.mean(obs))
        perms = []
        for _ in range(n_perm):
            gp = rng.permutation(gvec)
            vals = [x for x in (fst_one(c, gp) for c in poly_cols) if x is not None]
            perms.append(float(np.mean(vals)) if vals else 0.0)
        pval = (sum(1 for x in perms if x >= obs_fst) + 1) / (n_perm + 1)

        return {
            'fst_haplotype': hap_fst,      # Nei 单倍型频率法
            'fst_site': obs_fst, 'p': pval,  # 位点法 + permutation
            'groups': {g: int(tab[g]) for g in keep},
            'sites': len(obs), 'perm': n_perm,
        }

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

    if args.metadata:
        f = pg.fst(args.group, min_n=args.min_n, n_perm=args.perm)
        if f.get('fst_site') is not None:
            print(f"Fst(位点法) = {f['fst_site']:.4f}, p = {f['p']:.3f} (组: {f['groups']})")
            if f.get('fst_haplotype'):
                print(f"Fst(单倍型法, Nei) = {f['fst_haplotype']}")
        else:
            print(f"Fst: {f.get('error')}")


if __name__ == '__main__':
    main()
