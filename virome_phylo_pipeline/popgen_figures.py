#!/usr/bin/env python3
"""popgen_figures.py — 变异热点绘图 (论文风格)

  1. plot_entropy_profile: 逐位点 Shannon 熵谱 (可标注基因区)
  2. plot_gene_variation:  每基因变异密度/π 条形图

用法:
  from popgen_analysis import PopGenAnalyzer
  from popgen_figures import plot_entropy_profile, plot_gene_variation
"""
import os
import math
import sys
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── 发表级 rcParams (nature-figure 规范) ──
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['legend.frameon'] = False

GENE_COLORS = ['#0F4D92', '#3775BA', '#42949E', '#8BCF8B', '#E9A6A1', '#9A4D8E', '#767676']


def _shannon(col):
    valid = [c for c in col if c != '-']
    n = len(valid)
    if n < 2:
        return 0.0
    cnt = Counter(valid)
    return -sum((f / n) * math.log2(f / n) for f in cnt.values())


def compute_entropy(sequences, length):
    ent = np.zeros(length)
    for j in range(length):
        col = [s[j] for s, _ in sequences] if isinstance(sequences[0], tuple) else [s[j] for s in sequences]
        ent[j] = _shannon(col)
    return ent


def _ref_position_to_column_map(ref_seq):
    """参考序列坐标 (1-based, 含端) → 比对列下标 (0-based)。

    2026-09-15 (审查 P2-17) 修复依据: 熵谱的横轴是**比对列**, 而 `genes` 里的
    start/end 是**参考序列坐标**。MAFFT 比对含 gap 时两者不相等 —— 直接
    axvspan(start, end) 会把基因区标注错位, 且 gap 越多偏得越远。
    这里按参考序列逐列扫非 gap 字符重建映射。
    """
    m = {}
    k = 0
    for col, ch in enumerate(ref_seq):
        if ch != '-':
            k += 1
            m[k] = col
    return m


def _find_reference_index(seq_ids):
    """在序列名里找参考序列 (与 popgen_analysis.REFERENCE_ID_PATTERNS 同口径)。"""
    try:
        from popgen_analysis import looks_like_reference
    except ImportError:
        from utils.popgen_analysis import looks_like_reference  # type: ignore
    hits = [i for i, sid in enumerate(seq_ids) if looks_like_reference(sid)]
    return hits[0] if len(hits) == 1 else None


def plot_entropy_profile(analyzer, output, genes=None, title='', ref_name=None):
    """Shannon 熵谱图 (可标注基因区)

    genes : list of (gene, start, end, color), start/end 为 **参考序列坐标** (1-based)。
        只有能确定参考序列 (ref_name 指定或命名约定唯一命中) 时才做坐标→列映射;
        否则不标注基因区并告警 —— 宁可不画, 也不画错位置。
    """
    seqs = [s for _, s in analyzer.sequences_used]
    ids = [sid for sid, _ in analyzer.sequences_used]
    L = analyzer.L
    ent = np.zeros(L)
    for j in range(L):
        col = [s[j] for s in seqs if s[j] != '-']
        ent[j] = _shannon(col)

    fig, ax = plt.subplots(figsize=(7.5, 2.4))
    ax.fill_between(range(1, L + 1), ent, color='#0F4D92', alpha=0.7, lw=0.5)
    ax.axhline(0.5, color='#D9544D', ls='--', lw=0.7, alpha=0.7)
    ax.text(L * 1.005, 0.5, 'high-entropy threshold', fontsize=6, color='#D9544D', va='center')
    ax.set_ylabel('Shannon entropy')
    ax.set_xlim(0, L)
    ax.set_ylim(0, max(ent.max() * 1.15, 1.0))

    # ── 基因标注 (参考坐标 → 比对列, 见 P2-17) ────────────────────────
    _xlabel = 'Alignment column'   # 横轴严格来说是比对列; 无 gap 时才等于基因组坐标
    if genes:
        ref_idx = None
        if ref_name is not None:
            ref_idx = ids.index(ref_name) if ref_name in ids else None
            if ref_idx is None:
                print(f"[popgen_figures] 警告: ref_name='{ref_name}' 不在序列中, "
                      f"基因区标注已跳过", file=sys.stderr)
        else:
            ref_idx = _find_reference_index(ids)
            if ref_idx is None:
                print("[popgen_figures] 警告: 无法唯一确定参考序列 (命名约定未命中或"
                      "命中多条), 基因区标注已跳过 —— 避免把参考坐标画到比对列轴上。"
                      " 可传 ref_name= 显式指定", file=sys.stderr)

        if ref_idx is not None:
            _pos2col = _ref_position_to_column_map(seqs[ref_idx])
            if len(_pos2col) == L:
                _xlabel = 'Genome position (nt)'
            _missing = 0
            ymax = ax.get_ylim()[1]
            for g, start, end, color in genes:
                c0 = _pos2col.get(int(start))
                c1 = _pos2col.get(int(end))
                if c0 is None or c1 is None:
                    _missing += 1
                    continue
                ax.axvspan(c0 + 1, c1 + 1, color=color, alpha=0.12)
                ax.text((c0 + c1) / 2 + 1, ymax * 0.95, g, ha='center', fontsize=7,
                        color=color, fontweight='bold')
            if _missing:
                print(f"[popgen_figures] 警告: {_missing}/{len(genes)} 个基因区坐标超出"
                      f"参考序列长度, 未标注", file=sys.stderr)
    ax.set_xlabel(_xlabel)
    if title:
        ax.set_title(title, fontsize=10)
    fig.tight_layout()
    for fmt in ['svg', 'pdf', 'png']:
        fig.savefig(f'{output}.{fmt}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    return ent


def plot_gene_variation(gene_stats, output, title=''):
    """每基因变异密度/π 条形图"""
    genes = gene_stats['genes']
    names = [g['gene'] for g in genes]
    per_kb = [g['per_kb'] for g in genes]
    pis = [g['pi'] for g in genes]
    high = [g['high_freq'] for g in genes]

    fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.2))
    x = np.arange(len(names))
    colors = [GENE_COLORS[i % len(GENE_COLORS)] for i in range(len(names))]

    # a: 每kb变异密度
    axes[0].bar(x, per_kb, color=colors, alpha=0.85, edgecolor='white', lw=0.5)
    axes[0].set_xticks(x); axes[0].set_xticklabels(names, fontsize=7)
    axes[0].set_ylabel('Variable sites / kb')
    axes[0].set_title('a  Variation density', fontsize=8, loc='left')

    # b: π per site
    axes[1].bar(x, [p * 100 for p in pis], color=colors, alpha=0.85, edgecolor='white', lw=0.5)
    axes[1].set_xticks(x); axes[1].set_xticklabels(names, fontsize=7)
    axes[1].set_ylabel('π (%)')
    axes[1].set_title('b  Nucleotide diversity', fontsize=8, loc='left')

    # c: 高频变异位点
    axes[2].bar(x, high, color=colors, alpha=0.85, edgecolor='white', lw=0.5)
    axes[2].set_xticks(x); axes[2].set_xticklabels(names, fontsize=7)
    axes[2].set_ylabel('High-frequency sites (≥5)')
    axes[2].set_title('c  High-frequency variants', fontsize=8, loc='left')

    if title:
        fig.suptitle(title, fontsize=10, y=1.02)
    fig.tight_layout()
    for fmt in ['svg', 'pdf', 'png']:
        fig.savefig(f'{output}.{fmt}', dpi=300, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    import argparse
    import sys as _sys
    _here = os.path.dirname(os.path.abspath(__file__))
    _sys.path.insert(0, os.path.dirname(_here))
    _sys.path.insert(0, _here)
    from utils.dataset_config import dataset as get_dataset, known_dataset_names
    from popgen_analysis import PopGenAnalyzer

    ap = argparse.ArgumentParser(description='变异热点绘图 (配置驱动, datasets.yaml)')
    ap.add_argument('--virus', default=None, help='已知数据集名; 缺省遍历全部')
    args = ap.parse_args()

    names = [args.virus.upper()] if args.virus else known_dataset_names()
    for name in names:
        d = get_dataset(name)
        if d is None:
            print(f"跳过未知数据集: {name}")
            continue
        work = f'{d["dir"]}/{d["work"]}'
        fasta = f'{d["dir"]}/{d["fasta"]}'
        meta = f'{d["dir"]}/{d["metadata"]}'
        out_sub = d.get('popgen_out') or 'figures_main'
        genes = d.get('genes') or None
        print(f"=== {name} ===")
        pg = PopGenAnalyzer(fasta, meta)
        out_base = f'{d["dir"]}/{out_sub}/{name.lower()}_entropy_profile'
        plot_entropy_profile(pg, out_base, genes=genes,
                             title=f'{name} nucleotide diversity across the genome')
        print(f"  熵谱 → {out_base}.{{svg,pdf,png}}")
        # 基因变异图 (仅在配置 genes_csv 且文件存在时)
        gcsv = d.get('genes_csv')
        if gcsv:
            gcsv_path = gcsv if os.path.isabs(gcsv) else f'{d["dir"]}/{gcsv}'
            if os.path.exists(gcsv_path):
                gm = pg.gene_mapping(genes_csv=gcsv_path)
                gv_out = f'{d["dir"]}/{out_sub}/{name.lower()}_gene_variation'
                plot_gene_variation(gm, gv_out, title=f'{name} variation by gene')
                print(f"  基因变异 → {gv_out}.{{svg,pdf,png}}")
            else:
                print(f"  (跳过 gene_variation: 未找到 {gcsv_path})")
