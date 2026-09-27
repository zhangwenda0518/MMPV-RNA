#!/usr/bin/env python3
"""
distance_matrix.py — 遗传距离矩阵 + NJ 树 + 可视化

p-distance (未校正距离) 计算 + 距离矩阵 + Neighbor-Joining 树 + 热图。
参考: YR-MPE distance_utils (p-distance + IUPAC 简并碱基处理)

用法:
  python -m utils.distance_matrix --fasta mafft.aln.fasta --outdir dist_out
"""

import argparse
import os
import sys
from collections import OrderedDict

import numpy as np
from Bio import SeqIO

# IUPAC 简并碱基映射
IUPAC_DNA_MAP = {
    'A': 'A', 'T': 'T', 'C': 'C', 'G': 'G',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}


def _expand(base):
    b = base.upper()
    return list(IUPAC_DNA_MAP.get(b, b if b in 'ATCG' else ''))


def _base_diff(b1, b2):
    bases1 = _expand(b1)
    bases2 = _expand(b2)
    if not bases1 or not bases2:
        return 1.0
    return 0.0 if any(a == b for a in bases1 for b in bases2) else 1.0


def p_distance(seq1, seq2, gap_treatment='pairwise'):
    """两序列未校正 p-distance。

    gap_treatment: 仅 'pairwise' 有效（逐对跳过 gap 位点）。
    complete deletion 在 distance_matrix() 层实现（全比对剔列后统一计算），
    本函数收到 'complete' 时按 pairwise 处理（列已在上游剔除，无 gap 可遇）。

    注: 与 YR-MPE 的已知差异——'?' 等非 IUPAC 字符在 YR-MPE 中计入分母且算差异，
    本实现跳过该位点（不计分母）。N 双方一致（展开 ACGT 视为可匹配）。
    """
    min_len = min(len(seq1), len(seq2))
    diff = comp = 0
    for i in range(min_len):
        b1, b2 = seq1[i].upper(), seq2[i].upper()
        if gap_treatment == 'complete':
            # complete deletion: 任一方为 gap 则整个位点剔除，其余字符照常计
            if b1 == '-' or b2 == '-':
                continue
        else:
            # pairwise deletion: 跳过 gap 位点
            if b1 == '-' or b2 == '-':
                continue
        if b1 not in IUPAC_DNA_MAP or b2 not in IUPAC_DNA_MAP:
            continue
        comp += 1
        if _base_diff(b1, b2) > 0:
            diff += 1
    return diff / comp if comp > 0 else 0.0


def _complete_deletion_columns(sequences):
    """YR-MPE 口径: 全比对中只要任何序列该位为 gap（'-'）则剔除该列，
    其余列（含 N/简并碱基）保留。返回保留列索引列表。"""
    seqs = [str(s.seq).upper() for s in sequences]
    L = min(len(s) for s in seqs)
    return [i for i in range(L) if all(s[i] != '-' for s in seqs)]


def distance_matrix(sequences, gap_treatment='pairwise'):
    """n×n p-distance 矩阵。

    gap_treatment='complete' 时先在全比对层面剔除含 gap 列（complete deletion），
    所有配对都在同一组保留列上计算（与 YR-MPE complete_deletion_compare 同口径），
    而非逐对跳过。
    """
    n = len(sequences)
    mat = np.zeros((n, n))
    seqs = [str(s.seq) for s in sequences]
    if gap_treatment == 'complete':
        keep = _complete_deletion_columns(sequences)
        sub = ["".join(s[i] for i in keep) for s in seqs]
        for i in range(n):
            for j in range(i + 1, n):
                d = p_distance(sub[i], sub[j], gap_treatment='pairwise')
                mat[i, j] = mat[j, i] = d
        return mat
    for i in range(n):
        for j in range(i + 1, n):
            d = p_distance(seqs[i], seqs[j], gap_treatment)
            mat[i, j] = mat[j, i] = d
    return mat


def nj_tree(sequences, gap_treatment='pairwise'):
    """Neighbor-Joining 树 (Bio.Phylo DistanceTreeConstructor)"""
    from Bio.Phylo.TreeConstruction import DistanceMatrix as BDM, DistanceTreeConstructor
    mat = distance_matrix(sequences, gap_treatment)
    names = [s.id for s in sequences]
    n = len(names)
    # Bio.Phylo DistanceMatrix 需含对角线的下三角格式 (第 i 行含 i+1 个元素)
    lower = [mat[i, :i + 1].tolist() for i in range(n)]
    dm = BDM(names, lower)
    ctor = DistanceTreeConstructor()
    return ctor.nj(dm)


def plot_matrix(mat, names, out_path):
    """距离矩阵热图"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(max(6, len(names) * 0.5), max(6, len(names) * 0.5)))
    im = ax.imshow(mat, cmap='viridis')
    ax.set_xticks(range(len(names)))
    ax.set_yticks(range(len(names)))
    ax.set_xticklabels(names, rotation=90, fontsize=6)
    ax.set_yticklabels(names, fontsize=6)
    ax.set_title('p-distance matrix')
    fig.colorbar(im, ax=ax, label='p-distance')
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return out_path


def plot_tree(tree, out_path):
    """NJ 树图"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from Bio import Phylo
    fig, ax = plt.subplots(figsize=(10, max(4, len(tree.get_terminals()) * 0.3)))
    Phylo.draw(tree, axes=ax, do_show=False)
    ax.set_title('Neighbor-Joining tree (p-distance)')
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return out_path


def run_distance_analysis(fasta_path, outdir='dist_out', gap_treatment='pairwise', log=None):
    def emit(m):
        if log:
            log(m)
    os.makedirs(outdir, exist_ok=True)
    result = {'success': False}
    if not os.path.exists(fasta_path):
        result['error'] = f'文件不存在: {fasta_path}'
        return result

    sequences = list(SeqIO.parse(fasta_path, 'fasta'))
    if len(sequences) < 2:
        result['error'] = '序列数 < 2'
        return result
    emit(f'[distance] {len(sequences)} 条序列')

    names = [s.id for s in sequences]
    mat = distance_matrix(sequences, gap_treatment)

    # CSV
    csv_path = os.path.join(outdir, 'distance_matrix.csv')
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        f.write(',' + ','.join(names) + '\n')
        for i, name in enumerate(names):
            f.write(name + ',' + ','.join(f'{mat[i, j]:.6f}' for j in range(len(names))) + '\n')
    emit(f'[distance] 矩阵 → {csv_path}')

    # NJ 树
    tree = nj_tree(sequences, gap_treatment)
    tree_path = os.path.join(outdir, 'nj_tree.nwk')
    from Bio import Phylo
    Phylo.write(tree, tree_path, 'newick')
    emit(f'[distance] NJ 树 → {tree_path}')

    # 图
    heat = plot_matrix(mat, names, os.path.join(outdir, 'distance_heatmap.png'))
    treefig = plot_tree(tree, os.path.join(outdir, 'nj_tree.png'))
    emit(f'[distance] 图 → {heat}, {treefig}')

    result.update({'success': True, 'matrix': csv_path, 'tree': tree_path,
                   'heatmap': heat, 'tree_fig': treefig,
                   'n_sequences': len(sequences),
                   'mean_distance': float(np.mean(mat[np.triu_indices(len(names), 1)]))})
    return result


def main():
    ap = argparse.ArgumentParser(description='遗传距离矩阵 + NJ 树')
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--outdir', default='dist_out')
    ap.add_argument('--gap', default='pairwise', choices=['pairwise', 'complete'])
    args = ap.parse_args()
    r = run_distance_analysis(args.fasta, args.outdir, args.gap, log=print)
    if r.get('success'):
        print(f"✓ 距离分析完成: 平均 p-distance={r['mean_distance']:.6f}")
        return 0
    print(f"✗ 失败: {r.get('error')}")
    return 1


if __name__ == '__main__':
    sys.exit(main())
