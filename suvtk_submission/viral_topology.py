#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
viral_topology.py — 病毒基因组拓扑判断 (circular / linear) v1.0
==============================================================

结合 suvtk taxonomy 分类规则 + 序列末端重复检测, 判断每个病毒基因组是环形还是线性。

方法:
  1. 基于 ICTV 分类的已知规则 (某些科天然是环形的)
  2. 基于序列的直接末端重复 (DTR) 检测 (参考 Cenote-Taker3 terminal_repeats.py)
  3. 规则优先, 序列检测兜底

用法:
  python viral_topology.py \\
      --taxonomy suvtk.taxonomy_output/taxonomy.tsv \\
      --fasta sequences.fasta \\
      -o topology.tsv
"""

import argparse
import os
import sys
import re
from pathlib import Path
from collections import defaultdict

from Bio import SeqIO
from tqdm import tqdm

# ══════════════════════════════════════════════════════════════
# 已知环形/线性的病毒分类规则 (基于 ICTV)
# ══════════════════════════════════════════════════════════════

# 环状病毒科/属 (ssDNA, 某些 dsDNA, 类病毒)
CIRCULAR_TAXA = {
    # ssDNA 环状
    'geminiviridae', 'nanoviridae', 'circoviridae', 'anelloviridae',
    'bidnaviridae', 'parvoviridae',  # Parvoviridae 线性的, 排除
    'genomoviridae', 'smacoviridae', 'redondoviridae',
    # dsDNA 环状 (部分)
    'polyomaviridae', 'papillomaviridae',
    # 类病毒
    'pospiviroidae', 'avsunviroidae',
    # 环状 RNA 病毒 (部分)
    'deltavirus',
    # 噬菌体
    'microviridae', 'inoviridae', 'plasmaviridae',
}

# 明确线性的科
LINEAR_TAXA = {
    'rhabdoviridae', 'paramyxoviridae', 'filoviridae', 'bornaviridae',
    'pneumoviridae', 'orthomyxoviridae', 'bunyaviridae', 'arenaviridae',
    'coronaviridae', 'flaviviridae', 'togaviridae', 'picornaviridae',
    'potyviridae', 'tombusviridae', 'virgaviridae', 'bromoviridae',
    'closteroviridae', 'luteoviridae', 'secoviridae', 'tymoviridae',
    'alphaflexiviridae', 'betaflexiviridae', 'gammaflexiviridae',
    'partitiviridae', 'chrysoviridae', 'totiviridae', 'reoviridae',
    'phenuiviridae', 'tospoviridae', 'fimoviridae',
    'narnaviridae', 'mitoviridae', 'botourmiaviridae',
    'ourmiavirus', 'tobravirus', 'carlavirus', 'potexvirus',
    'vitivirus', 'foveavirus', 'capillovirus', 'trichovirus',
    'ampelovirus', 'badnavirus', 'caulimovirus', 'soymovirus',
    'prunevirus', 'solemoviridae', 'nodaviridae',
    'ilarvirus', 'alphavirus',
}


def taxonomy_to_topology(taxonomy_str):
    """根据 taxonomy 字符串推断拓扑"""
    tax_lower = taxonomy_str.lower()
    for circ in CIRCULAR_TAXA:
        if circ in tax_lower:
            return 'circular'
    for lin in LINEAR_TAXA:
        if lin in tax_lower:
            return 'linear'
    return 'undetermined'


# ══════════════════════════════════════════════════════════════
# 序列末端重复检测 (DTR) — 参考 Cenote-Taker3
# ══════════════════════════════════════════════════════════════

def find_dtr(seq, min_len=20, max_len=200, min_identity=0.9):
    """检测序列两端是否存在直接末端重复 (DTR)

    如果 DTR 长度 >= min_len 且 identity >= min_identity, 判定为环形。

    返回: (is_circular, dtr_length, dtr_identity)
    """
    seq_str = str(seq).upper()
    seq_len = len(seq_str)
    if seq_len < min_len * 2:
        return False, 0, 0.0

    best_len, best_id = 0, 0.0

    for dtr_len in range(min_len, min(max_len, seq_len // 4)):
        left = seq_str[:dtr_len]
        right = seq_str[seq_len - dtr_len:]
        matches = sum(1 for a, b in zip(left, right) if a == b)
        identity = matches / dtr_len

        if identity > best_id:
            best_id = identity
            best_len = dtr_len

        if identity >= min_identity and dtr_len >= min_len:
            return True, dtr_len, identity

    return False, best_len, best_id


def find_itr(seq, min_len=20, max_len=200):
    """检测倒转末端重复 (ITR) — 回文结构

    ITR 也暗示环形起源, 但不如 DTR 可靠。
    """
    from Bio.Seq import Seq
    seq_str = str(seq).upper()
    seq_len = len(seq_str)
    if seq_len < min_len * 2:
        return False, 0

    for dtr_len in range(min_len, min(max_len, seq_len // 4)):
        left = seq_str[:dtr_len]
        right_rev = str(Seq(seq_str[seq_len - dtr_len:]).reverse_complement())
        matches = sum(1 for a, b in zip(left, right_rev) if a == b)
        if matches / dtr_len >= 0.85:
            return True, dtr_len

    return False, 0


# ══════════════════════════════════════════════════════════════
# 主逻辑
# ══════════════════════════════════════════════════════════════

def determine_topology(taxonomy_file, fasta_file, output_file,
                       min_dtr=20, min_identity=0.9):
    """综合判断每条序列的基因组拓扑"""

    # 加载 taxonomy
    tax_map = {}
    with open(taxonomy_file, 'r') as f:
        header = f.readline().strip().split('\t')
        c_idx = 0
        t_idx = 1
        for i, h in enumerate(header):
            if h.lower() in ('contig', 'seq_id', 'sequence_id'):
                c_idx = i
            elif h.lower() in ('taxonomy', 'tax'):
                t_idx = i
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) > max(c_idx, t_idx):
                tax_map[cols[c_idx]] = cols[t_idx]

    # 加载序列
    records = list(SeqIO.parse(fasta_file, 'fasta'))

    results = []
    rule_count, dtr_count, undet_count = 0, 0, 0

    for rec in tqdm(records, desc="检测拓扑结构"):
        taxonomy = tax_map.get(rec.id, 'unclassified viruses')

        # 策略 1: 分类规则
        rule_topology = taxonomy_to_topology(taxonomy)

        # 策略 2: DTR 检测
        is_circ_dtr, dtr_len, dtr_id = find_dtr(rec.seq, min_dtr, min_identity)
        is_circ_itr, itr_len = find_itr(rec.seq, min_dtr)

        # 综合决策
        if rule_topology == 'circular':
            final = 'circular'
            evidence = f'taxonomy ({taxonomy.split(";")[-1] if ";" in taxonomy else taxonomy})'
            rule_count += 1
        elif rule_topology == 'linear':
            final = 'linear'
            evidence = f'taxonomy ({taxonomy.split(";")[-1] if ";" in taxonomy else taxonomy})'
            rule_count += 1
        elif is_circ_dtr:
            final = 'circular'
            evidence = f'DTR found: len={dtr_len}, identity={dtr_id:.2%}'
            dtr_count += 1
        elif is_circ_itr:
            final = 'circular'
            evidence = f'ITR found: len={itr_len}'
            dtr_count += 1
        else:
            # 默认: RNA 病毒大部分是线性的
            if 'RNA' in taxonomy.upper() or 'ssRNA' in taxonomy or 'dsRNA' in taxonomy:
                final = 'linear'
                evidence = 'assumed (RNA virus, no DTR detected)'
            else:
                final = 'linear'
                evidence = 'assumed (no evidence for circular)'
            undet_count += 1

        results.append({
            'contig': rec.id,
            'taxonomy': taxonomy,
            'rule_topology': rule_topology,
            'final_topology': final,
            'seq_length': len(rec.seq),
            'dtr_detected': is_circ_dtr,
            'dtr_length': dtr_len,
            'dtr_identity': f"{dtr_id:.3f}" if is_circ_dtr else '0',
            'evidence': evidence,
        })

    # 写结果
    import pandas as pd
    df = pd.DataFrame(results)
    df.to_csv(output_file, sep='\t', index=False)

    # 统计
    n_circ = (df['final_topology'] == 'circular').sum()
    n_linear = (df['final_topology'] == 'linear').sum()

    print(f"\n{'='*50}")
    print(f"  基因组拓扑检测结果")
    print(f"  {'='*50}")
    print(f"  总序列数:    {len(results)}")
    print(f"  环状:        {n_circ}")
    print(f"  线性:        {n_linear}")
    print(f"  {'='*50}")
    print(f"  分类规则确定: {rule_count}")
    print(f"  DTR/ITR 检测: {dtr_count}")
    print(f"  默认推断:     {undet_count}")
    print(f"  {'='*50}")
    print(f"  输出: {output_file}")

    return df


def main():
    parser = argparse.ArgumentParser(
        description="viral_topology.py — 病毒基因组拓扑判断 (circular/linear)")
    parser.add_argument('--taxonomy', required=True, help='suvtk taxonomy.tsv')
    parser.add_argument('--fasta', required=True, help='输入 FASTA (核酸序列)')
    parser.add_argument('-o', '--output', default='topology.tsv', help='输出文件')
    parser.add_argument('--min-dtr', type=int, default=20, help='DTR 最短长度 (默认 20)')
    parser.add_argument('--min-identity', type=float, default=0.85, help='DTR 最低 identity (默认 0.85)')
    args = parser.parse_args()

    if not os.path.exists(args.taxonomy):
        print(f"错误: taxonomy 文件不存在: {args.taxonomy}")
        sys.exit(1)
    if not os.path.exists(args.fasta):
        print(f"错误: FASTA 文件不存在: {args.fasta}")
        sys.exit(1)

    determine_topology(args.taxonomy, args.fasta, args.output,
                       args.min_dtr, args.min_identity)


if __name__ == '__main__':
    main()
