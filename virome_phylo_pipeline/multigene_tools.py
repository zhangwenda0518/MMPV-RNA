#!/usr/bin/env python3
"""
multigene_tools.py — 多基因工具 (序列拼接器 + TreeMerger)

sequence_concatenator:
  多个基因比对 (同样本集) → 串联 supermatrix (FASTA + 分区文件)
  参考: YR-MPE Sequence Concatenator 逻辑

tree_merger:
  多个树文件 → 合并成一个 NEXUS trees 块 (供 TreeAnnotator/MCC 查看)

用法:
  python -m multigene_tools concat --alns aln1.fa aln2.fa --out supermatrix.fa
  python -m multigene_tools merge --trees t1.nwk t2.nwk --out merged.nex
"""

import argparse
import os
import sys
from collections import OrderedDict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def concatenate_alignments(aln_paths, out_path, partition_path=None, log=None):
    """串联多个基因比对 → supermatrix FASTA (+ 分区)

    设计取舍（有意与 YR-MPE Sequence Concatenator 不同，非疏忽）:
    - 样本取交集而非并集：病毒 supermatrix 宁缺毋滥，避免大量 '?' 填充列
      人为拉低下游分辨率；如需并集策略（YR-MPE 默认）需另行实现
    - 键用 r.id（仅 FASTA 头第一段）而非 seq.description：避免带注释的
      头在拼接时产生同名不同键的静默丢失
    - 缺失填充用 'N' 而非 '?'：与本项目其他模块的缺失约定一致
    """
    def emit(msg):
        if log:
            log(msg)

    gene_seqs = []   # [{sample: SeqRecord}]
    for p in aln_paths:
        recs_all = {r.id: str(r.seq).upper() for r in SeqIO.parse(p, 'fasta')}
        # 同文件内序列等长校验 (YR-MPE 有此校验, 缺失会静默产出错位 supermatrix)
        _lens = {len(s) for s in recs_all.values()}
        if len(_lens) > 1:
            raise SystemExit(f'[concat] {Path(p).name}: 序列长度不等 {_lens}，'
                             f'请先检查比对是否完整')
        recs = recs_all
        gene_seqs.append(recs)
        emit(f'[concat] {Path(p).name}: {len(recs)} 序列')

    if not gene_seqs:
        raise SystemExit('无输入比对')

    # 样本交集
    common = set(gene_seqs[0].keys())
    for gs in gene_seqs[1:]:
        common &= set(gs.keys())
    common = sorted(common)
    emit(f'[concat] 共同样本: {len(common)}')

    if not common:
        raise SystemExit('无共同样本, 无法拼接')

    # 各基因长度
    lengths = [len(next(iter(gs.values()))) for gs in gene_seqs]
    total_len = sum(lengths)
    emit(f'[concat] 各基因长度: {lengths}, 总长: {total_len}')

    # 拼接 (缺失基因用 N)
    records = []
    for sample in common:
        concat = []
        for idx, gs in enumerate(gene_seqs):
            concat.append(gs.get(sample, 'N' * lengths[idx]))
        records.append(SeqRecord(Seq(''.join(concat)), id=sample, description=''))

    with open(out_path, 'w', encoding='utf-8') as f:
        SeqIO.write(records, f, 'fasta')
    emit(f'[concat] → {out_path}')

    # 分区文件 (RAxML/NEXUS 风格)
    if partition_path:
        parts = []
        start = 1
        for name, ln in zip([Path(p).stem for p in aln_paths], lengths):
            parts.append(f'DNA, {name} = {start}-{start + ln - 1}')
            start += ln
        with open(partition_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(parts) + '\n')
        emit(f'[concat] 分区 → {partition_path}')

    return {'n_samples': len(common), 'n_genes': len(aln_paths), 'total_len': total_len,
            'fasta': out_path, 'partition': partition_path}


def merge_trees(tree_paths, out_path, log=None):
    """合并多个树 → NEXUS trees 块"""
    def emit(msg):
        if log:
            log(msg)

    trees = []
    for p in tree_paths:
        content = Path(p).read_text(encoding='utf-8', errors='replace').strip()
        # 提取 newick (跳过注释/nexus 头)
        for line in content.splitlines():
            line = line.strip()
            # 2026-09-15 (审查 P2-20) 修复: 旧写法 `A or B and C` 因 and 优先级更高
            # 等价于 `A or (B and C)` —— 只要以 '(' 开头就无条件收下, 哪怕它不是
            # 完整 newick (没以 ';' 收尾); 而 ':' 分支又额外要求结尾 ';'。两条
            # 判据不对称, 容易把 nexus 头/注释行当树收进去。统一为:
            # 必须以 ';' 收尾, 且 (以 '(' 开头 或 前 50 字符内出现 ':')。
            if line.endswith(';') and (line.startswith('(') or ':' in line[:50]):
                trees.append(line)
                break
        else:
            # 整文件可能是 newick
            if content.endswith(';'):
                trees.append(content)

    if not trees:
        raise SystemExit('无有效树')

    emit(f'[merge] 合并 {len(trees)} 棵树')
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('#NEXUS\nBEGIN TREES;\n')
        for i, t in enumerate(trees):
            f.write(f'  TREE tree_{i+1} = {t}\n')
        f.write('END;\n')
    emit(f'[merge] → {out_path}')
    return {'n_trees': len(trees), 'nexus': out_path}


def main():
    ap = argparse.ArgumentParser(description='多基因工具: 拼接器 + TreeMerger')
    sub = ap.add_subparsers(dest='cmd', required=True)

    p1 = sub.add_parser('concat', help='序列拼接 (多基因 → supermatrix)')
    p1.add_argument('--alns', nargs='+', required=True, help='基因比对 FASTA (同样本集)')
    p1.add_argument('--out', required=True, help='输出 supermatrix FASTA')
    p1.add_argument('--partition', default=None, help='分区文件输出')

    p2 = sub.add_parser('merge', help='TreeMerger (多树 → NEXUS)')
    p2.add_argument('--trees', nargs='+', required=True, help='树文件 (newick)')
    p2.add_argument('--out', required=True, help='输出 NEXUS')

    args = ap.parse_args()
    if args.cmd == 'concat':
        r = concatenate_alignments(args.alns, args.out, args.partition, log=print)
        print(f"✓ 拼接完成: {r['n_samples']} 样本 × {r['n_genes']} 基因, 总长 {r['total_len']}")
    elif args.cmd == 'merge':
        r = merge_trees(args.trees, args.out, log=print)
        print(f"✓ 合并完成: {r['n_trees']} 棵树 → {r['nexus']}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
