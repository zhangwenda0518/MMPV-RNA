#!/usr/bin/env python3
"""
summarize_assemblies.py — 汇总 3_Virus_assemblies_final 中每个病毒的组装统计
=====================================================================
扫描所有病毒/样本的 Global_Evolution_Stats.tsv (或 13.Evolution_Stats.tsv)，
按病毒汇总最终组装结果。

输出:
  assembly_summary.tsv       每个样本的最终组装统计
  assembly_summary_by_virus.tsv  按病毒汇总

用法:
  python summarize_assemblies.py -d 3_Virus_assemblies_final/ -o assembly_summary/
"""

import argparse
import os
import sys
import re
from pathlib import Path
from collections import defaultdict
import csv

# RefSeq/GenBank 风格 accession 后缀 (如 NC_002634.1); 目录名以 accession 结尾
_ACC_RE = re.compile(r'([A-Z]{1,2}_[0-9]{5,9}\.[0-9]+)$')


def _int_or_none(value):
    """数值字段容错解析: 空串/非法值返回 None (由调用方 WARN 并跳过)。"""
    s = str(value).strip() if value is not None else ""
    if not s:
        return None
    try:
        return int(float(s))
    except (ValueError, OverflowError):
        return None


def parse_evolution_stats(tsv_path):
    """解析 Global_Evolution_Stats.tsv / 13.Evolution_Stats.tsv"""
    rows = []
    try:
        with open(tsv_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                rows.append(row)
    except Exception as e:
        print(f"  警告: 无法解析 {tsv_path}: {e}")
    return rows


def extract_virus_from_dirname(dirname):
    """从目录名提取 virus accession。

    目录命名: {virus_name}_{accession}。accession 为 RefSeq 风格时
    (如 NC_002634.1) 仅用 rsplit('_', 1) 会截成 '002634.1',
    故优先按 accession 正则取末尾; 取不到再回退旧的 rsplit 行为。
    """
    m = _ACC_RE.search(str(dirname))
    if m:
        return m.group(1)
    # 例如: Grapevine_associated_RNA_virus_4_MW648525.1 → MW648525.1
    parts = str(dirname).rsplit('_', 1)
    if len(parts) == 2 and '.' in parts[1]:
        return parts[1]
    return dirname


def main():
    parser = argparse.ArgumentParser(description="汇总病毒组装统计")
    parser.add_argument('-d', '--assembly-dir', required=True, help='3_Virus_assemblies_final 目录')
    parser.add_argument('-o', '--outdir', default='assembly_summary', help='输出目录')
    parser.add_argument('--min-non-n', type=float, default=0, help='最小 Non_N_Bases 阈值 (过滤)')
    args = parser.parse_args()

    asm_dir = Path(args.assembly_dir)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 收集所有 stats 文件
    stats_files = list(asm_dir.rglob('Global_Evolution_Stats.tsv'))
    if not stats_files:
        stats_files = list(asm_dir.rglob('13.Evolution_Stats.tsv'))
    if not stats_files:
        # 回退：搜索任何 *Evolution_Stats* 或 *evolution*stats*
        stats_files = list(asm_dir.rglob('*volution*tats*'))
        stats_files += list(asm_dir.rglob('*volution*tats*'))

    print(f"找到 {len(stats_files)} 个统计文件")

    # 解析所有数据
    all_samples = []  # [{sample, virus, virus_dir, step, total_len, non_n, n_count, n50, ...}]
    virus_stats = defaultdict(lambda: {
        'samples': 0, 'total_assembled': 0, 'complete': 0,
        'min_len': float('inf'), 'max_len': 0, 'avg_len': 0,
        'min_non_n': float('inf'), 'max_non_n': 0, 'avg_non_n': 0,
        'lengths': [], 'non_n_lengths': [],
        'sample_list': []
    })

    for tsv_path in sorted(stats_files):
        # 路径结构: virus_dir/sample_dir/sample_subdir/13.Evolution_Stats.tsv
        # 或 virus_dir/sample_dir/Global_Evolution_Stats.tsv
        parts = tsv_path.parts
        rel = tsv_path.relative_to(asm_dir)

        # 推断 virus_dir (第一级) 和 sample (第二级)
        if len(rel.parts) >= 1:
            virus_dir = rel.parts[0]
        else:
            virus_dir = 'Unknown'

        # 从路径推断 sample name
        if len(rel.parts) >= 2:
            sample_dir = rel.parts[1]
        else:
            sample_dir = tsv_path.parent.name

        rows = parse_evolution_stats(tsv_path)
        if not rows:
            continue

        # 取最后一步（通常叫 Ultimate_Result 或最后一行）
        final_row = None
        for row in reversed(rows):
            step = row.get('Step', '')
            if 'Ultimate' in step or 'Result' in step or 'final' in step.lower():
                final_row = row
                break
        if not final_row:
            final_row = rows[-1]  # 取最后一行

        sample = final_row.get('Sample', sample_dir)
        step = final_row.get('Step', '')
        # 数值字段容错: 半成品/截断的 stats 文件缺字段时 WARN 并跳过该样本,
        # 不再让单个文件的 ValueError 终止整个汇总
        total_len = _int_or_none(final_row.get('Total_Length'))
        non_n = _int_or_none(final_row.get('Non_N_Bases'))
        n_count = _int_or_none(final_row.get('N_Count'))
        n50 = _int_or_none(final_row.get('Contig_N50'))
        if None in (total_len, non_n, n_count, n50):
            print(f"  警告: {tsv_path} 最终行数值字段缺失/非法 "
                  f"(Total_Length={final_row.get('Total_Length')!r}, "
                  f"Non_N_Bases={final_row.get('Non_N_Bases')!r}, "
                  f"N_Count={final_row.get('N_Count')!r}, "
                  f"Contig_N50={final_row.get('Contig_N50')!r}), 跳过该样本")
            continue

        # 推断 virus accession
        virus_acc = extract_virus_from_dirname(virus_dir)

        if non_n < args.min_non_n:
            continue

        all_samples.append({
            'virus': virus_acc,
            'virus_dir': virus_dir,
            'sample': sample,
            'total_len': total_len,
            'non_n': non_n,
            'n_count': n_count,
            'n50': n50,
            'step': step,
        })

        vs = virus_stats[virus_acc]
        vs['samples'] += 1
        vs['lengths'].append(total_len)
        vs['non_n_lengths'].append(non_n)
        vs['sample_list'].append(sample)
        if n_count == 0:
            vs['complete'] += 1
        vs['total_assembled'] += 1

    # 计算病毒汇总
    for acc, vs in virus_stats.items():
        if vs['lengths']:
            vs['min_len'] = min(vs['lengths'])
            vs['max_len'] = max(vs['lengths'])
            vs['avg_len'] = sum(vs['lengths']) / len(vs['lengths'])
        else:
            vs['min_len'] = vs['max_len'] = vs['avg_len'] = 0
        if vs['non_n_lengths']:
            vs['min_non_n'] = min(vs['non_n_lengths'])
            vs['max_non_n'] = max(vs['non_n_lengths'])
            vs['avg_non_n'] = sum(vs['non_n_lengths']) / len(vs['non_n_lengths'])
        else:
            vs['min_non_n'] = vs['max_non_n'] = vs['avg_non_n'] = 0

    # ── 输出 1: 每个样本的最终组装统计 ──
    sample_out = out_dir / 'assembly_summary.tsv'
    with open(sample_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, delimiter='\t',
            fieldnames=['Virus', 'Virus_Dir', 'Sample', 'Step', 'Total_Length', 'Non_N_Bases', 'N_Count', 'N50'])
        writer.writeheader()
        for s in sorted(all_samples, key=lambda x: (x['virus'], x['sample'])):
            writer.writerow({
                'Virus': s['virus'], 'Virus_Dir': s['virus_dir'], 'Sample': s['sample'],
                'Step': s['step'], 'Total_Length': s['total_len'],
                'Non_N_Bases': s['non_n'], 'N_Count': s['n_count'], 'N50': s['n50'],
            })
    print(f"\n样本明细: {sample_out} ({len(all_samples)} 条)")

    # ── 输出 2: 按病毒汇总 ──
    virus_out = out_dir / 'assembly_summary_by_virus.tsv'
    with open(virus_out, 'w', newline='') as f:
        writer = csv.DictWriter(f, delimiter='\t',
            fieldnames=['Virus', 'Samples', 'Complete(0N)', 'Min_Len', 'Max_Len', 'Avg_Len',
                        'Min_NonN', 'Max_NonN', 'Avg_NonN', 'Sample_List'])
        writer.writeheader()
        for acc in sorted(virus_stats.keys()):
            vs = virus_stats[acc]
            writer.writerow({
                'Virus': acc, 'Samples': vs['samples'], 'Complete(0N)': vs['complete'],
                'Min_Len': vs['min_len'], 'Max_Len': vs['max_len'], 'Avg_Len': f"{vs['avg_len']:.0f}",
                'Min_NonN': vs['min_non_n'], 'Max_NonN': vs['max_non_n'], 'Avg_NonN': f"{vs['avg_non_n']:.0f}",
                'Sample_List': ', '.join(sorted(vs['sample_list'])),
            })
    print(f"\n病毒汇总: {virus_out} ({len(virus_stats)} 病毒)")

    # ── 打印摘要 ──
    print(f"\n{'Virus':<20} {'Samples':>8} {'Complete':>9} {'Avg_Len':>8} {'Avg_NonN':>9}")
    print('-' * 58)
    for acc in sorted(virus_stats.keys()):
        vs = virus_stats[acc]
        print(f"{acc:<20} {vs['samples']:>8} {vs['complete']:>9} {vs['avg_len']:>8.0f} {vs['avg_non_n']:>9.0f}")


if __name__ == '__main__':
    main()
