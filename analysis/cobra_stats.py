#!/usr/bin/env python3
"""
简版 COBRA 延伸统计 (无 CheckV 依赖)
用法: python cobra_stats.py -c 03_COBRA/ [-o cobra_summary.tsv]
产出:
  cobra_summary.tsv — 每样本延伸率/成功数/失败数
  cobra_contig_detail.tsv — 每条 contig 的延伸前后长度变化 + 孤儿末端
"""
import argparse, os, sys
from pathlib import Path
from collections import defaultdict

def parse_cobra_log(log_file):
    """解析 COBRA log → 返回 stats dict"""
    stats = {'total_queries': 0, 'extended_circular': 0, 'extended_partial': 0,
             'extended_failed': 0, 'orphan_end': 0, 'self_circular': 0,
             'extension_rate': 0.0}
    try:
        with open(log_file) as f:
            text = f.read()
        for line in text.split('\n'):
            s = line.strip()
            if s.startswith('# Total queries:'):
                stats['total_queries'] = int(s.split(':')[1].strip())
            elif 'Self_circular' in s:
                stats['self_circular'] = int(s.split(':')[1].strip().split()[0])
            elif 'Extended_circular' in s:
                v = s.split(':')[1].strip()
                stats['extended_circular'] = int(v.split()[0])
                if 'Unique:' in s:
                    stats['extended_circular_unique'] = int(s.split('Unique:')[1].strip().rstrip(')'))
            elif 'Extended_partial' in s:
                v = s.split(':')[1].strip()
                stats['extended_partial'] = int(v.split()[0])
                if 'Unique:' in s:
                    stats['extended_partial_unique'] = int(s.split('Unique:')[1].strip().rstrip(')'))
            elif 'Extended_failed' in s:
                stats['extended_failed'] = int(s.split(':')[1].strip())
            elif 'Orphan end' in s:
                stats['orphan_end'] = int(s.split(':')[1].strip())
        total = stats['extended_circular'] + stats['extended_partial']
        if stats['total_queries'] > 0:
            stats['extension_rate'] = round(total / stats['total_queries'] * 100, 1)
    except Exception as e:
        print(f"  [WARN] 解析 log 失败: {log_file}: {e}")
    return stats

def parse_fasta_lengths(fa_path):
    """返回 {contig_id: length}"""
    lens = {}; seq_id = None; seq_len = 0
    if not fa_path.is_file(): return lens
    for line in open(fa_path):
        l = line.strip()
        if l.startswith('>'):
            if seq_id: lens[seq_id] = seq_len
            seq_id = l[1:].split()[0]; seq_len = 0
        else: seq_len += len(l)
    if seq_id: lens[seq_id] = seq_len
    return lens

def analyze_cobra_dir(cobra_dir, out_dir):
    """分析 03_COBRA/ 目录"""
    cobra = Path(cobra_dir)
    if not cobra.is_dir():
        sys.exit(f"目录不存在: {cobra_dir}")

    sample_stats = []
    contig_details = []

    for sd in sorted(cobra.iterdir()):
        if not sd.is_dir(): continue
        sample = sd.name

        # 1. 查找 cobra.fa (延伸结果)
        cobra_fas = list(sd.rglob("*.cobra.fa"))
        # 2. 查找 virus.fa (原始输入, 用于计算延伸增量)
        virus_fas = list(sd.rglob("*virus*.fa"))
        virus_fas = [f for f in virus_fas if 'cobra' not in f.name]

        # 3. 查找 COBRA log
        log_files = list(sd.rglob("log"))
        # 过滤只取 COBRA 内部 log (在 *COBRA*/ 目录下的)
        cobra_logs = [f for f in log_files if 'COBRA' in str(f.parent.name)]
        if not cobra_logs:
            cobra_logs = [f for f in log_files if 'cobra' in str(f).lower()]

        row = {'Sample': sample, 'Has_Cobra': len(cobra_fas) > 0,
               'Has_Virus': len(virus_fas) > 0,
               'Cobra_FA': len(cobra_fas), 'Virus_FA': len(virus_fas),
               'Extension_Rate(%)': 0.0, 'Total_Queries': 0,
               'Extended_Circular': 0, 'Extended_Partial': 0,
               'Extended_Failed': 0, 'Orphan_End': 0, 'Self_Circular': 0}

        if cobra_logs:
            log_stats = parse_cobra_log(str(cobra_logs[0]))
            row.update({f'Total_Queries': log_stats['total_queries'],
                        f'Extended_Circular': log_stats['extended_circular'],
                        f'Extended_Partial': log_stats['extended_partial'],
                        f'Extended_Failed': log_stats['extended_failed'],
                        f'Orphan_End': log_stats['orphan_end'],
                        f'Self_Circular': log_stats['self_circular'],
                        f'Extension_Rate(%)': log_stats['extension_rate']})

        # 4. 对比 virus.fa vs cobra.fa → 计算每条 contig 的延伸 bp
        virus_lens = {}
        for vf in virus_fas:
            virus_lens.update(parse_fasta_lengths(vf))
        cobra_lens = {}
        for cf in cobra_fas:
            cobra_lens.update(parse_fasta_lengths(cf))

        n_extended = 0; n_new = 0; total_gain = 0
        for cid, clen in cobra_lens.items():
            # 查找匹配的原始 contig (去掉可能的后缀)
            orig_id = cid.replace('_extended', '').replace('.cobra', '')
            vlen = virus_lens.get(cid, virus_lens.get(orig_id, 0))
            gain = clen - vlen
            status = 'extended' if gain > 0 else ('same' if gain == 0 else 'new')
            if status == 'extended': n_extended += 1; total_gain += gain
            elif status == 'new': n_new += 1
            contig_details.append({
                'Sample': sample, 'Contig_ID': cid,
                'Virus_Len(bp)': vlen, 'Cobra_Len(bp)': clen,
                'Gain(bp)': max(gain, 0), 'Status': status,
            })

        row['Extended_Contigs'] = n_extended
        row['New_Contigs'] = n_new
        row['Total_Gain(bp)'] = total_gain
        # 孤儿末端比例
        if row['Total_Queries'] > 0:
            row['Orphan_Rate(%)'] = round(row['Orphan_End'] / row['Total_Queries'] * 100, 1)
        else:
            row['Orphan_Rate(%)'] = 0.0

        sample_stats.append(row)
        print(f"  {sample}: queries={row['Total_Queries']}, "
              f"extended={row['Extended_Contigs']}, "
              f"rate={row['Extension_Rate(%)']}%, "
              f"orphan={row['Orphan_End']}({row['Orphan_Rate(%)']}%)")

    # 写入 sample summary
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cols = ['Sample', 'Total_Queries', 'Extended_Circular', 'Extended_Partial',
            'Extended_Failed', 'Orphan_End', 'Self_Circular',
            'Extension_Rate(%)', 'Orphan_Rate(%)',
            'Extended_Contigs', 'New_Contigs', 'Total_Gain(bp)',
            'Has_Cobra', 'Has_Virus', 'Cobra_FA', 'Virus_FA']
    with open(out / "cobra_summary.tsv", 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in sample_stats:
            f.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')

    # 写入 contig detail
    if contig_details:
        dcols = ['Sample', 'Contig_ID', 'Virus_Len(bp)', 'Cobra_Len(bp)', 'Gain(bp)', 'Status']
        with open(out / "cobra_contig_detail.tsv", 'w') as f:
            f.write('\t'.join(dcols) + '\n')
            for r in contig_details:
                f.write('\t'.join(str(r.get(c, '')) for c in dcols) + '\n')

    # 打印汇总
    n_with_cobra = sum(1 for r in sample_stats if r['Has_Cobra'])
    n_extended_total = sum(r['Extended_Contigs'] for r in sample_stats)
    print(f"\n  汇总: {len(sample_stats)} 样本, {n_with_cobra} 有延伸结果, {n_extended_total} 条延伸")
    print(f"  cobra_summary.tsv → {out / 'cobra_summary.tsv'}")
    if contig_details:
        print(f"  cobra_contig_detail.tsv → {out / 'cobra_contig_detail.tsv'}")

def main():
    p = argparse.ArgumentParser(description="简版 COBRA 延伸统计 (无 CheckV)")
    p.add_argument('-c', '--cobra_dir', required=True, help='03_COBRA 目录')
    p.add_argument('-o', '--output', default='cobra_stats', help='输出目录 (默认: cobra_stats/)')
    args = p.parse_args()
    analyze_cobra_dir(args.cobra_dir, args.output)

if __name__ == '__main__':
    main()
