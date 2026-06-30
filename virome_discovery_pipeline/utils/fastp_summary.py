#!/usr/bin/env python3
"""
汇总 fastp JSON 报告 → data_summary.tsv
用法: python fastp_summary.py --clean_dir 00a_CleanData/
输出: 00a_CleanData/data_summary.tsv
"""
import argparse, json, os, sys
from pathlib import Path

def main():
    p = argparse.ArgumentParser(description="汇总 fastp JSON → data_summary.tsv")
    p.add_argument('--clean_dir', '-c', required=True, help='00a_CleanData 目录')
    p.add_argument('--output', '-o', help='输出路径 (默认: {clean_dir}/data_summary.tsv)')
    args = p.parse_args()

    clean = Path(args.clean_dir)
    log_dir = clean / 'logs'
    if not log_dir.is_dir():
        sys.exit(f"logs 目录不存在: {log_dir}")

    json_files = sorted(log_dir.glob('*_fastp_report.json'))
    if not json_files:
        sys.exit(f"未找到 *_fastp_report.json")

    rows = []
    for jf in json_files:
        sample = jf.name.replace('_fastp_report.json', '')
        try:
            d = json.load(open(jf))
        except Exception as e:
            print(f"  [WARN] 跳过 {jf.name}: {e}")
            continue

        s = d.get('summary', {})
        bef = s.get('before_filtering', {})
        aft = s.get('after_filtering', {})
        fil = d.get('filtering_result', {})
        dup = d.get('duplication', {})

        rows.append({
            'Sample': sample,
            'Raw_Reads': bef.get('total_reads', 0),
            'Raw_Bases': bef.get('total_bases', 0),
            'Raw_Q20(%)': round(bef.get('q20_rate', 0) * 100, 2),
            'Raw_Q30(%)': round(bef.get('q30_rate', 0) * 100, 2),
            'Raw_GC(%)': round(bef.get('gc_content', 0) * 100, 2),
            'Clean_Reads': aft.get('total_reads', 0),
            'Clean_Bases': aft.get('total_bases', 0),
            'Clean_Q20(%)': round(aft.get('q20_rate', 0) * 100, 2),
            'Clean_Q30(%)': round(aft.get('q30_rate', 0) * 100, 2),
            'Clean_GC(%)': round(aft.get('gc_content', 0) * 100, 2),
            'Reads_Retained(%)': round(
                aft.get('total_reads', 0) / max(bef.get('total_reads', 1), 1) * 100, 2),
            'Bases_Retained(%)': round(
                aft.get('total_bases', 0) / max(bef.get('total_bases', 1), 1) * 100, 2),
            'LowQ_Reads': fil.get('low_quality_reads', 0),
            'TooManyN_Reads': fil.get('too_many_N_reads', 0),
            'TooShort_Reads': fil.get('too_short_reads', 0),
            'Corrected_Reads': fil.get('corrected_reads', 0),
            'Duplication_Rate(%)': round(dup.get('rate', 0) * 100, 4),
        })

    out = Path(args.output) if args.output else clean / 'data_summary.tsv'
    cols = list(rows[0].keys()) if rows else []
    with open(out, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')

    # 汇总行
    if len(rows) > 1:
        totals = {c: sum(r[c] for r in rows) if isinstance(rows[0][c], (int, float)) and 'Rate' not in c and '%' not in c else '' for c in cols}
        totals['Sample'] = 'TOTAL'
        with open(out, 'a') as f:
            f.write('\t'.join(str(totals.get(c, '')) for c in cols) + '\n')

    print(f"  生成: {out} ({len(rows)} 样本)")
    for r in rows:
        print(f"    {r['Sample']}: {r['Raw_Reads']:,} → {r['Clean_Reads']:,} reads "
              f"({r['Reads_Retained(%)']}%), dup={r['Duplication_Rate(%)']}%")

if __name__ == '__main__':
    main()
