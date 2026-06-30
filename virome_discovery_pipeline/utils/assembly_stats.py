#!/usr/bin/env python3
"""
组装统计: N50/N90 + contig 长度分布 → assembly_summary.tsv
用法: python assembly_stats.py -a 01_Assembly/ [-o assembly_summary.tsv]
"""
import argparse, os, sys
from pathlib import Path

def stats_fasta(fasta_path):
    """返回 (n_contigs, total_bp, max_len, n50, n90, n_500bp, n_1000bp, ratio_500, ratio_1000)"""
    lens = []; seq = ""
    for line in open(fasta_path):
        l = line.strip()
        if l.startswith('>'):
            if seq: lens.append(len(seq))
            seq = ""
        else: seq += l
    if seq: lens.append(len(seq))
    if not lens: return (0,0,0,0,0,0,0,0.0,0.0)

    lens.sort(reverse=True)
    total = sum(lens); cum = 0
    half = total / 2; n90t = total * 0.9
    n50 = n90 = 0
    for la in lens:
        cum += la
        if n50 == 0 and cum >= half: n50 = la
        if n90 == 0 and cum >= n90t: n90 = la
        if n50 > 0 and n90 > 0: break
    c500 = sum(1 for la in lens if la > 500)
    c1000 = sum(1 for la in lens if la > 1000)
    r500 = round(c500 / len(lens) * 100, 1)
    r1000 = round(c1000 / len(lens) * 100, 1)
    return (len(lens), total, lens[0], n50, n90, c500, c1000, r500, r1000)

def main():
    p = argparse.ArgumentParser(description="组装统计: N50/N90 + 长度分布")
    p.add_argument('-a', '--assembly_dir', required=True, help='01_Assembly 目录')
    p.add_argument('-o', '--output', help='输出 TSV (默认: {assembly_dir}/assembly_summary.tsv)')
    p.add_argument('--size_unit', default='Mb', choices=['Mb','Kb','bp'], help='大小单位')
    args = p.parse_args()

    asm = Path(args.assembly_dir)
    if not asm.is_dir():
        sys.exit(f"目录不存在: {asm}")

    out = Path(args.output) if args.output else asm / "assembly_summary.tsv"

    div = {'Mb': 1e6, 'Kb': 1e3, 'bp': 1}[args.size_unit]
    rows = []
    for d in sorted(asm.iterdir()):
        if not d.is_dir(): continue
        for f in sorted(d.glob("*.contig.fasta")):
            sn = d.name
            at = f.stem.replace(f"{sn}_", "").replace(".contig", "")
            try:
                n, total_bp, mx, n50, n90, c500, c1000, r500, r1000 = stats_fasta(str(f))
            except: continue
            if n == 0: continue
            rows.append({
                'Sample': sn, 'Assembler': at,
                f'Size({args.size_unit})': round(total_bp/div, 2),
                'Contigs': n, 'Max_Len': mx, 'N50': n50, 'N90': n90,
                '>500bp': c500, '>500bp(%)': r500,
                '>1000bp': c1000, '>1000bp(%)': r1000,
            })

    if not rows:
        print("未找到组装结果"); return

    cols = ['Sample','Assembler',f'Size({args.size_unit})','Contigs','Max_Len','N50','N90',
            '>500bp','>500bp(%)','>1000bp','>1000bp(%)']
    with open(out, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')

    # 简要打印
    print(f"  {'Sample':<20} {'Size':>8} {'Contigs':>8} {'N50':>8} {'N90':>8} {'>500bp':>8}")
    print("  " + "-"*65)
    for r in rows:
        print(f"  {r['Sample']:<20} {r[f'Size({args.size_unit})']:>8} {r['Contigs']:>8} {r['N50']:>8} {r['N90']:>8} {r['>500bp']:>8}")
    print(f"\n  输出: {out}")

if __name__ == '__main__':
    main()
