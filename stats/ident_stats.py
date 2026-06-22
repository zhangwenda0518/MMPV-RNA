#!/usr/bin/env python3
"""
02_Identification 后处理: 工具一致性 + 过滤统计 + Venn 数据
用法: python ident_stats.py -i 02_Identification/ [-o ident_stats/]
产出:
  ident_summary.tsv      每样本每工具鉴定数量
  filter_summary.tsv     每样本 filter/strict 过滤前后对比
  tool_overlap.tsv       工具间重叠矩阵
"""
import argparse, os, sys
from pathlib import Path
from collections import defaultdict, Counter

def count_fasta(fpath):
    if not fpath.is_file(): return set()
    ids = set()
    for line in open(fpath):
        if line.startswith('>'): ids.add(line[1:].split()[0])
    return ids

def main():
    p = argparse.ArgumentParser(description="02_Identification 后处理统计")
    p.add_argument('-i', '--ident_dir', required=True, help='02_Identification 目录')
    p.add_argument('-o', '--output', default='ident_stats', help='输出目录')
    args = p.parse_args()

    ident = Path(args.ident_dir)
    if not ident.is_dir(): sys.exit(f"目录不存在: {ident}")

    out = Path(args.output); out.mkdir(parents=True, exist_ok=True)

    # Tool patterns
    tools_list = ['genomad','blast','metabuli','virsorter2','viralverify',
                  'virhunter','virbot','viralm','rdrpcatch']
    samples = [d for d in sorted(ident.iterdir()) if d.is_dir()]

    # ── 1. ident_summary.tsv ──
    ident_rows = []
    for sd in samples:
        all_fa = sd / f"{sd.name}_virus.all.candidate.fasta"
        all_ids = count_fasta(all_fa)
        row = {'Sample': sd.name, 'All_Candidate': len(all_ids)}
        for tool in tools_list:
            id_file = sd / f"{sd.name}_virus.{tool}.result.id"
            if id_file.is_file():
                tids = set(open(id_file).read().strip().split('\n'))
                row[f'{tool}'] = len(tids)
            else:
                row[f'{tool}'] = 0
        ident_rows.append(row)

    cols = ['Sample','All_Candidate'] + tools_list
    with open(out / "ident_summary.tsv", 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in ident_rows:
            f.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')

    # ── 2. filter_summary.tsv ──
    filter_rows = []
    for sd in samples:
        all_ids = count_fasta(sd / f"{sd.name}_virus.all.candidate.fasta")
        n_all = len(all_ids)
        for filt_name, filt_dir in [('filter','uniprot_filter_output_filter'),
                                     ('strict','uniprot_filter_output_strict'),
                                     ('comb','uniprot_filter_output')]:
            filt_fa = sd / filt_dir / f"{sd.name}_virus.uniprot_filtered.fasta"
            n_filt = len(count_fasta(filt_fa)) if filt_fa.is_file() else 0
            if n_filt > 0 or filt_dir == 'uniprot_filter_output':
                filter_rows.append({'Sample': sd.name, 'Mode': filt_name,
                    'All_Candidate': n_all, 'Passed': n_filt,
                    'Retained(%)': round(n_filt/max(n_all,1)*100, 1)})

    with open(out / "filter_summary.tsv", 'w') as f:
        cols2 = ['Sample','Mode','All_Candidate','Passed','Retained(%)']
        f.write('\t'.join(cols2) + '\n')
        for r in filter_rows:
            f.write('\t'.join(str(r[c]) for c in cols2) + '\n')

    # ── 3. tool_overlap.tsv (每个样本的工具重叠对) ──
    with open(out / "tool_overlap.tsv", 'w') as f:
        f.write("Sample\tTool_A\tTool_B\tA_Only\tB_Only\tOverlap\tUnion\tJaccard(%)\n")
        for sd in samples:
            tool_ids = {}
            for tool in tools_list:
                id_file = sd / f"{sd.name}_virus.{tool}.result.id"
                if id_file.is_file():
                    tool_ids[tool] = set(open(id_file).read().strip().split('\n'))
            active = list(tool_ids.keys())
            for i in range(len(active)):
                for j in range(i+1, len(active)):
                    a, b = active[i], active[j]
                    sa, sb = tool_ids[a], tool_ids[b]
                    intersect = len(sa & sb)
                    union = len(sa | sb)
                    jaccard = round(intersect/max(union,1)*100, 1) if union > 0 else 0
                    f.write(f"{sd.name}\t{a}\t{b}\t{len(sa-sb)}\t{len(sb-sa)}\t{intersect}\t{union}\t{jaccard}\n")

    # ── 打印 ──
    print(f"  样本: {len(samples)}")
    for r in ident_rows[:5]:
        top_tools = sorted([(t, r[t]) for t in tools_list if r.get(t,0) > 0], key=lambda x:-x[1])[:5]
        print(f"  {r['Sample']}: all={r['All_Candidate']}, " +
              ", ".join(f"{t}={c}" for t,c in top_tools))
    if len(ident_rows) > 5: print(f"  ... +{len(ident_rows)-5} 样本")
    print(f"\n  {out}/ident_summary.tsv")
    print(f"  {out}/filter_summary.tsv")
    print(f"  {out}/tool_overlap.tsv")

if __name__ == '__main__':
    main()
