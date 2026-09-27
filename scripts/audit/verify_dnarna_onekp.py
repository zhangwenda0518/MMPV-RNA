#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""onekp 两表 DNA/RNA 计数（补 verify 脚本 glob 漏掉的部分）"""
import csv, os, glob
from collections import Counter
ROOT = os.path.expanduser("~/MMPV-paper")

def rows(p):
    with open(p, newline="", encoding="utf-8", errors="surrogateescape") as f:
        return list(csv.reader(f, delimiter="\t"))

for rel, tag in [("09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv", "All_plant"),
                 ("09_Virome_Analysis/HQ_analysis/HQ_plant_viruses_info.tsv", "HQ"),
                 ("09_Virome_Analysis/all_plant_analysis/all_plant_viruses_genus_summary.tsv", "genus_summary")]:
    for q in sorted(glob.glob(os.path.join(ROOT, "*", "onekp-virus", rel))):
        r = rows(q)
        idx = {h: i for i, h in enumerate(r[0])}
        c = Counter((x[idx["Genome_Type"]] if idx["Genome_Type"] < len(x) else "") or "未定" for x in r[1:] if x and x[0])
        s = Counter(x[idx["Genome_Source"]] for x in r[1:] if x and x[0])
        print("%-14s onekp  DNA=%-5d RNA=%-6d 未定=%-5d | species=%d genus=%d tax=%d ambiguous=%d none=%d" % (
            tag, c["DNA"], c["RNA"], c["未定"], s["species"], s["genus"], s["tax"], s["genus_ambiguous"], s["none"]))

# genus_summary 全部项目（含 NA 属）
print()
tot = Counter()
for q in sorted(glob.glob(os.path.join(ROOT, "*", "*", "*", "09_Virome_Analysis/all_plant_analysis/all_plant_viruses_genus_summary.tsv"))) + \
         sorted(glob.glob(os.path.join(ROOT, "*", "onekp-virus", "09_Virome_Analysis/all_plant_analysis/all_plant_viruses_genus_summary.tsv"))):
    r = rows(q)
    idx = {h: i for i, h in enumerate(r[0])}
    c = Counter(x[idx["Genome_Type"]] or "未定" for x in r[1:] if x and x[0])
    tot.update(c)
    print("  %-58s DNA=%-3d RNA=%-4d 未定=%-3d" % (q.split("MMPV-paper/")[-1][-58:], c["DNA"], c["RNA"], c["未定"]))
print("  genus_summary 合计 DNA=%d RNA=%d 未定=%d" % (tot["DNA"], tot["RNA"], tot["未定"]))
