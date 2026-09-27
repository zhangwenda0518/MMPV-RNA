#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""看 09 新增的 DNA/RNA 三列在这 625 行冲突上取了哪一侧（种侧还是科侧）。"""
import csv, os
from collections import Counter

P = ("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/"
     "09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv")
VMR = "/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv"

sp2fam, fam2gnm = {}, {}
with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t"); next(rd)
    for r in rd:
        if len(r) < 27: continue
        fam, sp, gnm = r[13].strip().strip('"'), r[17].strip().strip('"'), r[25].strip().strip('"')
        if sp: sp2fam.setdefault(sp, fam)
        if fam: fam2gnm.setdefault(fam, gnm)

c = Counter()
with open(P, newline="", encoding="utf-8", errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        key = (r["Family"], r.get("Genome_Type", ""), r.get("Genome_Structure", ""))
        c[key] += 1

print("=== barbarum All_plant：Family 与新增 DNA/RNA 列的搭配（top 20）===")
print("%-24s %-10s %-12s %-8s %s" % ("Family", "Genome_Type", "Genome_Structure", "行数", "该科的VMR基因组"))
for k, v in sorted(c.items(), key=lambda kv: -kv[1])[:20]:
    print("%-24s %-10s %-12s %-8d %s" % (k[0], k[1], k[2], v, fam2gnm.get(k[0], "?")))

print("\n=== 只看 Family=Mimiviridae 的行 ===")
n = 0
with open(P, newline="", encoding="utf-8", errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["Family"].strip() != "Mimiviridae": continue
        n += 1
        if n <= 8:
            print("  %-46s Genus=%-12s Species=%-24s %s/%s" % (
                r["contig_id"][:46], r["Genus"], r["Species"],
                r.get("Genome_Type", ""), r.get("Genome_Structure", "")))
print("  Mimiviridae 行合计 %d" % n)
