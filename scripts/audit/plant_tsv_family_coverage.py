#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""12 个 ICTV 含植物科在 Plant.tsv 里的覆盖 (科级行数 / 属数 / top 属), 位置口径。"""
from collections import Counter, defaultdict

PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"
FAMS = ["Artoviridae", "Chrysoviridae", "Genomoviridae", "Kanorauviridae", "Mitoviridae",
        "Ourmiaviridae", "Pestiviridae", "Pseudoviridae", "Spiciviridae", "Tomosaviridae",
        "Virgaviridae", "Solemoviridae"]

cnt = defaultdict(int)
gen = defaultdict(Counter)
tot = 0
with open(PLANT, encoding="utf-8", errors="replace") as fh:
    fh.readline()
    for line in fh:
        fs = line.rstrip("\n").split("\t")
        if len(fs) < 8:
            continue
        tot += 1
        p = fs[3].split(";")
        if len(p) > 5:
            f = p[5].strip()
            cnt[f] += 1
            if len(p) > 6 and p[6].strip():
                gen[f][p[6].strip()] += 1

print("Plant.tsv 行 %d" % tot)
print("%-18s %7s %5s  %s" % ("Family", "rows", "genus", "top genera"))
for f in FAMS:
    top = ", ".join("%s %d" % (g, n) for g, n in gen[f].most_common(8)) or "(零覆盖)"
    print("%-18s %7d %5d  %s" % (f, cnt[f], len(gen[f]), top))
