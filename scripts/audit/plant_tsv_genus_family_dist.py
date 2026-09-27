#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""查指定属在 Plant.tsv 里的科归属分布 (位置口径), 用于核对嵌合行里植物属的来源。"""
import sys
from collections import Counter

PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"
TARGETS = sys.argv[1:] or ["Cilevirus", "Potyvirus", "Orthotospovirus", "Maculavirus",
                           "Totivirus", "Begomovirus", "Tungrovirus", "Caulimovirus", "Ourmiavirus"]
want = set(TARGETS)
dist = {t: Counter() for t in TARGETS}
with open(PLANT, encoding="utf-8", errors="replace") as fh:
    fh.readline()
    for line in fh:
        fs = line.rstrip("\n").split("\t")
        if len(fs) < 8:
            continue
        p = fs[3].split(";")
        if len(p) > 6 and p[6].strip() in want:
            dist[p[6].strip()][p[5].strip() if len(p) > 5 else "(空)"] += 1
for t in TARGETS:
    print("%-16s %s" % (t, ", ".join("%s %d" % (k, v) for k, v in dist[t].most_common(4)) or "(无记录)"))
