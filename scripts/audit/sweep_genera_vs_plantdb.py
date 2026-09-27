#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""反查属黑名单 (66) 与植物库的冲突: 黑名单属在 Plant.tsv 中是否有记录;
并给出 Plant.tsv 的属清单与记录数, 用于判断别名/改名导致的误登。"""
import importlib.util
import os
from collections import Counter

import pandas as pd

P = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

PLANT = os.path.expanduser("~/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv")
pl = pd.read_csv(PLANT, sep="\t", low_memory=False, usecols=["Virus_lineage"])
parts = [str(x).split(";") for x in pl["Virus_lineage"]]
fams, gens = [], []
for p in parts:
    f = g = ""
    for i, tok in enumerate(p):
        if tok.endswith("viridae"):
            f = tok
            g = p[i + 1] if i + 1 < len(p) else ""
            break
    fams.append(f)
    gens.append(g)
print("Plant.tsv %d 行 → 科 %d 个 / 属 %d 个" % (len(pl), len({f for f in fams if f}), len({g for g in gens if g})))
cg = Counter(g for g in gens if g)
cf = Counter(f for f in fams if f)
print("\n植物库科 Top: %s" % cf.most_common(8))
print("植物库属 Top10: %s" % cg.most_common(10))

bl = set(m.NON_PLANT_GENERA)
inter = {g: cg[g] for g in bl if g in cg}
print("\n【冲突】属黑名单中的属在植物库有记录 (%d 个): %s" % (len(inter), inter))
print("\n属黑名单中在植物库零记录的: %d / %d" % (len(bl) - len(inter), len(bl)))
print("\n植物库属清单 (%d):" % len(cg))
items = sorted(cg.items(), key=lambda x: -x[1])
for i in range(0, len(items), 5):
    print("  " + " ".join("%-24s" % ("%s(%d)" % (g, c)) for g, c in items[i:i + 5]))
