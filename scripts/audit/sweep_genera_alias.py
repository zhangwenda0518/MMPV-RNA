#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""属黑名单的别名风险检查: 植物库是否存在以该属名开头的物种名 (即植物库其实认得这个属,
只是没有属级记录)。同时给出属级记录数。"""
import importlib.util
import os
import re
from collections import Counter

import pandas as pd

P = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

PLANT = os.path.expanduser("~/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv")
pl = pd.read_csv(PLANT, sep="\t", low_memory=False, usecols=["Virus_lineage"])
gens, specs = Counter(), Counter()
for x in pl["Virus_lineage"]:
    parts = [t.strip() for t in str(x).split(";")]
    fi = next((i for i, t in enumerate(parts) if t.endswith("viridae")), None)
    if fi is None:
        continue
    if fi + 1 < len(parts) and parts[fi + 1]:
        gens[parts[fi + 1]] += 1
    if fi + 2 < len(parts) and parts[fi + 2]:
        specs[parts[fi + 2]] += 1
sp_list = list(specs)
print("植物库: 属级记录 %d 个属 / 物种名 %d 个" % (len(gens), len(sp_list)))

print("\n黑名单属逐条核查 (属级记录数 / 植物库是否存在以该属名开头的物种名):")
conflict = []
for g in sorted(m.NON_PLANT_GENERA):
    n_gen = gens.get(g, 0)
    pat = re.compile(r"^" + re.escape(g) + r"[\s_]", re.I)
    ex = [s for s in sp_list if pat.match(s)]
    flag = "冲突" if (n_gen or ex) else "ok"
    if n_gen or ex:
        conflict.append(g)
    print("   %-24s 属记录=%-4d 物种名匹配=%-3d %s %s"
          % (g, n_gen, len(ex), flag, ex[:2]))
print("\n【建议撤回】的属 (%d): %s" % (len(conflict), conflict))
