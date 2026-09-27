#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""属黑名单的系统性冲突检查: 每个黑名单属的参照科 (NCBI/VMR) 是否属于植物库的 46 科。
科属植物库已知 → 该属不应被登记为非植物属 (改名/别名导致的误登)。"""
import importlib.util
import os
from collections import Counter

import pandas as pd

P = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

PLANT = os.path.expanduser("~/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv")
REF = os.path.expanduser("~/database/taxonomy/genus_family_ref.tsv")
pl = pd.read_csv(PLANT, sep="\t", low_memory=False, usecols=["Virus_lineage"])
pf = set()
for x in pl["Virus_lineage"]:
    for tok in str(x).split(";"):
        if tok.endswith("viridae"):
            pf.add(tok)
            break
ref = pd.read_csv(REF, sep="\t")
g2n = dict(zip(ref["Genus"].astype(str), ref["NCBI_Family"].astype(str)))
g2v = dict(zip(ref["Genus"].astype(str), ref["VMR_Family"].astype(str)))
print("植物库科 %d 个; 属黑名单 %d 个" % (len(pf), len(m.NON_PLANT_GENERA)))
hits, missing = [], []
for g in sorted(m.NON_PLANT_GENERA):
    nf, vf = g2n.get(g, ""), g2v.get(g, "")
    if not nf and not vf:
        missing.append(g)
    if nf in pf or vf in pf:
        hits.append((g, nf, vf))
print("\n【冲突】黑名单属的参照科属于植物库已知科 (%d):" % len(hits))
for g, nf, vf in hits:
    print("   %-22s NCBI=%-24s VMR=%s" % (g, nf or "-", vf or "-"))
print("\n【参照表无记录】的属 (%d): %s" % (len(missing), missing))
print("\n全部黑名单属的参照科分布:")
for g in sorted(m.NON_PLANT_GENERA):
    print("   %-24s NCBI=%-26s VMR=%s" % (g, g2n.get(g, "-"), g2v.get(g, "-")))
