#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""校验服务器上 run_host_prediction.py 的名单现状 + 抽查 species 宿主表可疑条目。"""
import importlib.util
import os

import pandas as pd

P = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", P)
rhp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rhp)

print("科兜底名单 %d 条 ; 属名单 %d 条 ; 待复核清单 %s"
      % (len(rhp.NON_PLANT_FAMILIES_FALLBACK), len(rhp.NON_PLANT_GENERA),
         getattr(rhp, "NON_PLANT_FAMILIES_UNDER_REVIEW", "N/A")))
for f in ["Ourmiaviridae", "Mitoviridae", "Kanorauviridae", "Metaviridae",
          "Pithoviridae", "Mimiviridae", "Marseilleviridae", "Barnaviridae", "Narnaviridae"]:
    print("   %-18s in-list=%s" % (f, f in rhp.NON_PLANT_FAMILIES_FALLBACK))
print("\nis_blacklisted 抽样:")
for fam, gen, det in [("Mimiviridae", "", "Family(via Family)"),
                      ("Mimiviridae", "Potyvirus", "Species(via Species)"),
                      ("Ourmiaviridae", "", "Family(via Family)"),
                      ("Partitiviridae", "", "Family(via Family)"),
                      ("Potyviridae", "", "Family(via Family)")]:
    print("   %-16s %-10s %-20s -> %s" % (fam, gen or "-", det,
                                          rhp.is_blacklisted(fam, gen, det)))

PROB = "/home/zhangwenda/MMPV-RNA/database/cross_analysis/species_host_probability.tsv"
if os.path.isfile(PROB):
    sp = pd.read_csv(PROB, sep="\t")
    print("\nspecies 宿主表列: %s (%d 行)" % (list(sp.columns), len(sp)))
    key = sp.columns[0]
    for name in ["Bactericera cockerelli picorna-like virus",
                 "Alphambiguivirus magnaporthensis",
                 "Grapevine wood holobiome associated mycobunyavirales-like virus 2",
                 "Potyvirus rapae"]:
        hit = sp[sp[key].astype(str).str.strip() == name]
        print("   %-60s -> %s" % (name, hit.to_dict("records")[:1] if len(hit) else "无记录"))
