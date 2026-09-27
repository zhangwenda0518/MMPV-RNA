#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""追 V2 额外放行行的上游来源: C9 为什么把 Arenavirus/Simplexvirus 这类属判成 Plant?
查 (1) C9 结果表里这些 contig 的全列; (2) 属级概率表里这些属的 Plant_Records/Predicted_Host。"""
import glob
import os

import pandas as pd

OUT = "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus/06_HostPrediction"
C9 = os.path.join(OUT, "C9_ICTV_result/classification_result.tsv")
DB = "/home/zhangwenda/MMPV-RNA/database/cross_analysis"
TARGETS = ["Arenavirus", "Simplexvirus", "Lymphocystivirus", "Alphabaculovirus",
           "Chlorovirus", "Prasinovirus", "Klosneuvirus", "Iltovirus",
           "Potyvirus", "Macluravirus", "Badnavirus"]

c9 = pd.read_csv(C9, sep="\t", dtype=str)
print("C9 表 %d 行 列: %s" % (len(c9), list(c9.columns)))
print("Predicted_Host 分布: %s" % c9["Predicted_Host"].value_counts().to_dict())

for t in TARGETS:
    hit = c9[c9.apply(lambda r: r.astype(str).str.contains(t, case=False, regex=False).any(), axis=1)]
    if len(hit):
        print("\n--- C9 命中 %s: %d 行" % (t, len(hit)))
        print(hit.head(2).to_string(max_colwidth=28))

print("\n\n=== 属级概率表 ===")
for f in sorted(glob.glob(os.path.join(DB, "*.tsv")) + glob.glob(os.path.join(DB, "*.csv"))):
    try:
        d = pd.read_csv(f, sep="\t" if f.endswith(".tsv") else ",", dtype=str, low_memory=False)
    except Exception as e:
        print("%s 读取失败 %s" % (os.path.basename(f), e))
        continue
    cols = list(d.columns)
    gcol = next((c for c in cols if c.lower() in ("genus", "taxon", "name")), cols[0])
    if not any("plant" in c.lower() for c in cols):
        continue
    print("\n%s (%d 行) 列: %s" % (os.path.basename(f), len(d), cols[:12]))
    sub = d[d[gcol].astype(str).isin(TARGETS)]
    if len(sub):
        keep = [c for c in [gcol, "Plant_Records", "Total_Records", "Predicted_Host", "Host",
                            "Family", "top_host", "Host_Category"] if c in sub.columns]
        print(sub[keep].to_string(index=False, max_colwidth=24))
