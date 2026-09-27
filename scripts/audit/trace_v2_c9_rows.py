#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""只看 V2-V1 / V1-V0 差集行的 C9 原始判定 (直接按 contig_id 连接, 不做全表正则扫描)。"""
import os

import pandas as pd

C9 = "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus/06_HostPrediction/C9_ICTV_result/classification_result.tsv"
KEEP = ["contig_id", "primary_tool", "confidence", "Family", "Genus", "Species",
        "Family_agree", "Genus_agree", "Species_agree", "Predicted_Host",
        "Confidence_Level", "Integrated_Confidence", "Determination_Level", "_pass_conf"]

c9 = pd.read_csv(C9, sep="\t", dtype=str)
c9 = c9[[c for c in KEEP if c in c9.columns]]

for f, want in [("/tmp/modes_v2minusv1_onekp_83.tsv", "V2-V1(83)"),
                ("/tmp/modes_v1minusv0_onekp_407.tsv", "V1-V0(407)")]:
    if not os.path.isfile(f):
        print("%s 不存在 (先跑 compare_veto_modes2.py)" % f)
        continue
    d = pd.read_csv(f, sep="\t", dtype=str)
    j = d.merge(c9, on="contig_id", how="left", suffixes=("_05", "_c9"))
    print("\n=========== %s : %d 行 ===========" % (want, len(j)))
    for base in ("Predicted_Host", "Determination_Level"):
        col = base if base in j.columns else base + "_c9"
        if col in j.columns:
            print("C9 %s 分布: %s" % (base, j[col].value_counts().to_dict()))
    fc = "Family_c9" if "Family_c9" in j.columns else "Family"
    gc = "Genus_c9" if "Genus_c9" in j.columns else "Genus"
    print("C9 自身 (Family, Genus) 分布 top12:")
    print(j.groupby([fc, gc]).size().sort_values(ascending=False).head(12).to_string())
    cols = ["contig_id", fc, gc, "Species", "Predicted_Host",
            "Determination_Level", "primary_tool", "Genus_agree", "_pass_conf"]
    print(j[[c for c in cols if c in j.columns]].head(10).to_string(max_colwidth=24, index=False))
