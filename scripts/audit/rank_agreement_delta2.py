#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""修正版: 每棵树用各自的 C9 表; 并对放行行的 C9 Species 名去 Plant.tsv(Host 库) 里查证。

用法: python3 /tmp/rank_agreement_delta2.py
"""
import os
import re
from collections import Counter

import pandas as pd

TREES = {
    "goji": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out",
    "onekp": "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus",
}
PLANT_TSV = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"


def frac(s):
    m = re.match(r"^\s*(\d+)\s*/\s*(\d+)", str(s))
    return (int(m.group(1)), int(m.group(2))) if m else None


def load_plant_names():
    names = set()
    for chunk in pd.read_csv(PLANT_TSV, sep="\t", dtype=str, usecols=lambda c: c in
                             ("Virus_Name", "Host_Name", "Virus_lineage"), chunksize=50000):
        if "Virus_Name" in chunk.columns:
            names |= set(chunk["Virus_Name"].dropna().str.strip().str.lower())
        if "Host_Name" in chunk.columns:
            names |= set(chunk["Host_Name"].dropna().str.strip().str.lower())
    return names


PLANT_NAMES = load_plant_names()
print("Plant.tsv 名称集合 %d 个" % len(PLANT_NAMES))

for tag, tree in TREES.items():
    c9 = pd.read_csv(os.path.join(tree, "06_HostPrediction/C9_ICTV_result/classification_result.tsv"),
                     sep="\t", dtype=str)
    print("\n########## %s : C9 %d 行, Predicted_Host=Plant %d" %
          (tag, len(c9), (c9["Predicted_Host"] == "Plant").sum()))
    for f, sub in [("/tmp/modes_v1minusv0_%s" % tag, "V1-V0"), ("/tmp/modes_v2minusv1_%s" % tag, "V2-V1")]:
        import glob
        fs = sorted(glob.glob(f + "_*.tsv"))
        if not fs:
            print("  %s 明细缺失 (%s)" % (sub, f + "_*.tsv"))
            continue
        ids = pd.read_csv(fs[0], sep="\t", dtype=str)["contig_id"]
        j = c9[c9["contig_id"].isin(set(ids))].copy()
        print("\n  ==== %s (%s): %d 行, 匹配 C9 %d 行" % (sub, os.path.basename(fs[0]), len(ids), len(j)))
        if not len(j):
            continue
        for rank in ("Order", "Family", "Genus", "Species"):
            col = rank + "_agree"
            if col not in j.columns:
                continue
            vals = [frac(v) for v in j[col]]
            unan = sum(1 for v in vals if v and v[0] == v[1] and v[1] > 0)
            mino = sum(1 for v in vals if v and v[0] < v[1])
            zero = sum(1 for v in vals if v and v[1] == 0)
            print("    %-8s 一致 %3d / 少数 %3d / 0票 %3d" % (rank, unan, mino, zero))
        sp = j["Species"].dropna().str.strip().str.lower()
        hit = sum(1 for s in sp if s in PLANT_NAMES)
        print("    Species 名在 Plant.tsv 命中 %d/%d" % (hit, len(sp)))
        print("    Species top8: %s" % Counter(j["Species"]).most_common(8))
        print("    Family top8 : %s" % Counter(j["Family"]).most_common(8))
        print("    Genus top8  : %s" % Counter(j["Genus"]).most_common(8))
        print("    DetLevel    : %s" % Counter(j["Determination_Level"]).most_common(6))
        print("    ConfLevel   : %s" % Counter(j["Confidence_Level"]).most_common(6))
