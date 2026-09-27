#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""关键判据: V1/V2 放行行里, 是"科"票弱还是"种"票弱?
逐行取 C9 的 Family_agree / Genus_agree / Species_agree, 解析 "n/m: tools" 形式,
看族级票是否少数派, 而种级票是否一致。"""
import os
import re
from collections import Counter

import pandas as pd

C9 = "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus/06_HostPrediction/C9_ICTV_result/classification_result.tsv"
c9 = pd.read_csv(C9, sep="\t", dtype=str)


def frac(s):
    m = re.match(r"^\s*(\d+)\s*/\s*(\d+)", str(s))
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


SETS = [("/tmp/modes_v1minusv0_onekp_407.tsv", "V1-V0 onekp (407)"),
        ("/tmp/modes_v2minusv1_onekp_83.tsv", "V2-V1 onekp (83)"),
        ("/tmp/modes_v1minusv0_goji_13.tsv", "V1-V0 goji (13)"),
        ("/tmp/modes_v2minusv1_goji_19.tsv", "V2-V1 goji (19)")]

for f, tag in SETS:
    if not os.path.isfile(f):
        print("%s 缺失" % f)
        continue
    ids = pd.read_csv(f, sep="\t", dtype=str)["contig_id"]
    j = c9[c9["contig_id"].isin(set(ids))].copy()
    print("\n==== %s: 匹配到 C9 %d 行 ====" % (tag, len(j)))
    for rank in ("Family", "Genus", "Species"):
        col = rank + "_agree"
        if col not in j.columns:
            continue
        vals = [frac(v) for v in j[col]]
        unanimous = sum(1 for v in vals if v and v[0] == v[1] and v[1] > 0)
        minority = sum(1 for v in vals if v and v[0] < v[1])
        zero = sum(1 for v in vals if v and v[1] == 0)
        none = sum(1 for v in vals if v is None)
        print("  %-8s 票一致 %3d / 少数派 %3d / 0投票 %3d / 无值 %3d   示例: %s"
              % (rank, unanimous, minority, zero, none, list(j[col].head(3))))
    print("  primary_tool 分布: %s" % Counter(j["primary_tool"]).most_common(5))
    print("  Confidence_Level 分布: %s" % Counter(j["Confidence_Level"]).most_common(5))
    print("  这些行的科标签 top8: %s" % Counter(j["Family"]).most_common(8))
    # 种名里出现属/种关键词的数量
    print("  Species 示例: %s" % list(j["Species"].head(6)))
