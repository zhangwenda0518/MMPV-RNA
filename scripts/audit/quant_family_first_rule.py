#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""量化 "科级否决优先" 规则对 goji Plant 行的影响 (不改动服务器源码, 仅内存替换函数)。

现规则: 属/种级判定 (trusted) 直接放行 → 科级兜底被绕过;
       属非空且非黑名单属 → 也不查科。
候选规则 V2: 科在非植物科名单 → 直接否决 Plant (与判定层级无关), 其余逻辑不变。
"""
import importlib.util
import os
from collections import Counter

import pandas as pd

RHP = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", RHP)
rhp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rhp)

TREE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
TAX = os.path.join(TREE, "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")
OUT = os.path.join(TREE, "06_HostPrediction")
C9 = os.path.join(OUT, "C9_ICTV_result/classification_result.tsv")
RVH = os.path.join(OUT, "RVH_result/result.csv")
PB2 = os.path.join(OUT, "phabox2_output/final_prediction/cherry_prediction.tsv")

df = pd.read_csv(TAX, sep="\t")
c9 = pd.read_csv(C9, sep="\t")[["contig_id", "Predicted_Host", "Determination_Level"]].rename(
    columns={"Predicted_Host": "Host_ICTV"})
df = df.merge(c9, on="contig_id", how="left")
rvh = pd.read_csv(RVH)
idc = rvh.columns[0]
if "Unnamed" in idc or "y|virus order" in rvh.columns:
    rvh = rvh.rename(columns={idc: "contig_id"})
df = df.merge(rvh[[c for c in ["contig_id", "pred|L1", "pred|L2", "evidence"] if c in rvh.columns]],
              on="contig_id", how="left")
pb2 = pd.read_csv(PB2, sep="\t").rename(columns={"Accession": "contig_id"})
df = df.merge(pb2[[c for c in ["contig_id", "Host", "Host_NCBI_lineage", "Host_GTDB_lineage"]
                   if c in pb2.columns]], on="contig_id", how="left")

FALLBACK = set(rhp.NON_PLANT_FAMILIES_FALLBACK)
GENERA = set(rhp.NON_PLANT_GENERA)
TRUSTED = set(rhp.TRUSTED_LEVELS)


ORIG_BL = rhp.is_blacklisted          # 先抓住原始函数对象, 避免替换后自递归


def bl_now(fam, gen, det):
    return ORIG_BL(fam, gen, det)


def bl_v2(fam, gen, det):
    """科级否决优先: 科在名单内 → 直接否决, 不受判定层与属影响。"""
    f = "" if fam is None or pd.isna(fam) else str(fam).strip()
    if f and f in FALLBACK:
        return True
    if det is not None:
        d = str(det).split("(")[0].strip().lower()
        if d in TRUSTED:
            return False
    g = "" if gen is None or pd.isna(gen) else str(gen).strip()
    return bool(g and g not in ("NA", "nan") and g in GENERA)


def run(bl_fn, tag):
    orig = rhp.is_blacklisted
    rhp.is_blacklisted = bl_fn
    try:
        res = df.apply(rhp.decision_tree_cascade, axis=1)
    finally:
        rhp.is_blacklisted = orig
    h = [r[0] for r in res]
    n_bl = sum(1 for _, r in df.iterrows()
               if bl_fn(r.get("Family"), r.get("Genus"), r.get("Determination_Level")))
    print("  [%s] 否决命中 %d 行 | Plant=%d | %s"
          % (tag, n_bl, h.count("Plant"), dict(sorted(Counter(h).items(), key=lambda x: -x[1]))))
    return h


print("=== goji: 现规则 vs 候选 V2(科级否决优先) ===")
h_now = run(bl_now, "现规则")
h_v2 = run(bl_v2, "候选 V2")

flip = [(df.iloc[i]["contig_id"], df.iloc[i].get("Family"), df.iloc[i].get("Genus"),
         df.iloc[i].get("Species"), df.iloc[i].get("Determination_Level"), h_now[i], h_v2[i])
        for i in range(len(df)) if h_now[i] != h_v2[i]]
print("\n翻转行数 %d (其中 Plant→非Plant %d):"
      % (len(flip), sum(1 for f in flip if f[5] == "Plant")))
for f in flip:
    print("   %s | %s/%s | %s | %s | %s -> %s" % f)
