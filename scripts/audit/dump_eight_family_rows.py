#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""列出 8 个非植物科最终仍落在 Plant 的行及其证据链 (修复 C9 碰撞后的口径)。"""
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
FAMS = ["Metaviridae", "Draupnirviridae", "Iflaviridae", "Barnaviridae", "Discoviridae",
        "Hydriviridae", "Kanorauviridae", "Narnaviridae"]

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
res = df.apply(rhp.decision_tree_cascade, axis=1)
df["Final_Host"] = [r[0] for r in res]
df["Decision_Method"] = [r[1] for r in res]

sub = df[df["Family"].isin(FAMS)]
print("八科 contig 总数 %d → 最终 Plant %d"
      % (len(sub), int((sub["Final_Host"] == "Plant").sum())))
print("八科整体 Final_Host 分布: %s" % dict(Counter(sub["Final_Host"])))
print("\n按科统计 (总数 / 落 Plant 数 / 该科 Plant 行的判定层):")
for f in FAMS:
    s = sub[sub["Family"] == f]
    p = s[s["Final_Host"] == "Plant"]
    print("   %-18s %4d / %2d | %s" % (f, len(s), len(p),
                                      dict(Counter(p["Determination_Level"].fillna("NA")))))

cols = ["contig_id", "Family", "Genus", "Species", "Determination_Level", "Host_ICTV",
        "pred|L1", "evidence", "Host", "Decision_Method"]
cols = [c for c in cols if c in df.columns]
with pd.option_context("display.width", 260, "display.max_colwidth", 34):
    print("\n八科中最终为 Plant 的行:")
    print(sub[sub["Final_Host"] == "Plant"][cols].to_string(index=False))
