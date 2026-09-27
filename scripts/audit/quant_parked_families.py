#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Kanorauviridae / Metaviridae 是否该留在 "非植物科" 兜底名单: 用现网数据量化其影响。

ICTV 现行宿主列: Metaviridae=原生生物/真菌/植物/无脊椎/脊椎; Kanorauviridae 表列植物
(与 ViralZone "宿主未知" 冲突)。二者本不该简单标为非植物科 → 测其被移除后的 Plant 变化。
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
PARKED = ["Kanorauviridae", "Metaviridae"]

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

FULL = list(rhp.NON_PLANT_FAMILIES_FALLBACK)
print("当前兜底名单: %d 科; 其中待复核科: %s"
      % (len(FULL), [f for f in FULL if f in PARKED]))


def hosts_with(fam_list):
    rhp.NON_PLANT_FAMILIES_FALLBACK[:] = fam_list
    res = df.apply(rhp.decision_tree_cascade, axis=1)
    return [r[0] for r in res]


h_full = hosts_with(FULL)
h_cut = hosts_with([f for f in FULL if f not in PARKED])
cfull, ccut = Counter(h_full), Counter(h_cut)
print("\n保留两科: Plant=%d ; 移除两科: Plant=%d" % (cfull.get("Plant", 0), ccut.get("Plant", 0)))
diff = [(df.iloc[i]["contig_id"], df.iloc[i].get("Family"), df.iloc[i].get("Genus"),
         df.iloc[i].get("Determination_Level"), h_full[i], h_cut[i])
        for i in range(len(df)) if h_full[i] != h_cut[i]]
print("受影响行数: %d" % len(diff))
for d in diff:
    print("   %s | %s/%s | %s | %s -> %s" % d)

# 兜底名单对 Plant 行的实际贡献 (修复 C9 碰撞后)
rhp.NON_PLANT_FAMILIES_FALLBACK[:] = FULL
res = df.apply(rhp.decision_tree_cascade, axis=1)
df2 = df.assign(Final_Host=[r[0] for r in res])
plant_rows = df2[df2["Final_Host"] == "Plant"]
print("\n修复后 Plant 行 %d ; 其 Family 分布 Top15:" % len(plant_rows))
print(plant_rows["Family"].fillna("NA").value_counts().head(15).to_string())
print("\n八科在 Plant 行中的占比:")
for f in ["Metaviridae", "Draupnirviridae", "Iflaviridae", "Barnaviridae", "Discoviridae",
          "Hydriviridae", "Kanorauviridae", "Narnaviridae"]:
    sub = plant_rows[plant_rows["Family"] == f]
    if len(sub):
        print("   %-18s %3d 行 | 例: %s" % (f, len(sub), list(sub["contig_id"].head(2))))
