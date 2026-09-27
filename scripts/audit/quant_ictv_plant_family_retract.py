#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ICTV 复核发现的"宿主含植物"科, 从非植物科名单撤回后的影响 (goji 树实测, 内存模拟)。"""
import importlib.util
import os
from collections import Counter

import pandas as pd

RHP = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", RHP)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

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

PLANT_INCLUDED = ["Artoviridae", "Chrysoviridae", "Genomoviridae", "Pestiviridae",
                  "Pseudoviridae", "Spiciviridae", "Tomosaviridae", "Metaviridae",
                  "Kanorauviridae", "Mitoviridae", "Ourmiaviridae"]
BASE_F = list(m.NON_PLANT_FAMILIES_FALLBACK)


def run(tag, drop=()):
    m.NON_PLANT_FAMILIES_FALLBACK[:] = [f for f in BASE_F if f not in drop]
    res = df.apply(m.decision_tree_cascade, axis=1)
    h = [r[0] for r in res]
    m.NON_PLANT_FAMILIES_FALLBACK[:] = BASE_F
    print("  [%s] 科名单 %d | Plant=%d" % (tag, len(m.NON_PLANT_FAMILIES_FALLBACK) - len(drop), h.count("Plant")))
    return h


print("当前名单 %d 科" % len(BASE_F))
print("待处置科在各表中的出现情况:")
for f in PLANT_INCLUDED:
    n05 = int((df["Family"] == f).sum())
    c9p = int(((df["Family"] == f) & (df["Host_ICTV"].astype(str).str.strip() == "Plant")).sum())
    print("   %-18s in-list=%-5s 05 表 %5d 行 | C9 判 Plant %3d 行"
          % (f, f in BASE_F, n05, c9p))

h0 = run("现状")
h1 = run("撤回宿主含植物的 6 科",
         [f for f in ["Artoviridae", "Chrysoviridae", "Genomoviridae", "Pestiviridae",
                      "Pseudoviridae", "Spiciviridae"] if f in BASE_F])
d = [i for i in range(len(df)) if h0[i] != h1[i]]
print("  翻转 %d 行:" % len(d))
for i in d:
    r = df.iloc[i]
    print("    %s | %s/%s/%s | det=%s | %s -> %s"
          % (r["contig_id"], r.get("Family"), r.get("Genus"), r.get("Species"),
             r.get("Determination_Level"), h0[i], h1[i]))
