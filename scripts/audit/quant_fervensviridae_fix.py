#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""1) onekp 05 表两种拼写计数; 2) 用官方名 Fervensviridae 替换错拼条目对 goji Plant 的影响。"""
import importlib.util
import os
from collections import Counter

import pandas as pd

P = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

ON = "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"
if os.path.isfile(ON):
    d = pd.read_csv(ON, sep="\t", usecols=["Family"], low_memory=False)
    c = Counter(d["Family"].fillna("NA").astype(str))
    print("onekp 05 表 %d 行: Fervensiviridae=%d ; Fervensviridae=%d ; Microviridae=%d ; Autographiviridae=%d"
          % (len(d), c.get("Fervensiviridae", 0), c.get("Fervensviridae", 0),
             c.get("Microviridae", 0), c.get("Autographiviridae", 0)))
else:
    print("onekp 05 表不存在: %s" % ON)

TREE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
TAX = os.path.join(TREE, "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")
OUT = os.path.join(TREE, "06_HostPrediction")
df = pd.read_csv(TAX, sep="\t")
c9 = pd.read_csv(os.path.join(OUT, "C9_ICTV_result/classification_result.tsv"), sep="\t")[
    ["contig_id", "Predicted_Host", "Determination_Level"]].rename(columns={"Predicted_Host": "Host_ICTV"})
df = df.merge(c9, on="contig_id", how="left")
rvh = pd.read_csv(os.path.join(OUT, "RVH_result/result.csv"))
idc = rvh.columns[0]
if "Unnamed" in idc or "y|virus order" in rvh.columns:
    rvh = rvh.rename(columns={idc: "contig_id"})
df = df.merge(rvh[[c for c in ["contig_id", "pred|L1", "pred|L2", "evidence"] if c in rvh.columns]],
              on="contig_id", how="left")
pb2 = pd.read_csv(os.path.join(OUT, "phabox2_output/final_prediction/cherry_prediction.tsv"),
                  sep="\t").rename(columns={"Accession": "contig_id"})
df = df.merge(pb2[[c for c in ["contig_id", "Host", "Host_NCBI_lineage", "Host_GTDB_lineage"]
                   if c in pb2.columns]], on="contig_id", how="left")

BASE = list(m.NON_PLANT_FAMILIES_FALLBACK)
print("\n现名单: Fervensiviridae in-list=%s ; Fervensviridae in-list=%s"
      % ("Fervensiviridae" in BASE, "Fervensviridae" in BASE))


def run(tag, lst):
    m.NON_PLANT_FAMILIES_FALLBACK[:] = lst
    res = df.apply(m.decision_tree_cascade, axis=1)
    h = [r[0] for r in res]
    m.NON_PLANT_FAMILIES_FALLBACK[:] = BASE
    print("  [%s] 科 %d | Plant=%d" % (tag, len(lst), h.count("Plant")))
    return h


h0 = run("现状(含错拼)", BASE)
FIX = [f for f in BASE if f != "Fervensiviridae"] + ["Fervensviridae"]
h1 = run("改官方名", FIX)
d = [i for i in range(len(df)) if h0[i] != h1[i]]
print("  翻转 %d 行: %s" % (len(d), [(df.iloc[i]["contig_id"], df.iloc[i].get("Family"),
                                 df.iloc[i].get("Genus"), df.iloc[i].get("Determination_Level"),
                                 h0[i], h1[i]) for i in d][:10]))
fam_rows = df[df["Family"] == "Fervensviridae"]
print("  Fervensviridae 行: %d ; 其中 C9 判 Plant %d ; det 分布 %s"
      % (len(fam_rows), int((fam_rows["Host_ICTV"].astype(str).str.strip() == "Plant").sum()),
         dict(Counter(fam_rows["Determination_Level"].fillna("NA")))))
