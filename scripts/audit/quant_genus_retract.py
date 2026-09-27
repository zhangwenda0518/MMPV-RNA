#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""量化属黑名单撤回 (Betapartitivirus / Geminivirus) 对 goji Plant 行的影响 (内存模拟)。"""
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

BASE = list(m.NON_PLANT_GENERA)


def run(tag, drop=()):
    m.NON_PLANT_GENERA[:] = [g for g in BASE if g not in drop]
    res = df.apply(m.decision_tree_cascade, axis=1)
    h = [r[0] for r in res]
    print("  [%s] 属名单 %d | Plant=%d | %s"
          % (tag, len(m.NON_PLANT_GENERA), h.count("Plant"),
             dict(sorted(Counter(h).items(), key=lambda x: -x[1]))))
    m.NON_PLANT_GENERA[:] = BASE
    return h


h0 = run("现状")
for g in ["Betapartitivirus", "Geminivirus"]:
    idx = [i for i in range(len(df))
           if df.iloc[i].get("Genus") == g]
    print("  属 %s 涉及 %d 行; 其中 C9 判 Plant %d; 族分布 %s"
          % (g, len(idx), sum(1 for i in idx if str(df.iloc[i].get("Host_ICTV")).strip() == "Plant"),
             dict(Counter(df.iloc[i].get("Family") for i in idx).most_common(5))))
    h1 = run("撤回 " + g, (g,))
    d = [i for i in range(len(df)) if h0[i] != h1[i]]
    print("    翻转 %d 行: %s" % (len(d), [(df.iloc[i]["contig_id"], df.iloc[i].get("Family"),
                                        df.iloc[i].get("Genus"), df.iloc[i].get("Determination_Level"),
                                        h0[i], h1[i]) for i in d][:12]))
h2 = run("撤回两者", ("Betapartitivirus", "Geminivirus"))
print("  两者合计翻转 %d 行" % sum(1 for i in range(len(df)) if h0[i] != h2[i]))
