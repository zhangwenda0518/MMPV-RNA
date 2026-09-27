#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Sylvanvirus / Crucivirus 等 7 个 "VMR 无记录" 属的溯源与翻转行构成。
只读。用法: python3 /tmp/trace_novmr_genera.py
"""
import importlib.machinery
import importlib.util
import os
import subprocess

import pandas as pd

PIPE = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline"
CUR = os.path.join(PIPE, "run_host_prediction.py")
TREES = {
    "goji": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out",
    "onekp": "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus",
}
TARGET = ["Sylvanvirus", "Crucivirus", "Klosneuvirus", "Hokovirus", "Indivirus",
          "Clandestinovirus", "Arenavirus"]


def load_mod(path, name):
    loader = importlib.machinery.SourceFileLoader(name, path)
    spec = importlib.util.spec_from_loader(name, loader)
    m = importlib.util.module_from_spec(spec)
    loader.exec_module(m)
    return m


def load_frame(tree):
    tax = os.path.join(tree, "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")
    out = os.path.join(tree, "06_HostPrediction")
    df = pd.read_csv(tax, sep="\t")
    c9 = pd.read_csv(os.path.join(out, "C9_ICTV_result/classification_result.tsv"), sep="\t")
    c9 = c9[[c for c in ["contig_id", "Predicted_Host", "Determination_Level"] if c in c9.columns]]
    df = df.merge(c9.rename(columns={"Predicted_Host": "Host_ICTV"}), on="contig_id", how="left")
    rvp = os.path.join(out, "RVH_result/result.csv")
    if os.path.isfile(rvp):
        rvh = pd.read_csv(rvp)
        idc = rvh.columns[0]
        if "Unnamed" in idc:
            rvh = rvh.rename(columns={idc: "contig_id"})
        df = df.merge(rvh[[c for c in ["contig_id", "pred|L1", "pred|L2", "evidence"]
                           if c in rvh.columns]], on="contig_id", how="left")
    pbp = os.path.join(out, "phabox2_output/final_prediction/cherry_prediction.tsv")
    if os.path.isfile(pbp):
        pb = pd.read_csv(pbp, sep="\t").rename(columns={"Accession": "contig_id"})
        df = df.merge(pb[[c for c in ["contig_id", "Host", "Host_NCBI_lineage"]
                          if c in pb.columns]], on="contig_id", how="left")
    return df


print("=" * 92)
print("A. Plant.tsv 里是否出现这些属名 (溯源用)")
print("=" * 92)
PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"
hits = {t: 0 for t in TARGET}
with open(PLANT, encoding="utf-8", errors="replace") as fh:
    for line in fh:
        for t in TARGET:
            if t in line:
                hits[t] += 1
for t in TARGET:
    print("  %-18s Plant.tsv 命中行 %d" % (t, hits[t]))

print()
print("=" * 92)
print("B. 两棵树内这些属的行构成与翻转行所在科")
print("=" * 92)
m0 = load_mod(CUR, "rhp0")
for tag, tree in TREES.items():
    df = load_frame(tree)
    base = [r[0] for r in df.apply(m0.decision_tree_cascade, axis=1)]
    gen_s = df["Genus"].astype(str).str.strip()
    print()
    print("  ==== %s" % tag)
    for g in TARGET:
        sub = df[gen_s == g]
        if not len(sub):
            continue
        m2 = load_mod(CUR, "rhp_t_%s" % g)
        if g in m2.NON_PLANT_GENERA:
            m2.NON_PLANT_GENERA.remove(g)
        sub2 = sub.copy()
        sub2["_old"] = [base[i] for i in sub.index]
        sub2["_new"] = [r[0] for r in sub.apply(m2.decision_tree_cascade, axis=1)]
        fl = sub2[sub2["_old"] != sub2["_new"]]
        print("   [%s] 行 %d, 翻转 %d" % (g, len(sub), len(fl)))
        print("       翻转行 科/属/判定级/Host_ICTV/旧->新:")
        if len(fl):
            print(fl.groupby(["Family", "Genus", "Determination_Level", "Host_ICTV",
                              "_old", "_new"]).size()
                  .sort_values(ascending=False).head(8).to_string())
        print("       全部行 科分布 top5: %s" % {str(k): int(v) for k, v in
              sub["Family"].value_counts().head(5).items()})
print()
print("DONE")
