#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""残留洞: 被否决的 C9 Plant 行经 RVH/PB2 分支重新落回 Plant。
测 A: 这批行明细; B: 若把科级否决扩展到 RVH/PB2 分支, Plant 行数变化。
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
OUT = os.path.join(TREE, "06_HostPrediction")
C9_TREE = os.path.join(OUT, "C9_ICTV_result/classification_result.tsv")
C9_GATE = "/tmp/c9fix_rehearsal/C9_gate/classification_result.tsv"
RVH = os.path.join(OUT, "RVH_result/result.csv")
PB2 = os.path.join(OUT, "phabox2_output/final_prediction/cherry_prediction.tsv")


def build(tax_path, c9_path):
    df = pd.read_csv(tax_path, sep="\t")
    cols = ["contig_id", "Predicted_Host"] + [c for c in ("Determination_Level",)
                                             if c in pd.read_csv(c9_path, sep="\t", nrows=0).columns]
    c9 = pd.read_csv(c9_path, sep="\t")[cols].rename(columns={"Predicted_Host": "Host_ICTV"})
    df = df.merge(c9, on="contig_id", how="left")
    rvh = pd.read_csv(RVH)
    idc = rvh.columns[0]
    if "Unnamed" in idc or "y|virus order" in rvh.columns:
        rvh = rvh.rename(columns={idc: "contig_id"})
    df = df.merge(rvh[[c for c in ["contig_id", "pred|L1", "pred|L2", "evidence"] if c in rvh.columns]],
                  on="contig_id", how="left")
    pb2 = pd.read_csv(PB2, sep="\t").rename(columns={"Accession": "contig_id"})
    return df.merge(pb2[[c for c in ["contig_id", "Host", "Host_NCBI_lineage", "Host_GTDB_lineage"]
                         if c in pb2.columns]], on="contig_id", how="left")


def tree_all_branch(row):
    """科级否决适用于所有分支: 非植物科的行不允许被任何来源判为 Plant。"""
    if str(row.get("Class", "")).lower() in rhp.PHAGE_CLASSES:
        return "Bacteria", "Rule_Class_Taxonomy"
    fam, gen, det = rhp.FAMILY_FIRST_VETO and row.get("Family"), row.get("Genus"), row.get("Determination_Level")
    veto = rhp.is_blacklisted(fam, gen, det)
    h = rhp.normalize_c9(row.get("Host_ICTV", "Unknown"))
    if h != "Unknown":
        if not (h == "Plant" and veto):
            return h, "ICTV_Preferred"
        h = "Unknown"
    r = rhp.parse_rvh(row)
    if r != "Unknown":
        if not (r == "Plant" and veto):
            return r, "RVH_Preferred"
    p = rhp.parse_pb2(row)
    if p != "Unknown":
        if not (p == "Plant" and veto):
            return p, "PB2_Preferred"
    return "Unknown", "Unassigned"


for tag, tax, c9p in [("未闸门", os.path.join(TREE, "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"), C9_TREE),
                      ("闸门版", "/tmp/v62d_goji/final_integrated_classification.tsv", C9_GATE)]:
    df = build(tax, c9p)
    res_now = df.apply(rhp.decision_tree_cascade, axis=1)
    res_new = df.apply(tree_all_branch, axis=1)
    now = [r[0] for r in res_now]
    new = [r[0] for r in res_new]
    print("=== %s ===" % tag)
    print("  现部署 Plant=%d ; 扩展否决到 RVH/PB2 后 Plant=%d" % (now.count("Plant"), new.count("Plant")))
    idx = [i for i in range(len(df)) if now[i] != new[i]]
    print("  差异 %d 行:" % len(idx))
    for i in idx:
        r = df.iloc[i]
        print("    %s | %s/%s/%s | det=%s | pred|L1=%s | PB2=%s | %s -> %s"
              % (r["contig_id"], r.get("Family"), r.get("Genus"), r.get("Species"),
                 r.get("Determination_Level"), r.get("pred|L1"), r.get("Host"),
                 res_now[i], res_new[i]))
