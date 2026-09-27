#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""验证 C9 字段碰撞修复 + 预演 05 闸门表进入 host 决策后的结果。

Part 1: 现网(未闸门) 05 表 → 修复后黑名单命中数与 Plant 行数
Part 2: 闸门版 05 表 (/tmp/v62d_goji) 重跑 C9 → 三科嵌合行是否离开 Plant
"""
import importlib.util
import os
import subprocess
from collections import Counter

import pandas as pd

RHP = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
CG = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/utils/classify_contigs.py"
PROB = "/home/zhangwenda/MMPV-RNA/database/cross_analysis"
spec = importlib.util.spec_from_file_location("rhp", RHP)
rhp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rhp)

TREE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
TAX_OLD = os.path.join(TREE, "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")
TAX_GATE = "/tmp/v62d_goji/final_integrated_classification.tsv"
OUT = os.path.join(TREE, "06_HostPrediction")
C9_OLD = os.path.join(OUT, "C9_ICTV_result/classification_result.tsv")
RVH = os.path.join(OUT, "RVH_result/result.csv")
PB2 = os.path.join(OUT, "phabox2_output/final_prediction/cherry_prediction.tsv")
WORK = "/tmp/c9fix_rehearsal"
os.makedirs(WORK, exist_ok=True)


def build_merge(tax_path, c9_path):
    """按修复后的 run_ensemble 逻辑构建合并表。"""
    df = pd.read_csv(tax_path, sep="\t")
    if c9_path and os.path.isfile(c9_path):
        cols = ["contig_id", "Predicted_Host"]
        cols += [c for c in ("Determination_Level",)
                 if c in pd.read_csv(c9_path, sep="\t", nrows=0).columns]
        c9 = pd.read_csv(c9_path, sep="\t")[cols].rename(columns={"Predicted_Host": "Host_ICTV"})
        df = df.merge(c9, on="contig_id", how="left")
    if os.path.isfile(RVH):
        rvh = pd.read_csv(RVH)
        idc = rvh.columns[0]
        if "Unnamed" in idc or "y|virus order" in rvh.columns:
            rvh = rvh.rename(columns={idc: "contig_id"})
        use = ["contig_id"] + [c for c in ["pred|L1", "pred|L2", "evidence"] if c in rvh.columns]
        df = df.merge(rvh[use], on="contig_id", how="left")
    if os.path.isfile(PB2):
        pb2 = pd.read_csv(PB2, sep="\t").rename(columns={"Accession": "contig_id"})
        use = [c for c in ["contig_id", "Host", "Host_NCBI_lineage", "Host_GTDB_lineage"]
               if c in pb2.columns]
        df = df.merge(pb2[use], on="contig_id", how="left")
    return df


def run_tree(df, tag):
    res = df.apply(rhp.decision_tree_cascade, axis=1)
    hosts = [r[0] for r in res]
    methods = [r[1] for r in res]
    n_bl = sum(1 for _, row in df.iterrows()
               if rhp.is_blacklisted(row.get("Family"), row.get("Genus"),
                                     row.get("Determination_Level")))
    print("  [%s] Family 列在合并表中: %s ; 黑名单命中行数: %d"
          % (tag, "Family" in df.columns, n_bl))
    print("        Final_Host: %s" % dict(Counter(hosts)))
    return hosts, methods


print("=== Part 1: 未闸门 05 表 (现网树) ===")
m1 = build_merge(TAX_OLD, C9_OLD)
h1, me1 = run_tree(m1, "old-tax")

print("\n=== Part 2: 闸门版 05 表重跑 C9 ===")
gate_dir = os.path.join(WORK, "C9_gate")
os.makedirs(gate_dir, exist_ok=True)
c9_new = os.path.join(gate_dir, "classification_result.tsv")
if not os.path.isfile(c9_new):
    cmd = ["python3", CG, "-i", TAX_GATE, "--output_dir", gate_dir,
           "--prob_dir", PROB, "--mode", "high"]
    print("  $ %s" % " ".join(cmd))
    p = subprocess.run(cmd, capture_output=True, text=True)
    print("  rc=%d" % p.returncode)
    print("  " + "\n  ".join((p.stdout or "").strip().split("\n")[-6:]))
    if p.returncode != 0:
        print("  stderr: %s" % (p.stderr or "")[-800:])

if os.path.isfile(c9_new):
    m2 = build_merge(TAX_GATE, c9_new)
    h2, me2 = run_tree(m2, "gate-tax")
    # 三科嵌合行逐条对照
    key = ["contig_id", "Family", "Genus", "Species", "Determination_Level"]
    old = pd.read_csv(C9_OLD, sep="\t")
    new = pd.read_csv(c9_new, sep="\t")
    old_keep = [c for c in key if c in old.columns]
    new_keep = [c for c in key if c in new.columns]
    o = old[old_keep].add_suffix("_old").rename(columns={"contig_id_old": "contig_id"})
    n = new[new_keep].add_suffix("_new").rename(columns={"contig_id_new": "contig_id"})
    cm = o.merge(n, on="contig_id", how="inner")
    cm["host_old"] = h1
    cm["host_new"] = h2
    tgt = cm[cm["Family_old"].isin(["Mimiviridae", "Marseilleviridae", "Pithoviridae"])
             | cm["Family_new"].isin(["Mimiviridae", "Marseilleviridae", "Pithoviridae"])]
    print("\n  三科行数(合并表内): %d" % len(tgt))
    print("  其中 Plant→非Plant: %d ; 仍为 Plant: %d"
          % (int(((tgt["host_old"] == "Plant") & (tgt["host_new"] != "Plant")).sum()),
             int((tgt["host_new"] == "Plant").sum())))
    cols = ["contig_id", "Family_old", "Genus_old", "Determination_Level_old", "host_old",
            "Family_new", "Genus_new", "Determination_Level_new", "host_new"]
    cols = [c for c in cols if c in tgt.columns]
    with pd.option_context("display.width", 250, "display.max_colwidth", 46):
        print(tgt[cols].to_string(index=False))
