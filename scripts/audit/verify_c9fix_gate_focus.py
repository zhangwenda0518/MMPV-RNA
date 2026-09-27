#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""聚焦输出: 修复 C9 字段碰撞后 + 05 闸门表进入 host 决策的净效果 (只打摘要与三科 Plant 行)。"""
import importlib.util
import os
from collections import Counter

import pandas as pd

RHP = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", RHP)
rhp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rhp)

TREE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
TAX_OLD = os.path.join(TREE, "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")
TAX_GATE = "/tmp/v62d_goji/final_integrated_classification.tsv"
OUT = os.path.join(TREE, "06_HostPrediction")
C9_OLD = os.path.join(OUT, "C9_ICTV_result/classification_result.tsv")
C9_NEW = "/tmp/c9fix_rehearsal/C9_gate/classification_result.tsv"
RVH = os.path.join(OUT, "RVH_result/result.csv")
PB2 = os.path.join(OUT, "phabox2_output/final_prediction/cherry_prediction.tsv")
GFS = ["Mimiviridae", "Marseilleviridae", "Pithoviridae"]


def build_merge(tax_path, c9_path):
    df = pd.read_csv(tax_path, sep="\t")
    cols = ["contig_id", "Predicted_Host"]
    cols += [c for c in ("Determination_Level",)
             if c in pd.read_csv(c9_path, sep="\t", nrows=0).columns]
    c9 = pd.read_csv(c9_path, sep="\t")[cols].rename(columns={"Predicted_Host": "Host_ICTV"})
    df = df.merge(c9, on="contig_id", how="left")
    rvh = pd.read_csv(RVH)
    idc = rvh.columns[0]
    if "Unnamed" in idc or "y|virus order" in rvh.columns:
        rvh = rvh.rename(columns={idc: "contig_id"})
    use = ["contig_id"] + [c for c in ["pred|L1", "pred|L2", "evidence"] if c in rvh.columns]
    df = df.merge(rvh[use], on="contig_id", how="left")
    pb2 = pd.read_csv(PB2, sep="\t").rename(columns={"Accession": "contig_id"})
    use = [c for c in ["contig_id", "Host", "Host_NCBI_lineage", "Host_GTDB_lineage"]
           if c in pb2.columns]
    return df.merge(pb2[use], on="contig_id", how="left")


def run_tree(df, tag):
    res = df.apply(rhp.decision_tree_cascade, axis=1)
    h = [r[0] for r in res]
    n_bl = sum(1 for _, r in df.iterrows()
               if rhp.is_blacklisted(r.get("Family"), r.get("Genus"),
                                     r.get("Determination_Level")))
    hit_plant = sum(1 for i, r in enumerate(df.itertuples(index=False))
                    if rhp.is_blacklisted(df.iloc[i].get("Family"), df.iloc[i].get("Genus"),
                                          df.iloc[i].get("Determination_Level"))
                    and str(df.iloc[i].get("Host_ICTV")).strip() == "Plant")
    print("  [%s] 合并表含 Family 列=%s | 黑名单命中 %d 行 (其中 Host_ICTV=Plant 的 %d 行)"
          % (tag, "Family" in df.columns, n_bl, hit_plant))
    print("       Final_Host: %s" % dict(sorted(Counter(h).items(), key=lambda x: -x[1])))
    return h


print("=== Part 1  未闸门 05 表 (现网树) + 修复后的 run_ensemble ===")
h1 = run_tree(build_merge(TAX_OLD, C9_OLD), "old-tax")

print("\n=== Part 2  闸门版 05 表 (/tmp/v62d_goji) 重跑 C9 ===")
h2 = run_tree(build_merge(TAX_GATE, C9_NEW), "gate-tax")

o = pd.read_csv(C9_OLD, sep="\t")
n = pd.read_csv(C9_NEW, sep="\t")
key = ["contig_id", "Family", "Genus", "Species", "Determination_Level", "Predicted_Host"]
o = o[[c for c in key if c in o.columns]].add_suffix("_old").rename(columns={"contig_id_old": "contig_id"})
n = n[[c for c in key if c in n.columns]].add_suffix("_new").rename(columns={"contig_id_new": "contig_id"})
cm = o.merge(n, on="contig_id", how="inner")
cm["host_old"], cm["host_new"] = h1, h2
print("\n  三科 contig 总数: %d" % int(cm["Family_old"].isin(GFS).sum()))
tgt = cm[cm["Family_old"].isin(GFS) & (cm["host_old"] == "Plant")]
print("  其中旧表判定 Plant 的 %d 行, 闸门表后:" % len(tgt))
print("     " + str(dict(Counter(tgt["host_new"]))))
cols = ["contig_id", "Genus_old", "Species_old", "Determination_Level_old", "host_old",
        "Genus_new", "Determination_Level_new", "host_new"]
with pd.option_context("display.width", 220, "display.max_colwidth", 40):
    print(tgt[cols].to_string(index=False))
