#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""候选属加入黑名单的影响实测 (goji + onekp, 若 06 产物在)。

判据 (保守, 只收"非植物无疑义"的属):
  属级宿主 ∈ {Insecta, Human, Mammalia, Aves, Animal_other, Algae, Protist, Bacteria}
  且 Plant.tsv 无该属记录
  排除 Fungi 宿主属 —— ICTV 把这 9 个科的宿主写成含 plants, 而我们的库 Plant.tsv 对这些科零覆盖,
  属级"Fungi"很可能只是库缺植物记录造成的假非植物 (例如 Mitovirus / Alphachrysovirus / Gemycircularvirus)。
"""
import importlib.util
import os
import sys
from collections import Counter

import pandas as pd

TRI = "/tmp/PLANT_FAMILY_GENUS_TRIAGE.tsv"
t = pd.read_csv(TRI, sep="\t")
OK_HOST = {"Insecta", "Human", "Mammalia", "Aves", "Animal_other", "Algae", "Protist", "Bacteria"}
cand = t[(t["Verdict"].str.startswith("E3")) & (~t["Current_Blacklist"]) &
         (t["Genus_HP_Predicted_Host"].isin(OK_HOST))].copy()
cand = cand.sort_values("Rows_onekp", ascending=False)
print("保守候选 %d 个:" % len(cand))
for r in cand.itertuples(index=False):
    print("   %-18s %-22s goji=%-4d onekp=%-4d %-12s 总=%-8d %s"
          % (r.Family, r.Genus, r.Rows_goji, r.Rows_onekp, r.Genus_HP_Predicted_Host,
             r.Genus_HP_Total_Records, r.Genus_HP_Confidence))
print("排除的 Fungi 宿主属: %s"
      % sorted(t[(t["Verdict"].str.startswith("E3")) & (~t["Current_Blacklist"]) &
                 (t["Genus_HP_Predicted_Host"] == "Fungi")]["Genus"]))

P = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
BASE_G = list(m.NON_PLANT_GENERA)
BASE_F = list(m.NON_PLANT_FAMILIES_FALLBACK)
ADD = sorted(set(cand["Genus"]))

TREES = {
    "goji": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out",
    "onekp": "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus",
}


def load(tree):
    tax = os.path.join(tree, "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")
    out = os.path.join(tree, "06_HostPrediction")
    if not (os.path.isfile(tax) and os.path.isdir(os.path.join(out, "C9_ICTV_result"))):
        return None
    df = pd.read_csv(tax, sep="\t")
    c9 = pd.read_csv(os.path.join(out, "C9_ICTV_result/classification_result.tsv"), sep="\t")
    keep = [c for c in ["contig_id", "Predicted_Host", "Determination_Level"] if c in c9.columns]
    c9 = c9[keep].rename(columns={"Predicted_Host": "Host_ICTV"})
    df = df.merge(c9, on="contig_id", how="left")
    rvp = os.path.join(out, "RVH_result/result.csv")
    if os.path.isfile(rvp):
        rvh = pd.read_csv(rvp)
        idc = rvh.columns[0]
        if "Unnamed" in idc or "y|virus order" in rvh.columns:
            rvh = rvh.rename(columns={idc: "contig_id"})
        cols = [c for c in ["contig_id", "pred|L1", "pred|L2", "evidence"] if c in rvh.columns]
        df = df.merge(rvh[cols], on="contig_id", how="left")
    pbp = os.path.join(out, "phabox2_output/final_prediction/cherry_prediction.tsv")
    if os.path.isfile(pbp):
        pb = pd.read_csv(pbp, sep="\t").rename(columns={"Accession": "contig_id"})
        cols = [c for c in ["contig_id", "Host", "Host_NCBI_lineage", "Host_GTDB_lineage"]
                if c in pb.columns]
        df = df.merge(pb[cols], on="contig_id", how="left")
    return df


def label(df):
    m.NON_PLANT_GENERA[:] = BASE_G
    m.NON_PLANT_FAMILIES_FALLBACK[:] = BASE_F
    res = df.apply(m.decision_tree_cascade, axis=1)
    return [r[0] for r in res]


for tag, tree in TREES.items():
    df = load(tree)
    print("\n" + "=" * 92)
    if df is None:
        print("[%s] 06 产物不全, 跳过" % tag)
        continue
    h0 = label(df)
    m.NON_PLANT_GENERA[:] = BASE_G + [g for g in ADD if g not in BASE_G]
    res1 = df.apply(m.decision_tree_cascade, axis=1)
    h1 = [r[0] for r in res1]
    m.NON_PLANT_GENERA[:] = BASE_G
    print("[%s] 行数 %d | 现状 Plant=%d | 加入 %d 属后 Plant=%d"
          % (tag, len(df), h0.count("Plant"), len(ADD), h1.count("Plant")))
    flips = [i for i in range(len(df)) if h0[i] != h1[i]]
    print("  翻转 %d 行" % len(flips))
    cnt = Counter()
    for i in flips[:100000]:
        r = df.iloc[i]
        cnt[(r.get("Family"), r.get("Genus"), r.get("Determination_Level"), str(h0[i]) + "->" + str(h1[i]))] += 1
    for k, v in cnt.most_common(25):
        print("     %-22s %-20s det=%-22s %s ×%d" % (k[0], k[1], k[2], k[3], v))
    lost = [i for i in flips if h0[i] == "Plant"]
    print("  其中原判 Plant 被撤回: %d" % len(lost))
    for i in lost[:15]:
        r = df.iloc[i]
        print("     %s | F=%s G=%s det=%s ICTV=%s RVH=%s PHB=%s"
              % (str(r.get("contig_id"))[:28], r.get("Family"), r.get("Genus"),
                 r.get("Determination_Level"), r.get("Host_ICTV"), r.get("pred|L1"), r.get("Host")))
