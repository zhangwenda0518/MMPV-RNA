#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用服务器上现装的 run_host_prediction.py (C9 碰撞修复 + 名单核定 + 科级否决优先) 实测:
   未闸门 05 表 vs 闸门版 05 表, 各自的 Plant 行数与三科嵌合行去向。
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
GIANT = ["Mimiviridae", "Marseilleviridae", "Pithoviridae"]
print("部署态: 科名单 %d | 属名单 %d | FAMILY_FIRST_VETO=%s"
      % (len(rhp.NON_PLANT_FAMILIES_FALLBACK), len(rhp.NON_PLANT_GENERA),
         getattr(rhp, "FAMILY_FIRST_VETO", "N/A")))


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


for tag, tax, c9p in [("未闸门 05 表", os.path.join(TREE, "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"), C9_TREE),
                      ("闸门版 05 表", "/tmp/v62d_goji/final_integrated_classification.tsv", C9_GATE)]:
    df = build(tax, c9p)
    res = df.apply(rhp.decision_tree_cascade, axis=1)
    df["Final_Host"] = [r[0] for r in res]
    df["Decision_Method"] = [r[1] for r in res]
    print("\n=== %s ===" % tag)
    print("  Final_Host: %s" % dict(sorted(Counter(df["Final_Host"]).items(), key=lambda x: -x[1])))
    veto = df[[bool(rhp.is_blacklisted(r.Family, r.Genus, r.Determination_Level))
               for r in df.itertuples(index=False)]]
    plant_veto = veto[veto["Host_ICTV"].astype(str).str.strip() == "Plant"]
    print("  被科/属否决的行 %d, 其中 C9 原判 Plant 的 %d 行 → 去向 %s"
          % (len(veto), len(plant_veto), dict(Counter(plant_veto["Final_Host"]))))
    print("  否决涉及的科 Top8: %s"
          % dict(Counter(plant_veto["Family"].fillna("NA")).most_common(8)))
    g = df[df["Family"].isin(GIANT)]
    print("  三科 contig %d → 最终 Plant %d" % (len(g), int((g["Final_Host"] == "Plant").sum())))
    stay = g[g["Final_Host"] == "Plant"]
    if len(stay):
        print("    仍为 Plant 的行:")
        for r in stay.itertuples(index=False):
            print("      %s | %s/%s/%s | %s" % (r.contig_id, r.Family, r.Genus, r.Species,
                                                r.Decision_Method))
