#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""争议科 (VMR 宿主列含 plants 的科兜底条目) 在两棵树的行数与剥离后 Plant 增量。
只读。用法: python3 /tmp/measure_family_retract_contended.py
"""
import importlib.machinery
import importlib.util
import os

import pandas as pd

PIPE = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline"
CUR = os.path.join(PIPE, "run_host_prediction.py")
TREES = {
    "goji": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out",
    "onekp": "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus",
}
FAMS = ["Ouroboviridae", "Discoviridae", "Ambiguiviridae", "Steitzviridae", "Fiersviridae"]


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


m0 = load_mod(CUR, "rhp_f")
for tag, tree in TREES.items():
    df = load_frame(tree)
    base = [r[0] for r in df.apply(m0.decision_tree_cascade, axis=1)]
    print()
    print("  ==== %s (Plant %d / %d 行)" % (tag, base.count("Plant"), len(df)))
    fam_s = df["Family"].astype(str).str.strip()
    for f in FAMS:
        sub = df[fam_s == f]
        if not len(sub):
            print("     %-18s 树内 0 行 (不在名单或未被标到)" % f)
            continue
        inlist = f in m0.NON_PLANT_FAMILIES_FALLBACK
        m2 = load_mod(CUR, "rhp_f_%s" % f)
        if f in m2.NON_PLANT_FAMILIES_FALLBACK:
            m2.NON_PLANT_FAMILIES_FALLBACK.remove(f)
        old = [base[i] for i in sub.index]
        new = [r[0] for r in sub.apply(m2.decision_tree_cascade, axis=1)]
        fl = sum(1 for a, b in zip(old, new) if a != b)
        top = sum(1 for a, b in zip(old, new) if a != "Plant" and b == "Plant")
        print("     %-18s 在科兜底名单=%s, 树内 %6d 行, 剥离后翻转 %4d (->Plant %4d), 现判定 %s"
              % (f, inlist, len(sub), fl, top,
                 {str(k): int(v) for k, v in pd.Series(old).value_counts().items()}))
print()
print("DONE")
