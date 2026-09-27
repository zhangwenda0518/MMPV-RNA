#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""对拍三种口径, 直接回答"属级证据可信就直接采信"到底会放行哪些行。

V0 = 现状           : FAMILY_FIRST_VETO=True,  WHITELIST_OVERRIDES_FAMILY_VETO=False
V1 = 白名单越权(W3) : FAMILY_FIRST_VETO=True,  WHITELIST_OVERRIDES_FAMILY_VETO=True
V2 = 属级优先(全放) : FAMILY_FIRST_VETO=False (属/种级判定直接采信 C9, 科名单只兜底属缺失)

输出: 每棵树三种口径的 Plant 计数、两两差集规模, 以及 V2-V1 新增行的
      (Family, Genus, Determination_Level) 分布与是否命中植物属白名单。
结果写 /tmp/three_modes.log
"""
import importlib.machinery
import importlib.util
import os
import subprocess
import sys

import pandas as pd

PIPE = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline"
CUR = os.path.join(PIPE, "run_host_prediction.py")
TREES = {
    "goji": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out",
    "onekp": "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus",
}


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


def run_modes(df, m):
    res = {}
    for tag, ff, w3 in [("V0", True, False), ("V1", True, True), ("V2", False, False)]:
        m.FAMILY_FIRST_VETO = ff
        m.WHITELIST_OVERRIDES_FAMILY_VETO = w3
        res[tag] = [r[0] for r in df.apply(m.decision_tree_cascade, axis=1)]
        print("   %s done" % tag, flush=True)
    m.FAMILY_FIRST_VETO = True
    m.WHITELIST_OVERRIDES_FAMILY_VETO = False
    return res


def main():
    m = load_mod(CUR, "rhp")
    for tag, tree in TREES.items():
        print("\n================ %s ================" % tag, flush=True)
        df = load_frame(tree)
        print("行 %d" % len(df), flush=True)
        res = run_modes(df, m)
        for k in ("V0", "V1", "V2"):
            print("%s Plant=%d 非Plant=%d" % (k, res[k].count("Plant"), len(res[k]) - res[k].count("Plant")))
        for a, b in (("V0", "V1"), ("V1", "V2"), ("V0", "V2")):
            idx = [i for i in range(len(df)) if res[a][i] != res[b][i]]
            up = sum(1 for i in idx if res[a][i] != "Plant" and res[b][i] == "Plant")
            print("%s -> %s : 变化行 %d (升为 Plant %d, 其他 %d)" % (a, b, len(idx), up, len(idx) - up))
        # V1 vs V2 明细: V2 放行而 V1 否决的行
        idx2 = [i for i in range(len(df)) if res["V1"][i] != "V2" and res["V2"][i] == "Plant"]
        print("V2 比 V1 多判 Plant 的行: %d" % len(idx2))
        if idx2:
            sub = df.iloc[idx2].copy()
            sub["_V1"] = [res["V1"][i] for i in idx2]
            sub["_V2"] = [res["V2"][i] for i in idx2]
            wl = [i for i in idx2 if str(df.iloc[i].get("Genus")).strip() in m.PLANT_GENERA_WHITELIST]
            print("   其中属命中植物白名单: %d ; 未命中: %d" % (len(wl), len(idx2) - len(wl)))
            g = sub.groupby(["Family", "Genus", "Determination_Level"]).size().sort_values(ascending=False)
            print("   前 25 个 (Family, Genus, Determination_Level):")
            for (f, gg, d), n in g.head(25).items():
                print("     %-28s %-22s %-14s %5d" % (str(f)[:28], str(gg)[:22], str(d)[:14], n))
            cols = [c for c in ["contig_id", "Family", "Genus", "Determination_Level",
                                "Host_ICTV", "pred|L1", "pred|L2", "Host"] if c in sub.columns]
            sub[cols].to_csv("/tmp/modes_v2minusv1_%s.tsv" % tag, sep="\t", index=False)
            print("   明细写 /tmp/modes_v2minusv1_%s.tsv" % tag)
            print(sub[cols].head(8).to_string(max_colwidth=26))


with open("/tmp/three_modes.log", "w", encoding="utf-8") as fh:
    old = sys.stdout
    sys.stdout = fh
    try:
        main()
    finally:
        sys.stdout = old
print(open("/tmp/three_modes.log", encoding="utf-8").read())
