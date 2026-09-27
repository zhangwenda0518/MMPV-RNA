#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""三种否决口径的差集明细 (修正上一版字符串笔误: 上一版把 "V2" 当标签比, 导致明细表其实是全量 Plant)。
本版: 明确算 V1-V0 与 V2-V1 两个差集, 各自给 (科,属,判定层) 分布与前几行的工具输入,
并把三份标签存盘 /tmp/labels_<tree>_<mode>.tsv 以便复用。
"""
import importlib.machinery
import importlib.util
import os
import sys

import pandas as pd

PIPE = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline"
CUR = os.path.join(PIPE, "run_host_prediction.py")
TREES = {
    "goji": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out",
    "onekp": "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus",
}
COLS = ["contig_id", "Family", "Genus", "Determination_Level", "Host_ICTV",
        "pred|L1", "pred|L2", "evidence", "Host", "Host_NCBI_lineage"]


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


def report(df, labels, idx, tag, m, wl_check=True):
    if not idx:
        print("   (空差集)")
        return
    sub = df.iloc[idx].copy()
    if wl_check:
        inwl = [str(g).strip() in m.PLANT_GENERA_WHITELIST for g in sub.get("Genus")]
        print("   属命中植物白名单 %d / 未命中 %d" % (sum(inwl), len(inwl) - sum(inwl)))
    g = sub.groupby(["Family", "Genus", "Determination_Level"]).size().sort_values(ascending=False)
    print("   (科, 属, 判定层) 分布:")
    for (f, gg, d), n in g.items():
        print("     %-26s %-22s %-20s %5d" % (str(f)[:26], str(gg)[:22], str(d)[:20], n))
    keep = [c for c in COLS if c in sub.columns]
    out = "/tmp/modes_%s_%s.tsv" % (tag, len(idx))
    sub[keep].to_csv(out, sep="\t", index=False)
    print("   明细 -> %s" % out)
    print(sub[keep].head(6).to_string(max_colwidth=24))


def main():
    m = load_mod(CUR, "rhp")
    for tag, tree in TREES.items():
        print("\n================ %s ================" % tag, flush=True)
        df = load_frame(tree)
        res = {}
        for mode, ff, w3 in [("V0", True, False), ("V1", True, True), ("V2", False, False)]:
            m.FAMILY_FIRST_VETO = ff
            m.WHITELIST_OVERRIDES_FAMILY_VETO = w3
            res[mode] = [r[0] for r in df.apply(m.decision_tree_cascade, axis=1)]
            pd.DataFrame({"contig_id": df["contig_id"], "host": res[mode]}).to_csv(
                "/tmp/labels_%s_%s.tsv" % (tag, mode), sep="\t", index=False)
            print("   %s Plant=%d" % (mode, res[mode].count("Plant")), flush=True)
        m.FAMILY_FIRST_VETO = True
        m.WHITELIST_OVERRIDES_FAMILY_VETO = False

        idx_a = [i for i in range(len(df)) if res["V0"][i] != res["V1"][i]]
        print(" [V1-V0] 白名单越权(仅植物库认证属) 放行 %d 行" % len(idx_a))
        report(df, res["V1"], idx_a, "v1minusv0_%s" % tag, m)
        idx_b = [i for i in range(len(df)) if res["V1"][i] != res["V2"][i]]
        print(" [V2-V1] 科名单让位于属/种级判定 额外放行 %d 行" % len(idx_b))
        report(df, res["V2"], idx_b, "v2minusv1_%s" % tag, m)


with open("/tmp/three_modes2.log", "w", encoding="utf-8") as fh:
    old = sys.stdout
    sys.stdout = fh
    try:
        main()
    finally:
        sys.stdout = old
print(open("/tmp/three_modes2.log", encoding="utf-8").read())
