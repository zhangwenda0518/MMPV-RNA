#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""开启 W3 之后的验证: 现装模块的判定结果必须 (a) 与开跑前的 V1 预测完全一致;
(b) 与 V0 的差集恰好是先前测到的那些行, 且每一行属都在植物白名单里。"""
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


def load_mod(path, name):
    loader = importlib.machinery.SourceFileLoader(name, path)
    spec = importlib.util.spec_from_loader(name, loader)
    m = importlib.util.module_from_spec(spec)
    loader.exec_module(m)
    return m


def load_frame(tree):
    """必须与 compare_veto_modes2.py 的 load_frame 完全一致: cascade 会读 RVH/PB2 的证据列。"""
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


m = load_mod(CUR, "rhp_installed")
print("现装 switch: FAMILY_FIRST_VETO=%s  WHITELIST_OVERRIDES_FAMILY_VETO=%s"
      % (m.FAMILY_FIRST_VETO, m.WHITELIST_OVERRIDES_FAMILY_VETO))
assert m.WHITELIST_OVERRIDES_FAMILY_VETO is True, "开关没开到位"

for tag, tree in TREES.items():
    df = load_frame(tree)
    now = pd.DataFrame({"contig_id": df["contig_id"],
                        "host": [r[0] for r in df.apply(m.decision_tree_cascade, axis=1)]})
    print("\n==== %s : 现装 Plant=%d (行数 %d)" % (tag, now["host"].eq("Plant").sum(), len(now)))
    for ref in ("V0", "V1"):
        p = "/tmp/labels_%s_%s.tsv" % (tag, ref)
        if not os.path.isfile(p):
            print("  %s 缺失, 跳过" % p)
            continue
        old = pd.read_csv(p, sep="\t", dtype=str)
        j = now.merge(old, on="contig_id", suffixes=("_now", "_ref"))
        d = j[j["host_now"] != j["host_ref"]]
        verdict = "PASS" if (ref == "V1" and len(d) == 0) else ("记录" if ref == "V0" else "FAIL")
        print("  与 %s 比对: contig 匹配 %d, 判定不同 %d  -> %s" % (ref, len(j), len(d), verdict))
        if ref == "V0" and len(d):
            sub = df[df["contig_id"].isin(set(d["contig_id"]))].copy()
            sub = sub.drop_duplicates("contig_id")
            inwl = sub["Genus"].astype(str).str.strip().isin(m.PLANT_GENERA_WHITELIST)
            print("    差集行属命中植物白名单 %d / 未命中 %d -> %s"
                  % (inwl.sum(), (~inwl).sum(), "PASS" if inwl.all() else "FAIL"))
            print("    差集行 (Family, Genus, Det) top10:")
            print(sub.groupby(["Family", "Genus", "Determination_Level"]).size()
                  .sort_values(ascending=False).head(10).to_string())
