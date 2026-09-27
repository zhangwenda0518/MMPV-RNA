#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""白名单语义 W3 实测: "属级/种级 C9=Plant 且属在植物白名单" 时, 是否允许越过科级否决。

现状: is_blacklisted 在 FAMILY_FIRST_VETO=True 下先查科名单, 科命中即否决 (不管判定层)。
      后果: Family=Orthoherpesviridae + Genus=Macluravirus + det=Species 这类嵌合行,
      C9 依据植物属给出 Plant, 却被科级否决。
W3:   科命中时, 若判定层为属/种级 且 该属在植物白名单 (Plant.tsv 认证) → 放行 Plant。

实测两棵树: 现状 Plant 数, W3 后 Plant 数, 净增行数, 以及被救回行的证据明细 (抽样)。
注意: is_blacklisted 是全局函数, 这里用运行时替换实现 W3 (不落盘), 便于对拍。
"""
import importlib.util
import os

import pandas as pd

P = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
FAMS, GENS = set(m.NON_PLANT_FAMILIES_FALLBACK), set(m.NON_PLANT_GENERA)
WL_G = set(m.PLANT_GENERA_WHITELIST)

TREES = {
    "goji": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out",
    "onekp": "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus",
}


def load(tree):
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
        if "Unnamed" in idc or "y|virus order" in rvh.columns:
            rvh = rvh.rename(columns={idc: "contig_id"})
        df = df.merge(rvh[[c for c in ["contig_id", "pred|L1", "pred|L2", "evidence"]
                           if c in rvh.columns]], on="contig_id", how="left")
    pbp = os.path.join(out, "phabox2_output/final_prediction/cherry_prediction.tsv")
    if os.path.isfile(pbp):
        pb = pd.read_csv(pbp, sep="\t").rename(columns={"Accession": "contig_id"})
        df = df.merge(pb[[c for c in ["contig_id", "Host", "Host_NCBI_lineage"]
                          if c in pb.columns]], on="contig_id", how="left")
    return df


def label_all(df):
    return [r[0] for r in df.apply(m.decision_tree_cascade, axis=1)]


orig_black = m.is_blacklisted


def w3(family, genus, det_level):
    fam = str(family).strip() if family is not None else ""
    gen = str(genus).strip() if genus is not None else ""
    if fam in FAMS and m.is_trusted_level(det_level) and gen in WL_G:
        return False
    return orig_black(family, genus, det_level)


for tag, tree in TREES.items():
    df = load(tree)
    h0 = label_all(df)
    m.is_blacklisted = w3
    h1 = label_all(df)
    m.is_blacklisted = orig_black
    grp = df.groupby("Family")
    print("\n" + "=" * 90)
    print("[%s] 行 %d | 现状 Plant=%d | W3 Plant=%d | 净增 %d"
          % (tag, len(df), h0.count("Plant"), h1.count("Plant"), h1.count("Plant") - h0.count("Plant")))
    rows = []
    for i in range(len(df)):
        if h0[i] != h1[i]:
            r = df.iloc[i]
            rows.append({"contig_id": r.get("contig_id"), "Family": r.get("Family"),
                         "Genus": r.get("Genus"), "Det": r.get("Determination_Level"),
                         "C9": r.get("Host_ICTV"), "RVH": r.get("pred|L1"),
                         "PB2": r.get("Host"), "old": h0[i], "W3": h1[i]})
    d = pd.DataFrame(rows)
    if len(d):
        print("  翻转 %d 行, 按科统计:" % len(d))
        print(d.groupby(["Family", "Det"]).size().sort_values(ascending=False).head(12).to_string())
        print("  抽样前 8 行证据:")
        for r in d.head(8).itertuples(index=False):
            print("     %s F=%-22s G=%-18s det=%-20s C9=%s RVH=%s PB2=%s %s->%s"
                  % (str(r.contig_id)[:24], r.Family, r.Genus, str(r.Det), r.C9, r.RVH, r.PB2, r.old, r.W3))
        d.to_csv("/tmp/w3_rescued_%s.tsv" % tag, sep="\t", index=False)
        print("  明细写出 /tmp/w3_rescued_%s.tsv" % tag)
