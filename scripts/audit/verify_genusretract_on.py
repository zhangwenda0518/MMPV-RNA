#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E1+E2 补丁上线后的验收:
 (1) 模块内 NON_PLANT_GENERA 已无 Miraophiovirus / Oleurovirus, 撒回表含 4 名;
 (2) normalize_c9 单一定义且语义为 strip/NA 归一版;
 (3) 两棵树 Plant 计数 = goji 1042 / onekp 19990;
 (4) 与 W3 开启时(V1)基线逐行对拍: goji 0 差异, onekp 恰 23 行且全部属 ∈ 撒回集。
只读。用法: python3 /tmp/verify_genusretract_on.py
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
EXPECT = {"goji": 1042, "onekp": 19990}
RETRACTED = {"Miraophiovirus", "Oleurovirus"}


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


m = load_mod(CUR, "rhp_after")
ok = True

print("[1] 名单状态")
for g in RETRACTED:
    a = g in m.NON_PLANT_GENERA
    b = g in m.NON_PLANT_GENERA_UNDER_REVIEW
    print("    %-16s in NON_PLANT_GENERA=%s (应为 False), in UNDER_REVIEW=%s (应为 True)" % (g, a, b))
    ok &= (not a) and b
print("    UNDER_REVIEW =", m.NON_PLANT_GENERA_UNDER_REVIEW)
print("    开关: FAMILY_FIRST_VETO=%s WHITELIST_OVERRIDES_FAMILY_VETO=%s"
      % (m.FAMILY_FIRST_VETO, m.WHITELIST_OVERRIDES_FAMILY_VETO))
print("    白名单∩黑名单 = %s (应为空)"
      % (set(m.PLANT_GENERA_WHITELIST) & set(m.NON_PLANT_GENERA)))

print()
print("[2] normalize_c9 语义 (合并后单一定义)")
src = open(CUR, encoding="utf-8").read()
n_def = src.count("\ndef normalize_c9")
print("    def normalize_c9 出现次数(含缩进外): %d" % n_def)
for probe in ["NA", "nan", " Plant", "", None, "Unknown", "Oomycetes"]:
    print("      normalize_c9(%-10r) -> %r" % (probe, m.normalize_c9(probe)))

print()
print("[3][4] 两棵树判定与 V1 基线对拍")
for tag, tree in TREES.items():
    df = load_frame(tree)
    now = pd.DataFrame({"contig_id": df["contig_id"],
                        "host_now": [r[0] for r in df.apply(m.decision_tree_cascade, axis=1)]})
    n = int(now["host_now"].eq("Plant").sum())
    p = "/tmp/labels_%s_V1.tsv" % tag
    old = pd.read_csv(p, sep="\t", dtype=str)
    j = now.merge(old, on="contig_id")
    d = j[j["host_now"] != j["host"]]
    print("    ==== %s: Plant=%d (期望 %d) %s" % (tag, n, EXPECT[tag], "PASS" if n == EXPECT[tag] else "FAIL"))
    ok &= (n == EXPECT[tag])
    sub = df[df["contig_id"].isin(set(d["contig_id"]))].drop_duplicates("contig_id")
    print("         与 V1 差异 %d 行 (期望 goji 0 / onekp 23)" % len(d))
    if len(sub):
        print(sub.groupby(["Family", "Genus", "Determination_Level", "Host_ICTV"]).size().to_string())
        gens = set(sub["Genus"].astype(str).str.strip())
        same = gens <= RETRACTED
        print("         差异行属集合 = %s -> %s" % (gens, "PASS 全在撒回集" if same else "FAIL"))
        ok &= same
print()
print("总判定:", "PASS" if ok else "FAIL")
