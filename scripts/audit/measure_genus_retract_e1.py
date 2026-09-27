#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E1/E2/E3/S1 核实: 属黑名单里的真植物病毒属 (Miraophiovirus / Oleurovirus)
造成的召回损失, 以及 normalize_c9 双定义、TRUSTED_LEVELS 死元素、否决被 RVH/PB2 复活的计数。

只读: 只在内存里改模块对象, 不写任何管线文件。
用法: python3 /tmp/measure_genus_retract_e1.py
"""
import copy
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
TARGET_GENERA = ["Miraophiovirus", "Oleurovirus"]


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
    return df, out


m = load_mod(CUR, "rhp_inst")

print("=" * 78)
print("[E2] normalize_c9 双定义 / [E3] TRUSTED_LEVELS 死元素")
print("=" * 78)
src = open(CUR, encoding="utf-8").read().splitlines()
defs = [i + 1 for i, l in enumerate(src) if l.startswith("def normalize_c9")]
print("  def normalize_c9 出现行号:", defs)
live = m.normalize_c9
print("  运行时 normalize_c9 来自第 %s 行 (后定义者生效)"
      % (len(defs) and defs[-1]))
for probe in ["NA", "nan", " Plant", "Plant", "", None]:
    try:
        r = live(probe)
    except Exception as e:  # noqa: BLE001
        r = "EXC:%s" % e
    print("    normalize_c9(%-8r) -> %r" % (probe, r))
print("  TRUSTED_LEVELS =", m.TRUSTED_LEVELS)

print()
print("=" * 78)
print("[E1] 属黑名单误否决的真植物病毒属: 现装 vs 内存中撤回该属")
print("=" * 78)
for g in TARGET_GENERA:
    print("  %s in NON_PLANT_GENERA=%s  in PLANT_GENERA_WHITELIST=%s"
          % (g, g in m.NON_PLANT_GENERA, g in m.PLANT_GENERA_WHITELIST))

res = {}
for tag, tree in TREES.items():
    df, outdir = load_frame(tree)
    m_cur = m
    base = df.apply(m_cur.decision_tree_cascade, axis=1)
    base_host = [r[0] for r in base]
    base_method = [r[1] for r in base]

    # 内存中撤回 (不动文件)
    m2 = load_mod(CUR, "rhp_retract_test")
    removed = []
    for g in TARGET_GENERA:
        if g in m2.NON_PLANT_GENERA:
            m2.NON_PLANT_GENERA.remove(g)
            removed.append(g)
    new = df.apply(m2.decision_tree_cascade, axis=1)
    new_host = [r[0] for r in new]
    new_method = [r[1] for r in new]

    diff = [i for i in range(len(df)) if base_host[i] != new_host[i]]
    print()
    print("  ==== %s (%d 行), 内存撤回 %s" % (tag, len(df), removed))
    print("       Plant: 现装 %d -> 撤回后 %d (%+d)"
          % (base_host.count("Plant"), new_host.count("Plant"),
             new_host.count("Plant") - base_host.count("Plant")))
    print("       受影响行 %d" % len(diff))
    if diff:
        sub = df.iloc[diff].copy()
        sub["_old"] = [base_host[i] for i in diff]
        sub["_new"] = [new_host[i] for i in diff]
        sub["_old_m"] = [base_method[i] for i in diff]
        sub["_new_m"] = [new_method[i] for i in diff]
        print(sub.groupby(["Family", "Genus", "Determination_Level", "Host_ICTV",
                           "_old", "_new", "_old_m", "_new_m"]).size().to_string())

    # 目标行现状: 全部行 (不只 diff), 看是否有复活
    tgt = df[df["Genus"].astype(str).str.strip().isin(TARGET_GENERA)].copy()
    tgt["_host"] = [base_host[i] for i in tgt.index]
    tgt["_method"] = [base_method[i] for i in tgt.index]
    tgt["_bl"] = df.loc[tgt.index].apply(
        lambda r: m.is_blacklisted(r.get("Family"), r.get("Genus"), r.get("Determination_Level")),
        axis=1)
    print("       目标属行共 %d, 现装判定分布:" % len(tgt))
    print(tgt.groupby(["Family", "Genus", "Determination_Level", "Host_ICTV",
                       "_bl", "_host", "_method"]).size().to_string())

    # 树里已存的宿主产物 (用于与真实落盘值对照)
    for fn in sorted(os.listdir(outdir)):
        p = os.path.join(outdir, fn)
        if os.path.isfile(p) and fn.endswith((".tsv", ".csv")):
            print("       [产物] %s" % fn)

    # S1: 被否决但最终仍是 Plant 的行 (RVH/PB2 复活)
    vet = [i for i in range(len(df))
           if m.is_blacklisted(df.iloc[i].get("Family"), df.iloc[i].get("Genus"),
                               df.iloc[i].get("Determination_Level"))]
    if vet:
        sub = df.iloc[vet].copy()
        sub["_host"] = [base_host[i] for i in vet]
        sub["_method"] = [base_method[i] for i in vet]
        res[tag] = sub
        print("       [S1] 命中否决行 %d, 其中最终 Plant %d"
              % (len(vet), int((sub["_host"] == "Plant").sum())))
        print(sub[sub["_host"] == "Plant"].groupby(
            ["Family", "Genus", "_method"]).size().sort_values(ascending=False).head(12).to_string())

print()
print("DONE")
