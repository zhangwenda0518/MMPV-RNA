#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""争议条目最终取证:
 A. VMR 里 host=plants 的具体种名 (Ouroboviridae / Discoviridae / Rimosavirus)
 B. 两棵树里 "VMR 无记录" 的 7 个属黑名单条目的行数与释放后的 Plant 增量 (逐条内存撤回)
只读。用法: python3 /tmp/vmr_contested_final.py
"""
import csv
import importlib.machinery
import importlib.util
import os
import re

import pandas as pd

PIPE = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline"
CUR = os.path.join(PIPE, "run_host_prediction.py")
VMR = "/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv"
TREES = {
    "goji": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out",
    "onekp": "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus",
}

rows = list(csv.reader(open(VMR, encoding="utf-8", errors="replace"), delimiter="\t"))
SPEC, FAM, GEN, VNAME, HOST = 17, 13, 15, 20, 26
print("=" * 92)
print("A. VMR host=plants 的科级争议条目, 具体种/病毒名")
print("=" * 92)
want_fam = {"Ouroboviridae", "Discoviridae"}
for r in rows[1:]:
    if len(r) <= HOST:
        continue
    if r[FAM].strip() in want_fam and "plant" in r[HOST].lower():
        print("  %-16s %-22s %-40s host=%s" % (r[FAM].strip(), r[GEN].strip(),
                                               r[VNAME].strip()[:40], r[HOST].strip()))
print()
print("  Rimosavirus (Tombusviridae) 全部记录:")
for r in rows[1:]:
    if len(r) > HOST and r[GEN].strip() == "Rimosavirus":
        print("     %-20s %-42s host=%s" % (r[FAM].strip(), r[VNAME].strip()[:42], r[HOST].strip()))
print()
print("  Sylvanvirus 在 VMR 的近似名 (含 sylv):")
for r in rows[1:]:
    if len(r) > HOST and re.search(r"sylv", r[GEN] + r[VNAME], re.I):
        print("     %-16s %-22s %-40s host=%s" % (r[FAM].strip(), r[GEN].strip(),
                                                  r[VNAME].strip()[:40], r[HOST].strip()))


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


m0 = load_mod(CUR, "rhp0")
NO_VMR = ["Sylvanvirus", "Klosneuvirus", "Hokovirus", "Indivirus", "Clandestinovirus",
          "Crucivirus", "Arenavirus"]
CONTENDED = NO_VMR + ["Miraophiovirus", "Oleurovirus", "Rimosavirus"]
print()
print("=" * 92)
print("B. 逐条撤回的召回效果 (内存, 不动文件)")
print("=" * 92)

for tag, tree in TREES.items():
    df = load_frame(tree)
    base = [r[0] for r in df.apply(m0.decision_tree_cascade, axis=1)]
    n0 = base.count("Plant")
    print()
    print("  ==== %s: 行 %d, Plant %d" % (tag, len(df), n0))
    gen_s = df["Genus"].astype(str).str.strip()
    for g in CONTENDED:
        sub = df[gen_s == g]
        if len(sub) == 0:
            print("     %-18s 树内 0 行" % g)
            continue
        host = [base[i] for i in sub.index]
        m2 = load_mod(CUR, "rhp_tmp_%s" % g)
        if g in m2.NON_PLANT_GENERA:
            m2.NON_PLANT_GENERA.remove(g)
        new = [r[0] for r in sub.apply(m2.decision_tree_cascade, axis=1)]
        delta = sum(1 for a, b in zip(host, new) if a != b)
        flip = sum(1 for a, b in zip(host, new) if a != "Plant" and b == "Plant")
        dist = {str(k): int(v) for k, v in pd.Series(host).value_counts().items()}
        print("     %-18s 树内 %4d 行, 剥离该条后翻转 %3d (其中 ->Plant %3d), 现判定分布 %s"
              % (g, len(sub), delta, flip, dist))
print()
print("DONE")
