#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""C9 输入域实测: Predicted_Host / Determination_Level 的取值集合,
以及 normalize_c9 两个版本 (497 行带 strip/NA 归一 vs 655 行现行) 在两棵树上的判定差异。
只读。用法: python3 /tmp/probe_c9_domain.py
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


def load_mod(path, name):
    loader = importlib.machinery.SourceFileLoader(name, path)
    spec = importlib.util.spec_from_loader(name, loader)
    m = importlib.util.module_from_spec(spec)
    loader.exec_module(m)
    return m


def norm_v1(h):  # 497 行版 (带 strip / NA 归一)
    h = str(h).strip() if not pd.isna(h) else 'nan'
    if h in ('', 'Unknown', 'None', 'NA', 'nan'):
        return 'Unknown'
    if h in ('Insecta', 'Arachnida', 'Aves', 'Human', 'Animal_other'):
        return 'Animal'
    if h == 'Oomycetes':
        return 'Protist'
    return h


print("=" * 88)
print("C9 输入域")
print("=" * 88)
for tag, tree in TREES.items():
    p = os.path.join(tree, "06_HostPrediction/C9_ICTV_result/classification_result.tsv")
    c9 = pd.read_csv(p, sep="\t")
    print()
    print("  ==== %s (%d 行)" % (tag, len(c9)))
    print("   Predicted_Host 取值:")
    for k, v in c9["Predicted_Host"].value_counts(dropna=False).items():
        print("      %-24r %d" % (k, v))
    print("   Determination_Level 取值:")
    for k, v in c9["Determination_Level"].value_counts(dropna=False).items():
        print("      %-40r %d" % (k, v))
    raw = c9["Predicted_Host"]
    odd = raw[raw.astype(str).str.strip().isin(['NA', 'nan', '', 'None']) | raw.isna()]
    print("   对 normalize 敏感的原值行数: %d" % len(odd))
    lvl = c9["Determination_Level"].astype(str)
    print("   判定级以 species* 开头的行数: %d" % int(lvl.str.lower().str.startswith("species*").sum()))

print()
print("=" * 88)
print("normalize_c9 两版本对最终判定的影响 (整树对拍)")
print("=" * 88)
m = load_mod(CUR, "rhp_probe")
m_alt = load_mod(CUR, "rhp_probe_alt")
m_alt.normalize_c9 = norm_v1  # 内存替换, 不动文件


def load_frame(tree, mod):
    """用给定模块的 normalize_c9 语义重算, 其余接线与 compare_veto_modes2 一致。"""
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


for tag, tree in TREES.items():
    a = load_frame(tree, m)
    b = a.copy()
    ha = [r[0] for r in a.apply(m.decision_tree_cascade, axis=1)]
    hb = [r[0] for r in b.apply(m_alt.decision_tree_cascade, axis=1)]
    diff = sum(1 for x, y in zip(ha, hb) if x != y)
    print("  %-6s Plant 现行 %d  vs  strip/NA 归一版 %d, 判定差异 %d  -> %s"
          % (tag, ha.count("Plant"), hb.count("Plant"), diff,
             "PASS 无差异" if diff == 0 else "FAIL 有差异"))
print()
print("DONE")
