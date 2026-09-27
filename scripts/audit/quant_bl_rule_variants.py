#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""比较 Plant 否决口径 (内存替换决策树, 不改源码):
  无否决 : 等价现网代码 (名单恒不生效) → Plant 上界
  V2     : 科在非植物科名单 → 直接否决 Plant
  V3     : 科在名单 → 比较科级与深层级投票支持数 (k/n), 科级不弱于深层才否决
名单已按 ICTV 宿主核定剔除 Ourmiaviridae / Mitoviridae / Kanorauviridae / Metaviridae。
"""
import importlib.util
import os
import re
from collections import Counter

import pandas as pd

RHP = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", RHP)
rhp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rhp)

TREE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
TAX = os.path.join(TREE, "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")
OUT = os.path.join(TREE, "06_HostPrediction")
C9 = os.path.join(OUT, "C9_ICTV_result/classification_result.tsv")
RVH = os.path.join(OUT, "RVH_result/result.csv")
PB2 = os.path.join(OUT, "phabox2_output/final_prediction/cherry_prediction.tsv")

df = pd.read_csv(TAX, sep="\t")
c9 = pd.read_csv(C9, sep="\t")[["contig_id", "Predicted_Host", "Determination_Level"]].rename(
    columns={"Predicted_Host": "Host_ICTV"})
df = df.merge(c9, on="contig_id", how="left")
rvh = pd.read_csv(RVH)
idc = rvh.columns[0]
if "Unnamed" in idc or "y|virus order" in rvh.columns:
    rvh = rvh.rename(columns={idc: "contig_id"})
df = df.merge(rvh[[c for c in ["contig_id", "pred|L1", "pred|L2", "evidence"] if c in rvh.columns]],
              on="contig_id", how="left")
pb2 = pd.read_csv(PB2, sep="\t").rename(columns={"Accession": "contig_id"})
df = df.merge(pb2[[c for c in ["contig_id", "Host", "Host_NCBI_lineage", "Host_GTDB_lineage"]
                   if c in pb2.columns]], on="contig_id", how="left")

FALLBACK = set(rhp.NON_PLANT_FAMILIES_FALLBACK)
GENERA = set(rhp.NON_PLANT_GENERA)
TRUSTED = set(rhp.TRUSTED_LEVELS)
REMOVE = {"Ourmiaviridae", "Mitoviridae", "Kanorauviridae", "Metaviridae"}
FALLBACK -= REMOVE
print("模拟名单 %d 科 (剔除 %s)" % (len(FALLBACK), sorted(REMOVE)))
_AGREE = re.compile(r"^\s*(\d+)\s*/\s*(\d+)")


def _norm(x):
    return "" if x is None or (isinstance(x, float) and pd.isna(x)) else str(x).strip()


def _bl_common(fam, gen, det):
    """属/科通用后半段: 判定层受信则放行, 否则查属名单。"""
    d = str(det).split("(")[0].strip().lower() if det is not None else ""
    if d in TRUSTED:
        return False
    g = _norm(gen)
    return bool(g and g not in ("NA", "nan") and g in GENERA)


def bl_none(row):
    return False


def bl_v2(row):
    f = _norm(row.get("Family"))
    if f and f in FALLBACK:
        return True
    return _bl_common(f, row.get("Genus"), row.get("Determination_Level"))


def _support(row, det):
    def k(col):
        m = _AGREE.match(_norm(row.get(col)))
        return int(m.group(1)) if m else 0
    d = str(det).split("(")[0].strip().lower() if det is not None else ""
    deep = k("Species_agree") if "species" in d else (
        k("Genus_agree") if "genus" in d else (k("Order_agree") if "order" in d else 0))
    return k("Family_agree"), deep


def bl_v3(row):
    f, det = _norm(row.get("Family")), row.get("Determination_Level")
    if f and f in FALLBACK:
        fam_k, deep_k = _support(row, det)
        return fam_k >= max(deep_k, 1)
    return _bl_common(f, row.get("Genus"), det)


def make_tree(bl_row):
    def tree(row):
        if str(row.get("Class", "")).lower() in rhp.PHAGE_CLASSES:
            return "Bacteria", "Rule_Class_Taxonomy"
        h = rhp.normalize_c9(row.get("Host_ICTV", "Unknown"))
        if h != "Unknown":
            if h == "Plant" and bl_row(row):
                h = "Unknown"
            else:
                return h, "ICTV_Preferred"
        r = rhp.parse_rvh(row)
        if r != "Unknown":
            return r, "RVH_Preferred"
        p = rhp.parse_pb2(row)
        if p != "Unknown":
            return p, "PB2_Preferred"
        return "Unknown", "Unassigned"
    return tree


def run(bl_row, tag):
    res = df.apply(make_tree(bl_row), axis=1)
    h = [r[0] for r in res]
    print("  [%s] Plant=%d | %s" % (tag, h.count("Plant"),
                                    dict(sorted(Counter(h).items(), key=lambda x: -x[1]))))
    return h


h0 = run(bl_none, "无否决 (现网等价)")
h2 = run(bl_v2, "V2 科级优先")
h3 = run(bl_v3, "V3 支持数比较")

d23 = [i for i in range(len(df)) if h2[i] != h3[i]]
print("\nV2 vs V3 差异 %d 行 (V3 保留 Plant 的 %d 行)"
      % (len(d23), sum(1 for i in d23 if h3[i] == "Plant")))
cols = ["contig_id", "Family", "Family_agree", "Genus", "Genus_agree", "Species",
        "Species_agree", "Determination_Level"]
for i in d23:
    r = df.iloc[i]
    print("   V2=%-8s V3=%-8s | %s" % (h2[i], h3[i],
                                       " ".join("%s=%s" % (c, _norm(r.get(c))) for c in cols)))
d02 = [i for i in range(len(df)) if h0[i] != h2[i]]
print("\nV2 相对现网剔除的 Plant 行 %d; 其中 V3 也剔除 %d"
      % (sum(1 for i in d02 if h0[i] == "Plant"),
         sum(1 for i in d02 if h0[i] == "Plant" and h3[i] != "Plant")))
