#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""诊断: 科级否决 (FAMILY_FIRST_VETO) 实际拦下了哪些 C9=Plant 的行, 其中有多少属是植物属。

回答两个问题:
  Q1 科级否决在 goji / onekp 上共拦下多少行, 拦下的是什么 (科/属/判定层)。
  Q2 其中有多少行的属是植物属 (白名单命中) —— 即 "科级优先" 口径的代价面, 供是否放行的决策。
"""
import importlib.util
import os
from collections import Counter

import pandas as pd

P = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
FAMS, GENS = set(m.NON_PLANT_FAMILIES_FALLBACK), set(m.NON_PLANT_GENERA)
WL_F = set(m.PLANT_FAMILIES_WHITELIST)
WL_G = set(m.PLANT_GENERA_WHITELIST)
print("导入成功: 科黑名单 %d / 属黑名单 %d / 科白名单 %d / 属白名单 %d"
      % (len(FAMS), len(GENS), len(WL_F), len(WL_G)))

TREES = {
    "goji": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out",
    "onekp": "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus",
}

for tag, tree in TREES.items():
    tax = os.path.join(tree, "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")
    c9p = os.path.join(tree, "06_HostPrediction/C9_ICTV_result/classification_result.tsv")
    if not (os.path.isfile(tax) and os.path.isfile(c9p)):
        print("\n[%s] C9 产物缺失, 跳过" % tag)
        continue
    df = pd.read_csv(tax, sep="\t")
    c9 = pd.read_csv(c9p, sep="\t")[["contig_id", "Predicted_Host", "Determination_Level"]]
    df = df.merge(c9.rename(columns={"Predicted_Host": "Host_ICTV"}), on="contig_id", how="left")
    plant = df[df["Host_ICTV"] == "Plant"]
    fam_veto = plant[plant["Family"].isin(FAMS)]
    gen_veto = plant[(~plant["Family"].isin(FAMS)) &
                     (~plant["Determination_Level"].map(m.is_trusted_level)) &
                     (plant["Genus"].astype(str).isin(GENS))]
    print("\n%s: 总 %d | C9=Plant %d | 科级否决 %d | 属级否决 %d"
          % (tag, len(df), len(plant), len(fam_veto), len(gen_veto)))
    print("  科级否决明细 (科 ×属×判定层 top15):")
    c = Counter()
    for r in fam_veto.itertuples(index=False):
        c[(r.Family, getattr(r, "Genus"), r.Determination_Level, getattr(r, "Genus") in WL_G)] += 1
    for (f, g, d, wg), n in c.most_common(15):
        print("     %-24s %-22s det=%-18s 属在白名单=%s ×%d" % (f, g, d, wg, n))
    hit = [k for k in c if k[3]]
    print("  科级否决行中属命中植物白名单的组数: %d (行数 %d)"
          % (len(hit), sum(c[k] for k in hit)))
    if gen_veto.shape[0]:
        print("  属级否决明细:")
        for (g, f, d), n in Counter(
                (r.Genus, r.Family, r.Determination_Level) for r in gen_veto.itertuples(index=False)).most_common(15):
            print("     %-22s %-24s det=%-18s ×%d" % (g, f, d, n))
