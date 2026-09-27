#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""1) 核对拼写差异 (Fervensiviridae vs Fervensviridae) 与 6 个非现行科名在 05 表中的实际出现;
2) 断言: 现装名单中不含任何 ICTV 宿主含植物的科 (对照 ICTV_PLANT_HOST_FLAG 口径的 12 科)。"""
import importlib.util
import os
from collections import Counter

import pandas as pd

P = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
REG = set(m.NON_PLANT_FAMILIES_FALLBACK)
GREG = set(m.NON_PLANT_GENERA)
print("现装: 科 %d / 属 %d" % (len(REG), len(GREG)))

TAXES = {
    "goji-旧 05 表": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/"
                 "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv",
    "goji-闸门版": "/tmp/v62d_goji/final_integrated_classification.tsv",
}
CHECK = ["Fervensiviridae", "Fervensviridae", "Microviridae", "Autographiviridae",
         "Pandoraviridae", "Cruciviridae", "Adomaviridae", "Reoviridae",
         "Kanorauviridae", "Tomosaviridae"]

for tag, path in TAXES.items():
    if not os.path.isfile(path):
        print("\n%s: 文件不存在 %s" % (tag, path))
        continue
    df = pd.read_csv(path, sep="\t", usecols=["Family", "Genus"], low_memory=False)
    fc = Counter(df["Family"].fillna("NA").astype(str))
    print("\n%s (共 %d 行, %d 个 Family 值)" % (tag, len(df), len(fc)))
    for c in CHECK:
        print("   %-20s 行数=%-6d in-list=%s" % (c, fc.get(c, 0), c in REG))
    misspelled = fc.get("Fervensiviridae", 0)
    official = fc.get("Fervensviridae", 0)
    print("   -> 拼写: 结果表用 Fervensiviridae %d 行 / 官方 Fervensviridae %d 行"
          % (misspelled, official))

PLANT_INCLUDED = {"Artoviridae", "Chrysoviridae", "Genomoviridae", "Kanorauviridae",
                  "Mitoviridae", "Ourmiaviridae", "Pestiviridae", "Pseudoviridae",
                  "Spiciviridae", "Tomosaviridae", "Virgaviridae", "Solemoviridae",
                  "Metaviridae"}
print("\n断言: 名单内残留的宿主含植物科 -> %s"
      % (sorted(REG & PLANT_INCLUDED) or "无 (通过)"))
print("非现行科名仍在名单内 (无植物风险, 仅标注): %s"
      % sorted(REG & {"Microviridae", "Autographiviridae", "Pandoraviridae",
                      "Cruciviridae", "Adomaviridae", "Reoviridae"}))
print("名单内查不到官方条目的科 (拼写疑点): %s" % sorted(REG & {"Fervensiviridae"}))
