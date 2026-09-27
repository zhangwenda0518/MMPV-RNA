#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""系统性反查: 登记为非植物的科里, 是否存在"其成员属在植物病毒库中出现"的科。

数据源:
  A. 现装 run_host_prediction.py 的 NON_PLANT_FAMILIES_FALLBACK (233 条)
  B. ~/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv (植物病毒库)
  C. ~/database/taxonomy/genus_family_ref.tsv (属->科参照)
判定: 属在 Plant.tsv 出现 → 该属属植物病毒 → 其所属科可疑 (可能是改名/别名, 例如
      Ophioviridae -> Aspiviridae 之类的 ICTV 改名会让"零记录"判据失效)。
"""
import importlib.util
import os
from collections import Counter, defaultdict

import pandas as pd

RHP = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", RHP)
rhp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rhp)
REG = sorted(set(rhp.NON_PLANT_FAMILIES_FALLBACK))
print("登记科数: %d" % len(REG))

PLANT = os.path.expanduser("~/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv")
REF = os.path.expanduser("~/database/taxonomy/genus_family_ref.tsv")
print("Plant.tsv 存在: %s ; 参照表存在: %s" % (os.path.isfile(PLANT), os.path.isfile(REF)))

pl = pd.read_csv(PLANT, sep="\t", low_memory=False)
print("Plant.tsv 列: %s (%d 行)" % (list(pl.columns)[:12], len(pl)))
lin = "Virus_lineage" if "Virus_lineage" in pl.columns else pl.columns[0]


def fam_of(lineage):
    """Plant.tsv 的 Virus_lineage 是无前缀的 Realm;Kingdom;...;Family;Genus;Species;Name,
    按 'viridae' 结尾定位科名。"""
    parts = [p.strip() for p in str(lineage).split(";")]
    for i, p in enumerate(parts):
        if p.endswith("viridae"):
            return p, (parts[i + 1] if i + 1 < len(parts) else "")
    return "", ""


_pairs = pl[lin].map(fam_of)
pl["_fam"] = [p[0] for p in _pairs]
pl["_gen"] = [p[1] for p in _pairs]
plant_fams = sorted({f for f in pl["_fam"] if f})
plant_gens = {g for g in pl["_gen"] if g}
print("植物库: %d 科 / %d 属" % (len(plant_fams), len(plant_gens)))
print("植物库科表: %s" % ", ".join(plant_fams))

ref = pd.read_csv(REF, sep="\t")
print("\n参照表列: %s" % list(ref.columns))
gcol = [c for c in ref.columns if "genus" in c.lower()][0]
fcol = [c for c in ref.columns if "family" in c.lower()][0]
fam2gen = defaultdict(set)
for g, f in zip(ref[gcol].astype(str), ref[fcol].astype(str)):
    fam2gen[f.strip()].add(g.strip())

hits = []
for f in REG:
    gens = fam2gen.get(f, set())
    inter = sorted(gens & plant_gens)
    if inter:
        hits.append((f, len(gens), inter))
print("\n【告警】登记为非植物科, 但其成员属出现在植物病毒库中 (%d 科):" % len(hits))
for f, n, inter in hits:
    print("   %-22s 参照属数 %4d | 与植物库交集 %s" % (f, n, inter))

ref_fams = {f for f in fam2gen}
print("\n登记科中在参照表里查不到属映射的: %d" % sum(1 for f in REG if f not in ref_fams))
# 反向: 植物库的科里, 有哪些被登记为非植物
both = sorted(set(plant_fams) & set(REG))
print("植物库科与登记名单直接重名: %s" % (both or "无"))
