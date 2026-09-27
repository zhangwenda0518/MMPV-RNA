#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成 run_host_prediction.py 的植物白名单代码块 (文本), 写到 /tmp/plant_whitelist_block.py。

白名单构成:
  A. Plant.tsv (位置口径 index5=Family / index6=Genus) 出现的全部科与属
  B. ICTV 官方 Virus Properties by Family 中 Host 含 plants, 但 Plant.tsv 零覆盖的科
     (这 9 科正是 2026-09-15 从非植物科名单里撒回的那些, 必须进白名单, 否则下次又被误加)

模板用 @@TOKEN@@ 替换 (不用 str.format, 避免与生成代码里的花括号冲突)。
"""
import os
from collections import Counter

PLANT = os.path.expanduser("~/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv")
OUT = "/tmp/plant_whitelist_block.py"

ICTV_PLANT_INCLUDING = {
    "Artoviridae": "plants, invertebrates",
    "Chrysoviridae": "fungi, plants, invertebrates",
    "Genomoviridae": "fungi, plants, invertebrates, vertebrates",
    "Kanorauviridae": "plants, invertebrates, vertebrates",
    "Mitoviridae": "fungi, plants",
    "Ourmiaviridae": "invertebrates, plants",
    "Pestiviridae": "vertebrates, invertebrates, plants",
    "Pseudoviridae": "protists, fungi, plants, invertebrates",
    "Spiciviridae": "plants, invertebrates",
    "Tomosaviridae": "plants",
    "Virgaviridae": "plants",
    "Solemoviridae": "plants",
}

fam, gen = Counter(), Counter()
with open(PLANT, encoding="utf-8", errors="replace") as fh:
    fh.readline()
    for line in fh:
        fs = line.rstrip("\n").split("\t")
        if len(fs) < 8:
            continue
        p = fs[3].split(";")
        if len(p) > 5 and p[5].strip():
            fam[p[5].strip()] += 1
        if len(p) > 6 and p[6].strip():
            gen[p[6].strip()] += 1

zd = sorted(set(ICTV_PLANT_INCLUDING) - set(fam))
families = sorted(set(fam) | set(ICTV_PLANT_INCLUDING))
genera = sorted(set(gen))


def wrap(names, per_line=5, indent="    "):
    return "\n".join(indent + ", ".join("'%s'" % n for n in names[i:i + per_line]) + ","
                     for i in range(0, len(names), per_line))


TPL = """# ── 2026-09-15 植物白名单 (防误否决, 只做一致性断言, 不参与否决决策) ──
# 构造 A: PlantVirusDB 权威植物库 Plant.tsv 按位置口径解析 (Virus_lineage 9 段,
#         index5=Family / index6=Genus, 与 3_host_probability.py 的 levels 一致):
#         @@NFAM@@ 科 / @@NGEN@@ 属 (Plant.tsv md5 见 HOST_REGISTRY_AUDIT_20260915.md)
# 构造 B: ICTV 官方 Virus Properties by Family Host 列含 plants, 但 Plant.tsv 零覆盖的 @@NZD@@ 科:
#         @@ZDLIST@@
#         这 @@NZD@@ 科即 2026-09-15 从 NON_PLANT_FAMILIES_FALLBACK 撒回的全部科名
#         (审计: scripts/audit/HOST_REGISTRY_AUDIT_20260915.md, 证据: ICTV_HOST_SWEEP_233.tsv)。
# 用途: (1) 机器化不变量 "白名单 交 黑名单 = 空", 在导入时 fail-loud;
#       (2) 供审计脚本/人复核时直接检索, 无需重跑 Plant.tsv 解析。
#       注意: 白名单【不覆盖科级否决】—— 科级否决是刻意的嵌合拦截
#       (例: Potyvirus 属的嵌合行被标 Mimiviridae 时按科否决),
#       若让白名单越过否决会反向制造假阳性, 故此处仅登记, 不改变否决行为。
PLANT_FAMILIES_WHITELIST = frozenset([
@@FAMS@@
])

PLANT_GENERA_WHITELIST = frozenset([
@@GENS@@
])

# ICTV 现行宿主列含 plants 的 12 科 (科名 -> ICTV Host 列原文), 审计取数用
ICTV_PLANT_INCLUDING_FAMILIES = {
@@ICTV@@
}

_pl_wl_bad_f = sorted(set(NON_PLANT_FAMILIES_FALLBACK) & PLANT_FAMILIES_WHITELIST)
_pl_wl_bad_g = sorted(set(NON_PLANT_GENERA) & PLANT_GENERA_WHITELIST)
if _pl_wl_bad_f or _pl_wl_bad_g:
    raise RuntimeError(
        '宿主名单自相矛盾: 植物白名单与黑名单相交 '
        '科=' + repr(_pl_wl_bad_f) + ' 属=' + repr(_pl_wl_bad_g) +
        ' (历史上出现过 Plant.tsv 零记录即判非植物科的误否决, 见 HOST_REGISTRY_AUDIT_20260915.md)')
"""

block = (TPL
         .replace("@@NFAM@@", str(len(families)))
         .replace("@@NGEN@@", str(len(genera)))
         .replace("@@NZD@@", str(len(zd)))
         .replace("@@ZDLIST@@", ", ".join(zd))
         .replace("@@FAMS@@", wrap(families))
         .replace("@@GENS@@", wrap(genera))
         .replace("@@ICTV@@", "\n".join("    '%s': '%s'," % (k, ICTV_PLANT_INCLUDING[k])
                                        for k in sorted(ICTV_PLANT_INCLUDING))))

with open(OUT, "w", encoding="utf-8") as fh:
    fh.write(block)
print("写入 %s: 科 %d (Plant.tsv %d + 零覆盖 %d) / 属 %d"
      % (OUT, len(families), len(fam), len(zd), len(genera)))
print("Plant.tsv 零覆盖的 ICTV 含植物科: %s" % zd)
