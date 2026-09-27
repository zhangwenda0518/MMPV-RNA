#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""465 行 Biavirus 裁决第四层：ACVirus 自身置信度分档。

ACVirus 置信度定义（cl i.py:898 assign_confidence_scores_dynamic）：
  对每个分类等级，若 contig 的 Max_Coverage(%) >= 该 taxon 的 Min_Coverage(%) 阈值
  → High-confidence，否则 Low-confidence。
Min_Coverage 阈值（cli.py:543 calculate_min_coverage_final_robust）：
  取「库内该 taxon 成员之间」参考序列互比的最小覆盖度。
  Biavirus 库里只有 HG999358 一个基因组 → 自比 ≈ 99.46% → 属级达标线近乎不可达。
"""
import csv
import os
from collections import Counter, defaultdict

BASE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
VOT = os.path.join(BASE, "05_Taxonomy", "Votus.integrated")
BAD = os.path.join(VOT, "calibration_20260914", "biavirus_arbitration")
AC = os.path.join(BASE, "05_Taxonomy", "Votus.classed", "ACVirus_results", "Votus.acvirus")


def p(*a):
    print(*a, flush=True)


ids = [l.strip() for l in open(os.path.join(BAD, "selected_ids.txt"), encoding="utf-8") if l.strip()]
idset = set(ids)
p("目标 contig 数: %d" % len(ids))

p("\n=== A. medium_result.csv 里目标行的完整分类字段（前 3 行）===")
with open(os.path.join(AC, "medium_result.csv"), encoding="utf-8", errors="replace") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    cols = rd.fieldnames
    p("列: %s" % cols)
    shown = 0
    for d in rd:
        if (d.get("Nucleotide") or "").strip() in idset:
            p("   " + " | ".join("%s=%s" % (c, (d.get(c) or "").strip()) for c in cols))
            shown += 1
            if shown >= 3:
                break

p("\n=== B. final_result_with_confidence.tsv 目标行置信档 ===")
fp = os.path.join(AC, "final_result_with_confidence.tsv")
rd = csv.DictReader(open(fp, newline="", encoding="utf-8", errors="replace"), delimiter="\t")
cols = rd.fieldnames
p("列数 %d: %s" % (len(cols), cols))
conf_cols = [c for c in cols if c.endswith("_Confidence")]
lvls = [c.replace("_Confidence", "") for c in conf_cols]
per_level = defaultdict(Counter)
n_in = 0
allc = defaultdict(Counter)
sample = None
for d in rd:
    key = (d.get("Nucleotide") or "").strip()
    for c in conf_cols:
        allc[c][(d.get(c) or "(空)").strip()] += 1
    if key not in idset:
        continue
    n_in += 1
    if sample is None:
        sample = d
    for lv, c in zip(lvls, conf_cols):
        per_level[lv][(d.get(c) or "(空)").strip()] += 1

p("目标行命中 %d / %d" % (n_in, len(ids)))
p("\n-- 目标行各级置信档（含该级分类单元出现次数）--")
for lv in lvls:
    p("   %-12s %s" % (lv, dict(per_level[lv])))
p("\n-- 全表（所有 contig）各级置信档对照 --")
for c in conf_cols:
    p("   %-24s %s" % (c, dict(allc[c])))

p("\n=== C. 目标行的分类单元取值分布（来自 with_confidence 表）===")
tax = defaultdict(Counter)
for d in csv.DictReader(open(fp, newline="", encoding="utf-8", errors="replace"), delimiter="\t"):
    if (d.get("Nucleotide") or "").strip() in idset:
        for lv in lvls:
            tax[lv][(d.get(lv) or "(空)").strip()] += 1
for lv in lvls:
    p("   %-12s %s" % (lv, dict(tax[lv].most_common(4))))

p("\n=== D. 目标行 × 属级置信 交叉 ===")
cross = Counter()
for d in csv.DictReader(open(fp, newline="", encoding="utf-8", errors="replace"), delimiter="\t"):
    if (d.get("Nucleotide") or "").strip() in idset:
        g = (d.get("Genus") or "").strip()
        gc = (d.get("Genus_Confidence") or "").strip() or "(无)"
        cross[(g, gc)] += 1
for k, v in cross.most_common():
    p("   Genus=%-14s Confidence=%-16s %d" % (k[0], k[1], v))

p("\n=== E. 分类表里这一行的 Family 列 (对比 integrated 表) ===")
mn = {}
with open(os.path.join(BASE, "05_Taxonomy", "Votus.integrated", "final_integrated_classification.tsv"),
          encoding="utf-8", errors="replace") as fh:
    rd2 = csv.DictReader(fh, delimiter="\t")
    for d in rd2:
        cid = (d.get("contig_id") or "").strip().strip('"')
        if cid in idset:
            mn[(d.get("Family") or "").strip().strip('"')] = mn.get((d.get("Family") or "").strip().strip('"'), 0) + 1
p("   integrated 表 Family 分布: %s" % mn)
