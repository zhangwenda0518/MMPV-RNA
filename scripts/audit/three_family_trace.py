#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""逐 contig 追踪: 旧植物表里 Family=Mimiviridae/Marseilleviridae/Pithoviridae 的行,
在 v62d 闸门版 05 表里变成了什么 (Family/Genus/是否还在)。

只读; 输出 旧(Family,Genus) -> 新(Family,Genus) 的计数矩阵。
"""
import csv
import os
from collections import Counter

TARGETS = {"Mimiviridae", "Marseilleviridae", "Pithoviridae"}

PAIRS = [
    ("goji",
     "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/"
     "RNA-Lycium_barbarum_out/10_Reports/All_plant.viruses_info.tsv",
     "/tmp/v62d_goji/final_integrated_classification.tsv"),
    ("onekp",
     "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus/"
     "09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv",
     "/tmp/v62d_onekp/final_integrated_classification.tsv"),
]


def load_new(path):
    out = {}
    with open(path, encoding="utf-8", errors="replace", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        key = "_id" if "_id" in (r.fieldnames or []) else r.fieldnames[0]
        for row in r:
            out[(row.get(key) or "").strip()] = (
                (row.get("Family") or "").strip(), (row.get("Genus") or "").strip())
    return out


for tag, old_path, new_path in PAIRS:
    print("=== %s ===" % tag)
    if not (os.path.isfile(old_path) and os.path.isfile(new_path)):
        print("  缺文件, 跳过"); continue
    newmap = load_new(new_path)
    print("  新表载入 %d 个 contig" % len(newmap))
    pairs = Counter()
    missing = 0
    with open(old_path, encoding="utf-8", errors="replace", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        first = r.fieldnames[0]
        for row in r:
            fam = (row.get("Family") or "").strip()
            if fam not in TARGETS:
                continue
            cid = (row.get(first) or "").strip()
            gen_old = (row.get("Genus") or "").strip() or "NA"
            if cid in newmap:
                nf, ng = newmap[cid]
                pairs[((fam, gen_old), (nf or "NA", ng or "NA"))] += 1
            else:
                missing += 1
    print("  旧表三科行: %d 行; 在新表找不到 contig: %d 行" % (sum(pairs.values()) + missing, missing))
    for (o, n), cnt in sorted(pairs.items(), key=lambda kv: (-kv[1], kv[0])):
        print("    %-18s/%-24s -> %-18s/%-24s  %d" % (o[0], o[1], n[0], n[1], cnt))
