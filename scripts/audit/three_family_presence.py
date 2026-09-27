#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""核查 Mimiviridae / Marseilleviridae / Pithoviridae 三个科在
   ① 旧下游植物表  ② 新 05 表(v62d 闸门版)  里的行数与属分布。

只读；输出键值对，便于留档。
"""
import csv
import os
import sys

TARGETS = ["Mimiviridae", "Marseilleviridae", "Pithoviridae"]

FILES = [
    ("goji_植物表_旧", "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/"
                   "RNA-Lycium_barbarum_out/10_Reports/All_plant.viruses_info.tsv"),
    ("goji_05表_新", "/tmp/v62d_goji/final_integrated_classification.tsv"),
    ("onekp_植物表_旧", "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus/"
                    "09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv"),
    ("onekp_05表_新", "/tmp/v62d_onekp/final_integrated_classification.tsv"),
]


def scan(tag, path):
    if not os.path.isfile(path):
        print("[%s] 文件不存在: %s" % (tag, path))
        return
    hit = {}
    total = 0
    with open(path, encoding="utf-8", errors="replace", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        fam_col = "Family" if "Family" in (r.fieldnames or []) else None
        gen_col = "Genus" if "Genus" in (r.fieldnames or []) else None
        if fam_col is None:
            print("[%s] 无 Family 列; 列名: %s" % (tag, (r.fieldnames or [])[:12]))
            return
        for row in r:
            total += 1
            fam = (row.get(fam_col) or "").strip()
            if fam in TARGETS:
                gen = (row.get(gen_col) or "").strip()
                hit.setdefault((fam, gen), 0)
                hit[(fam, gen)] += 1
    print("[%s] %s" % (tag, path))
    print("    总行数 %d; 三科命中 %d 行" % (total, sum(hit.values())))
    for (fam, gen), n in sorted(hit.items(), key=lambda kv: (-kv[1], kv[0])):
        print("      %-18s Genus=%-24s %d" % (fam, gen or "(空)", n))
    if not hit:
        print("      三科零记录")


for tag, path in FILES:
    scan(tag, path)
