#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""核实：产物里是否有占位值漏网 + 各工具原始表里 environmental 类取值分布。"""
import csv, os, re
from collections import defaultdict

BASE = ("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/"
        "RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated")
SHIP = os.path.join(BASE, "final_integrated_classification.tsv")

PAT = re.compile(r"environmental|unclassified|uncultured|unknown|"
                 r"\bsamples?\b|no rank|unplaced|^\s*$", re.I)

cnt = defaultdict(lambda: defaultdict(int))
n = 0
with open(SHIP, newline="", encoding="utf-8", errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        n += 1
        for k in ("Family", "Genus", "Species"):
            v = (r.get(k) or "").strip()
            if v and PAT.search(v):
                cnt[k][v] += 1
print("产物行数 %d" % n)
print("=== 产物里带占位/可疑字样的取值 ===")
for k in ("Family", "Genus", "Species"):
    if not cnt[k]:
        print("  [%s] 无" % k)
        continue
    print("  [%s]" % k)
    for v, c in sorted(cnt[k].items(), key=lambda kv: -kv[1])[:25]:
        print("     %6d  %s" % (c, v))

print("=== 各工具 standardized 表 Species 里含 environmental 的条数 ===")
for t in ["ACVirus", "VITAP", "mmseqs", "metabuli", "CAT", "genomad", "diamond_lca"]:
    p = os.path.join(BASE, "standardized_%s.tsv" % t)
    if not os.path.exists(p):
        print("  %-12s 缺文件" % t)
        continue
    hit = tot = 0
    samples = defaultdict(int)
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            v = (r.get("Species") or "").strip()
            if v:
                tot += 1
                if re.search(r"environmental", v, re.I):
                    hit += 1
                    samples[v] += 1
    print("  %-12s Species 有值 %6d, 含 environmental %4d" % (t, tot, hit))
    if samples:
        for v, c in sorted(samples.items(), key=lambda kv: -kv[1])[:5]:
            print("        %6d  %s" % (c, v))
