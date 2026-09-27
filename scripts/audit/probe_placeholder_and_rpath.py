#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""核实占位值漏网 + 定位 R 脚本真实路径。"""
import csv, os, re, subprocess
from collections import defaultdict

BASE = ("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/"
        "RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated")
SHIP = os.path.join(BASE, "final_integrated_classification.tsv")

PAT = re.compile(r"environmental|unclassified|uncultured|unknown|^\s*$|"
                 r"\bsamples?\b|no rank|unplaced", re.I)
TAX = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

# 1) 产物里的占位值分布
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
print("\n=== 产物里带占位/可疑字样的取值 ===")
for k in ("Family", "Genus", "Species"):
    if not cnt[k]:
        print("  %-8s 无" % k)
        continue
    print("  [%s]" % k)
    for v, c in sorted(cnt[k].items(), key=lambda kv: -kv[1])[:25]:
        print("     %6d  %s" % (c, v))

# 2) 复刻侧：P1 结果里出现 environmental samples 的来源
print("\n=== 各工具 standardized 表里 Species 含 environmental samples 的条数 ===")
for t in ["ACVirus", "VITAP", "mmseqs", "metabuli", "CAT", "genomad", "diamond_lca"]:
    p = os.path.join(BASE, "standardized_%s.tsv" % t)
    if not os.path.exists(p):
        print("  %-12s 缺文件" % t)
        continue
    hit = tot = 0
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            v = (r.get("Species") or "").strip()
            if v:
                tot += 1
                if re.search(r"environmental", v, re.I):
                    hit += 1
    print("  %-12s Species 有值 %6d, 含 environmental %4d" % (t, tot, hit))

# 3) 定位 R 脚本
print("\n=== 定位 virus_classifier_analysis.R ===")
for root in ("/home/zhangwenda", "/home/zhangwenda/MMPV-RNA",
             "/home/zhangwenda/MMPV-paper"):
    if not os.path.isdir(root):
        continue
    for dp, dn, fn in os.walk(root):
        dn[:] = [d for d in dn if d not in (".git", "node_modules", "__pycache__")]
        for x in fn:
            if x == "virus_classifier_analysis.R":
                fp = os.path.join(dp, x)
                print("  %s  (mtime %s, %d bytes)"
                      % (fp, __import__("time").strftime(
                          "%Y-%m-%d %H:%M", __import__("time").localtime(os.path.getmtime(fp))),
                         os.path.getsize(fp)))
