#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""只读：VMR Host source 列的 (S) 语义实测。
思路：若 (S) = 证据来自序列/环境来源而非确认宿主，则公认植物病毒属应几乎全是 "plants"（无 S），
      而 (S) 会集中在环境/未知类别上。"""
import csv
from collections import Counter

VMR = "/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv"

rows = []
with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t")
    next(rd)
    for r in rd:
        if len(r) >= 27:
            rows.append(r)

print("VMR 总行数 =", len(rows))
print("\n=== Host source 全量取值分布 (top 40) ===")
c = Counter(r[26].strip().strip('"') for r in rows)
for k, v in c.most_common(40):
    print("  %6d  %s" % (v, k))
print("  ... 取值种类总数 =", len(c))

print("\n=== 带 (S) 的取值 vs 不带的取值 ===")
withs = sum(v for k, v in c.items() if "(S)" in k)
nos = sum(v for k, v in c.items() if "(S)" not in k)
print("  含 (S) 行数 = %d ; 不含 = %d" % (withs, nos))

print("\n=== 公认植物病毒属的 Host source 分布（对照） ===")
by_genus = {}
for r in rows:
    by_genus.setdefault(r[15].strip().strip('"'), []).append(r)
for g in ["Potyvirus", "Tobamovirus", "Cucumovirus", "Tombusvirus", "Begomovirus",
          "Caulimovirus", "Sobemovirus", "Rimosavirus", "Miraophiovirus"]:
    rs = by_genus.get(g, [])
    cc = Counter(x[26].strip().strip('"') for x in rs)
    print("  %-16s 行数=%3d  %s" % (g, len(rs), dict(cc.most_common(6))))

print("\n=== 同时出现 plants 与 plants (S) 的属（(S) 是证据等级而非宿主类别的迹象） ===")
both = []
for g, rs in by_genus.items():
    vals = set(x[26].strip().strip('"') for x in rs)
    if any(v.startswith("plants") and "(S)" not in v for v in vals) and any("plants (S)" == v for v in vals):
        both.append((g, len(rs)))
print("  属数 =", len(both), "| 前 15:", both[:15])
