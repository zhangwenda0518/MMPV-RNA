#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VMR Host source 列 (col26) 取值分布 + 含 plant 的原值, 用于确认 (S) 标记语义。只读。"""
import collections

P = "/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv"
raw = open(P, encoding="utf-8", errors="replace").read().splitlines()
hdr = raw[0].split("\t")
print("列数:", len(hdr), "| col26 表头:", repr(hdr[26]) if len(hdr) > 26 else "?")
c = collections.Counter()
for ln in raw[1:]:
    p = ln.split("\t")
    if len(p) > 26:
        c[p[26]] += 1
print("不同取值数:", len(c), "| 总行:", sum(c.values()))
print()
print("含 plant 的原值:")
tot = 0
for k, v in sorted(c.items(), key=lambda x: (-x[1], x[0])):
    if "lant" in k:
        print("   %-70s %d" % (k[:70], v))
        tot += v
print("   合计:", tot)
print()
print("出现频次 TOP 12 取值:")
for k, v in c.most_common(12):
    print("   %-70s %d" % (k[:70], v))
print()
print("含 '(S)' 的原值 (前 15 个):")
n = 0
for k, v in sorted(c.items(), key=lambda x: -x[1]):
    if "(S)" in k:
        print("   %-70s %d" % (k[:70], v))
        n += 1
        if n >= 15:
            break
