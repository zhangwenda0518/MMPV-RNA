#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""诊断：对比 Alternaria All_plant 原表(备份) 与已追加版本，找出不满足前缀的行。"""
import csv

D = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Alternaria_alternata_out/09_Virome_Analysis/all_plant_analysis/"
P = D + "All_plant.viruses_info.tsv"
B = P + ".bak_dnarna_20260830"

def load(p):
    with open(p, encoding="utf-8", errors="surrogateescape", newline="") as f:
        return f.read().split("\n")

old, new = load(B), load(P)
print("old lines=%d new lines=%d" % (len(old), len(new)))
print("old tail repr: %r" % old[-1][:80])
print("new tail repr: %r" % new[-1][:80])
print("old line0 repr: %r" % old[0][:60])
print("CRLF? old has \\r\\n: %d ; new: %d" % (sum("\r" in l for l in old), sum("\r" in l for l in new)))

diff = []
for i, (a, b) in enumerate(zip(old, new)):
    if not b.startswith(a + "\t"):
        diff.append(i)
print("不满足前缀的行号: %s (共 %d)" % (diff[:10], len(diff)))
for i in diff[:5]:
    print("--- line %d" % i)
    print("  old[%d]: %r" % (i, old[i][:200]))
    print("  new[%d]: %r" % (i, new[i][:260]))
    if i + 1 < len(new):
        print("  new[%d]: %r" % (i + 1, new[i + 1][:200]))

# 解析层对比：原列是否逐字段一致
def parsed(p):
    with open(p, encoding="utf-8", errors="surrogateescape", newline="") as f:
        return list(csv.reader(f, delimiter="\t"))
ro, rn = parsed(B), parsed(P)
print("\n解析行数 old=%d new=%d" % (len(ro), len(rn)))
bad = 0
for i, (a, b) in enumerate(zip(ro, rn)):
    if a != b[:len(a)]:
        bad += 1
        if bad <= 3:
            print("字段不一致 line %d:\n  old=%r\n  new=%r" % (i, a[:6], b[:9]))
print("字段层不一致行数=%d" % bad)
print("列数 old=%d new=%d" % (len(ro[0]), len(rn[0])))
