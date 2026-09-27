#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v6.6 vs v6.7 成品对照：逐阶元差异方向（有->空 / 空->有 / 值变）+ 非空行数变化。"""
import csv
import os

LV = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
EMPTY = {"", "na", "n/a", "nan", "-", "none", "null"}


def nz(x):
    x = (x or "").strip()
    return "" if x.lower() in EMPTY else x


A = os.environ.get("A", "/tmp/rc66_count/final_integrated_classification.tsv")
B = os.environ.get("B", "/tmp/rc67_count/final_integrated_classification.tsv")


def load(p):
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        return {r["contig_id"].strip(): r for r in csv.DictReader(f, delimiter="\t")}


a, b = load(A), load(B)
print("%s\nvs %s" % (A, B))
print("行数 %d / %d" % (len(a), len(b)))
print()
print("| 阶元 | 非空 v66 | 非空 v67 | 变动行 | 有->空 | 空->有 | 值变 |")
print("|---|---|---|---|---|---|---|")
tot = set()
for lv in LV:
    ne_a = sum(1 for c in a if nz(a[c].get(lv)))
    ne_b = sum(1 for c in b if nz(b[c].get(lv)))
    down = up = chg = 0
    for c in a:
        x = nz(a[c].get(lv))
        y = nz(b.get(c, {}).get(lv))
        if x == y:
            continue
        tot.add(c)
        if x and not y:
            down += 1
        elif y and not x:
            up += 1
        else:
            chg += 1
    print("| %s | %d | %d | %d | %d | %d | %d |" % (lv, ne_a, ne_b, down + up + chg, down, up, chg))
print()
print("涉及 contig 合计 %d" % len(tot))
print()
print("变化行全清单（%d 行 × 8 阶元，只列有差异的阶元）" % len(tot))
for c in sorted(tot):
    parts = []
    for lv in LV:
        x = nz(a[c].get(lv))
        y = nz(b.get(c, {}).get(lv))
        if x != y:
            parts.append("%s: [%s] -> [%s]" % (lv, x or "空", y or "空"))
    print("  %s" % c)
    for p in parts:
        print("      %s" % p)
