#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""成品层「世系自洽性」体检：成品里的每一格取值，是否真被某个工具的原始报告支持。

动机：逐级投票是「每个阶元各自定值」，如果某个阶元的票首来自 A 阵营、下一阶元的票首来自 B 阵营，
成品就会出现一条没有任何工具报过的拼接世系（例如 Realm=巨型 DNA 病毒、Family=人丙肝病毒）。
闸门只管科-属、科-种两组配对，管不到这种跨阶元错配。

判据：对成品每一行，取相邻阶元对 (L_i, L_i+1)，检查工具原始报告里是否存在
      某个工具、在同一个 contig 上、同时报出这两个值。存在=有出处，不存在=拼接。
用法: python3 check_chimera.py <combined.tsv> <成品.tsv> [成品2.tsv ...]
"""
import csv
import sys
from collections import Counter

LV = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
PAIRS = [(LV[i], LV[i + 1]) for i in range(len(LV) - 1)]
EMPTY = {"", "na", "n/a", "nan", "-", "none", "null", "unclassified", "unknown"}


def nz(x):
    if x is None:
        return None
    x = x.strip().strip('"').strip()
    return None if x.lower() in EMPTY else x


def read_tsv(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


combined, products = sys.argv[1], sys.argv[2:]

# 工具原始报告：每个相邻阶元对，建 (contig, 上值, 下值) 的支持集合
support = {p: set() for p in PAIRS}
for r in read_tsv(combined):
    c = nz(r.get("seq_name"))
    if not c:
        continue
    for a, b in PAIRS:
        va, vb = nz(r.get(a)), nz(r.get(b))
        if va and vb:
            support[(a, b)].add((c, va, vb))

print("工具原始报告: %s" % combined)
print("支持集合规模（相邻阶元对 -> 不同组合数）:")
for p in PAIRS:
    print("  %-20s %d" % ("%s-%s" % p, len(support[p])))
print()

hdr = ["成品", "行数", "有值对数"] + ["%s-%s" % p for p in PAIRS] + ["任意一对不支持的行数"]
print("| " + " | ".join(hdr) + " |")
print("|" + "---|" * len(hdr))

for prod in products:
    rows = read_tsv(prod)
    bad_by_pair = Counter()
    tot_by_pair = Counter()
    bad_rows = set()
    for r in rows:
        c = nz(r.get("contig_id"))
        if not c:
            continue
        for a, b in PAIRS:
            va, vb = nz(r.get(a)), nz(r.get(b))
            if not (va and vb):
                continue
            tot_by_pair[(a, b)] += 1
            if (c, va, vb) not in support[(a, b)]:
                bad_by_pair[(a, b)] += 1
                bad_rows.add(c)
    cells = []
    for p in PAIRS:
        t, b = tot_by_pair[p], bad_by_pair[p]
        cells.append("%d/%d (%.1f%%)" % (b, t, 100.0 * b / t) if t else "-")
    name = prod.split("/")[-2] + "/" + prod.split("/")[-1]
    print("| %s | %d | %d | %s | %d |" % (
        name, len(rows), sum(tot_by_pair.values()), " | ".join(cells), len(bad_rows)))
