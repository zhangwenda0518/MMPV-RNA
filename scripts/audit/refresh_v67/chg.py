#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""独立复核：各阶元取值变化行数 + Order/Family 有值->空 行数。
用法: python3 chg.py <old.tsv> <new.tsv>
"""
import csv
import sys
from collections import Counter

RANKS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None
    return v


def load(path):
    d = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = (r.get("contig_id") or "").strip().strip('"')
            d[k] = {c.strip(): norm(v) for c, v in r.items()}
    return d


old = load(sys.argv[1])
new = load(sys.argv[2])
keys = [k for k in old if k in new]
print("COMMON\t%d\told=%d\tnew=%d" % (len(keys), len(old), len(new)))
chg = Counter()
lost = Counter()
for k in keys:
    o, n = old[k], new[k]
    for L in RANKS:
        if (o.get(L) or "") != (n.get(L) or ""):
            chg[L] += 1
    for L in ("Order", "Family"):
        if o.get(L) and not n.get(L):
            lost[L] += 1
        elif not o.get(L) and n.get(L):
            lost[L + "_gain"] += 1
for L in RANKS:
    print("CHG\t%s\t%d" % (L, chg[L]))
print("LOST\tOrder\t%d" % lost["Order"])
print("LOST\tFamily\t%d" % lost["Family"])
print("GAIN\tOrder\t%d" % lost["Order_gain"])
print("GAIN\tFamily\t%d" % lost["Family_gain"])
