#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""独立复核：新旧产物 contig_id 集合是否完全相同（不是只比计数）。
用法: python3 setchk.py <A.tsv> <B.tsv>
"""
import csv
import sys


def ids(path):
    out = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out.append((r.get("contig_id") or "").strip().strip('"'))
    return out


a = ids(sys.argv[1])
b = ids(sys.argv[2])
sa, sb = set(a), set(b)
print("A rows=%d unique=%d  dup=%d" % (len(a), len(sa), len(a) - len(sa)))
print("B rows=%d unique=%d  dup=%d" % (len(b), len(sb), len(b) - len(sb)))
print("A_only=%d  B_only=%d" % (len(sa - sb), len(sb - sa)))
if sa - sb:
    print("A_ONLY_SAMPLE\t" + "\t".join(sorted(sa - sb)[:10]))
if sb - sa:
    print("B_ONLY_SAMPLE\t" + "\t".join(sorted(sb - sa)[:10]))
print("SET_EQUAL\t%s" % (sa == sb))
print("ORDER_IDENTICAL\t%s" % (a == b))
