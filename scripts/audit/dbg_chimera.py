#!/usr/bin/env python3
"""诊断：同一份表、同一份 dmp，两种实现为何给出不同矛盾数。
打印 chimera_multi.py 口径下被判矛盾的实例（含参照条目），供人工核对。
用法: python3 /tmp/dbg_chimera.py <tsv> [limit]
"""
import csv
import sys
from collections import Counter

RANKED = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
RANKS = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm"]
REF_IDX = {"Species": 2, "Genus": 3, "Family": 4, "Order": 5,
           "Class": 6, "Phylum": 7, "Kingdom": 8, "Realm": 9}
PAIRS = [("Realm", "Kingdom", "Realm"), ("Kingdom", "Phylum", "Kingdom"),
         ("Phylum", "Class", "Phylum"), ("Class", "Order", "Class"),
         ("Order", "Family", "Order"), ("Family", "Genus", "Family"),
         ("Genus", "Species", "Genus")]

TSV = sys.argv[1]
LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else 15


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None
    return v


rows = []
need = set()
with open(TSV, newline="", encoding="utf-8", errors="replace") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        d = {k.strip(): norm(v) for k, v in r.items()}
        rows.append(d)
        for L in RANKS:
            if d.get(L):
                need.add(d[L].lower())
print("行数 %d | need name %d" % (len(rows), len(need)))

ref = {}
with open(RANKED, encoding="utf-8", errors="replace") as fh:
    for line in fh:
        parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
        if len(parts) < 10:
            continue
        nm = parts[1].lower()
        if nm in need and nm not in ref:
            ref[nm] = {L: parts[REF_IDX[L]] for L in RANKS}
print("ref 命中 %d" % len(ref))

shown = Counter()
cnt = Counter()
for i, r in enumerate(rows):
    for h, l, refcol in PAIRS:
        hv, lv = r.get(h), r.get(l)
        if not hv or not lv:
            continue
        e = ref.get(lv.lower())
        if e is None:
            continue
        pv = (e.get(refcol) or "").strip()
        if pv and pv != hv:
            cnt["%s-%s" % (h, l)] += 1
            key = "%s-%s" % (h, l)
            if shown[key] < 3:
                shown[key] += 1
                print("MISS row=%d %s: 表内 %s=%s %s=%s | ref 中 %s 的 %s=%s"
                      % (i + 1, key, h, hv, l, lv, lv, refcol, pv))
print("---")
for k, v in cnt.items():
    print("%s %d" % (k, v))
