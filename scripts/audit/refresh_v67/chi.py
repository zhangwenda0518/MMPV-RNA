#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""独立复核：相邻阶元参照相容性（矛盾率），判据对齐 /tmp/chimera_multi.py。
用法: python3 chi.py <tag1:path1> [<tag2:path2> ...]
对每个文件输出：行数、7 对相邻阶元矛盾数、矛盾行数、矛盾率。
参照库: rankedlineage.dmp  列序 tax_id|tax_name|species|genus|family|order|class|phylum|kingdom|superkingdom
"""
import csv
import sys
from collections import Counter

RANKED = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
RANKS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
REF_IDX = {"Species": 2, "Genus": 3, "Family": 4, "Order": 5,
           "Class": 6, "Phylum": 7, "Kingdom": 8, "Realm": 9}
PAIRS = [
    ("Realm", "Kingdom", "Realm"),
    ("Kingdom", "Phylum", "Kingdom"),
    ("Phylum", "Class", "Phylum"),
    ("Class", "Order", "Class"),
    ("Order", "Family", "Order"),
    ("Family", "Genus", "Family"),
    ("Genus", "Species", "Genus"),
]


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
    rows = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append({k.strip(): norm(v) for k, v in r.items()})
    return rows


def main():
    inputs = []
    for arg in sys.argv[1:]:
        tag, _, p = arg.partition(":")
        inputs.append((tag, p))
    tables = [(tag, p, load(p)) for tag, p in inputs]

    need = set()
    for _, _, rows in tables:
        for r in rows:
            for L in RANKS:
                if r.get(L):
                    need.add(r[L].lower())

    ref = {}
    with open(RANKED, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
            if len(parts) < 10:
                continue
            nm = parts[1].lower()
            if nm in need and nm not in ref:
                ref[nm] = {L: parts[REF_IDX[L]] for L in RANKS}
    print("REF_HITS\t%d\tNEED\t%d" % (len(ref), len(need)))

    for tag, p, rows in tables:
        bc = Counter()
        nconf = 0
        for r in rows:
            bad = 0
            for h, l, refcol in PAIRS:
                hv, lv = r.get(h), r.get(l)
                if not hv or not lv:
                    continue
                e = ref.get(lv.lower())
                if e is None:
                    continue
                pv = (e.get(refcol) or "").strip()
                if pv and pv != hv:
                    bc["%s-%s" % (h, l)] += 1
                    bad = 1
            nconf += bad
        print("---\t%s\t%s" % (tag, p))
        print("ROWS\t%d" % len(rows))
        for h, l, _ in PAIRS:
            print("PAIR\t%s-%s\t%d" % (h, l, bc["%s-%s" % (h, l)]))
        print("CONF_ROWS\t%d" % nconf)
        print("CONF_RATE\t%.4f%%" % (100.0 * nconf / len(rows)))


if __name__ == "__main__":
    main()
