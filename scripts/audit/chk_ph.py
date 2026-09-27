#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""扫描 8 个阶元列里残留的占位值（限定词子串 + 明确的占位整串），
用于判断「占位值清理只覆盖 Family/Genus/Species」的漏网影响面。
用法: python3 chk_ph.py P1.tsv [P2.tsv ...]
"""
import csv, io, re, sys, os
from collections import defaultdict

TAX = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
NA_TOKENS = {"", "na", "n/a", "nan", "null", "none", "undefined"}
QUAL = re.compile("environmental|uncultured|unclassified|unidentified|unassigned|not assigned|no hit|unknown|undefined")
EXACT = {"virus", "viruses", "other", "samples", "sample", "unclassified sequences", "root"}


def norm(v):
    v = (v or "").strip()
    return "" if v.lower() in NA_TOKENS else v


for path in sys.argv[1:]:
    if not os.path.isfile(path):
        print("缺失 %s" % path)
        continue
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        raw = f.read()
    lines = [ln for ln in raw.splitlines() if not ln.startswith("#")]
    rd = csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t")
    hits = defaultdict(lambda: defaultdict(int))
    for r in rd:
        for k in TAX:
            v = norm(r.get(k))
            if not v:
                continue
            lv = v.lower()
            if QUAL.search(lv) or lv in EXACT:
                hits[k][v] += 1
    print("\n=== %s" % path)
    tot = 0
    for k in TAX:
        n = sum(hits[k].values())
        tot += n
        if n:
            top = sorted(hits[k].items(), key=lambda kv: -kv[1])[:4]
            print("  %-9s %5d 格  例: %s" % (k, n, "; ".join("%s(%d)" % t for t in top)))
    if tot == 0:
        print("  8 个阶元全部无占位残留")
