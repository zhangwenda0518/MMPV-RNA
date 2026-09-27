#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""打印指定 contig 在 05 表中的完整分类行 (含 Species / confidence / primary_tool)。"""
import csv
import sys

path = sys.argv[1]
ids = set(sys.argv[2:])
with open(path, encoding="utf-8", errors="replace", newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    names = r.fieldnames
    key = names[0]
    cols = [c for c in ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family",
                        "Genus", "Species", "confidence", "primary_tool"]
            if c in names]
    for row in r:
        cid = (row.get(key) or "").strip()
        if cid in ids:
            print("%s | %s" % (cid, " | ".join(
                "%s=%s" % (c, (row.get(c) or "NA").strip()) for c in cols)))
