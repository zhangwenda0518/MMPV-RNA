#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""打印指定 contig 在 7 个工具里的逐阶元原始取值（矩阵视图）。"""
import csv
import glob
import os
import sys

OUT = os.environ.get("DIAG_OUT", "/tmp/pregate_cascade")
LEVELS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
EMPTY = {"", "na", "n/a", "nan", "-", "none", "null"}

targets = sys.argv[1:] or ["CRR1440126_clean_NODE_914_length_1792_cov_3.429319"]

vals = {}
tools = []
for p in sorted(glob.glob(os.path.join(OUT, "standardized_*.tsv"))):
    tool = os.path.basename(p)[len("standardized_"):-len(".tsv")]
    tools.append(tool)
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            vals.setdefault(r["contig_id"].strip(), {})[tool] = {
                lv: (r.get(lv) or "").strip() for lv in LEVELS
            }

for cid in targets:
    print("=" * 78)
    print(cid)
    print("=" * 78)
    d = vals.get(cid)
    if d is None:
        print("  该 contig 不在标准化表里")
        continue
    print("  %-9s %s" % ("阶元", "".join("%-22s" % t for t in tools)))
    for lv in LEVELS:
        cells = "".join("%-22s" % (d.get(t, {}).get(lv, "") or "·") for t in tools)
        print("  %-9s %s" % (lv, cells))
