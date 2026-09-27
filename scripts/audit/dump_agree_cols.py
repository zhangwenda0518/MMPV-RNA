#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""打印指定 contig 的逐 rank 一致度列 (_agree), 用于判断哪个 rank 的标签有多工具支持。"""
import csv
import sys

path = sys.argv[1]
ids = set(sys.argv[2:])
want = ["Realm", "Kingdom", "Family", "Family_agree", "Genus", "Genus_agree",
        "Species", "Species_agree", "confidence", "primary_tool"]
with open(path, encoding="utf-8", errors="replace", newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    names = r.fieldnames or []
    cols = [c for c in want if c in names]
    print("可用列: %s" % [c for c in names if "agree" in c])
    for row in r:
        if (row.get("contig_id") or "").strip() in ids:
            print("  %s\n     %s" % (row["contig_id"], " | ".join(
                "%s=%s" % (c, (row.get(c) or "").strip()) for c in cols)))
