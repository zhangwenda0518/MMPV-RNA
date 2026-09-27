#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""打印指定 contig 在 C9 classification_result.tsv 中的宿主判定关键字段。"""
import csv
import sys

path = sys.argv[1]
ids = set(sys.argv[2:])
cols = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
        "Species", "Predicted_Host", "Confidence_Level", "Determination_Level"]
with open(path, encoding="utf-8", errors="replace", newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        cid = (row.get("contig_id") or "").strip()
        if cid in ids:
            print("%s\n    %s" % (cid, " | ".join(
                "%s=%s" % (c, (row.get(c) or "").strip()) for c in cols)))
