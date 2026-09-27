#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""列出旧植物表里 Family 属于三个巨病毒科的行 (含 contig_id), 便于逐票追踪。"""
import csv
import sys

TARGETS = {"Mimiviridae", "Marseilleviridae", "Pithoviridae"}

path = sys.argv[1]
with open(path, encoding="utf-8", errors="replace", newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    key = r.fieldnames[0]
    for row in r:
        fam = (row.get("Family") or "").strip()
        if fam in TARGETS:
            print("%s\t%s\t%s\t%s" % (row[key], fam, row.get("Genus"),
                                      row.get("Species")))
