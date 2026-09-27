#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""逐工具票源追踪: 对给定 contig 列表, 打印 7 个工具各自给的 (Family, Genus)。

用法: python3 tool_vote_trace.py <int_dir> <contig_id> [<contig_id> ...]
"""
import csv
import glob
import os
import sys


def main():
    int_dir, ids = sys.argv[1], set(sys.argv[2:])
    files = sorted(glob.glob(os.path.join(int_dir, "standardized_*.tsv")))
    print("工具文件: %s" % [os.path.basename(f) for f in files])
    for f in files:
        tool = os.path.basename(f)[len("standardized_"):-len(".tsv")]
        with open(f, encoding="utf-8", errors="replace", newline="") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            names = r.fieldnames or []
            key = names[0]
            for row in r:
                cid = (row.get(key) or "").strip().strip('"')
                if cid in ids:
                    fam = (row.get("Family") or "NA").strip().strip('"')
                    gen = (row.get("Genus") or "NA").strip().strip('"')
                    sp = (row.get("Species") or "").strip().strip('"')
                    print("  %-16s %-34s | %-18s | %-22s | %s" % (tool, cid[:34], fam, gen, sp[:40]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
