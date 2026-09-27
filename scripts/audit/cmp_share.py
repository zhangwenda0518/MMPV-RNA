#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cmp_share.py —— 把各票首门槛下的最终产物与 0.5 基准逐格比对

输出：每个门槛相对基准的「不同行数 / 逐列不同格数」，用来判断这个门槛到底值多少行。
"""
import csv
import os
import sys

BASE = "/tmp/rc64b_cascade/final_integrated_classification.tsv"
VARIANTS = [("0.3", "/tmp/share_0.3/final_integrated_classification.tsv"),
            ("0.4", "/tmp/share_0.4/final_integrated_classification.tsv"),
            ("0.6", "/tmp/share_0.6/final_integrated_classification.tsv"),
            ("0.7", "/tmp/share_0.7/final_integrated_classification.tsv")]


def load(p):
    with open(p, "r", encoding="utf-8", newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        return rdr.fieldnames, {r["contig_id"]: r for r in rdr}


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    if not os.path.exists(BASE):
        print("缺基准: %s" % BASE)
        return
    cols, base = load(BASE)
    print("基准(门槛 0.50): %s 行 %d" % (BASE, len(base)))
    for tag, p in VARIANTS:
        if not os.path.exists(p):
            print("  [跳过] 缺 %s" % p)
            continue
        vcols, var = load(p)
        use = [c for c in cols if c in vcols and c != "contig_id"]
        rows_diff = set()
        coldiff = {c: 0 for c in use}
        for cid, b in base.items():
            v = var.get(cid)
            if v is None:
                rows_diff.add(cid)
                continue
            for c in use:
                if (b.get(c) or "").strip() != (v.get(c) or "").strip():
                    coldiff[c] += 1
                    rows_diff.add(cid)
        tot = sum(coldiff.values())
        print("  门槛 %s : 不同行 %d (%.2f%% of %d) ; 不同格总数 %d"
              % (tag, len(rows_diff), 100.0 * len(rows_diff) / max(1, len(base)), len(base), tot))
        hot = sorted([(n, c) for c, n in coldiff.items() if n > 0], reverse=True)
        for n, c in hot:
            print("        %-22s %6d 格" % (c, n))
        if not hot:
            print("        与基准逐格相同")
        print()


if __name__ == "__main__":
    main()
