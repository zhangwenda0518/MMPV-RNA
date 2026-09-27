#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""逐行对比两份 final_integrated_classification.tsv：列出改动量、改动分布与样本。"""
import argparse
import csv
import sys
from collections import Counter, defaultdict


def norm(v):
    if v is None:
        return ""
    v = str(v).strip()
    while v.startswith('"') and v.endswith('"') and len(v) >= 2:
        v = v[1:-1]
    return v.strip()


def load(path, key="contig_id"):
    rows = {}
    header = None
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        header = [c.strip() for c in (rd.fieldnames or [])]
        for r in rd:
            k = norm(r.get(key))
            if k:
                rows[k] = {c.strip(): norm(r.get(c)) for c in header}
    return header, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--sample", type=int, default=10)
    args = ap.parse_args()

    ho, ro = load(args.old)
    hn, rn = load(args.new)
    print("old: %s 行=%d  列=%d" % (args.old.split("/")[-1], len(ro), len(ho)))
    print("new: %s 行=%d  列=%d" % (args.new.split("/")[-1], len(rn), len(hn)))
    only_old = [k for k in ro if k not in rn]
    only_new = [k for k in rn if k not in ro]
    print("仅 old 有: %d ; 仅 new 有: %d" % (len(only_old), len(only_new)))
    if only_old[:3]:
        print("  only_old 例:", only_old[:3])
    if only_new[:3]:
        print("  only_new 例:", only_new[:3])

    cols = [c for c in hn if c in ho]
    changed = defaultdict(int)
    per_col_examples = defaultdict(list)
    row_changed = 0
    for k in ro:
        if k not in rn:
            continue
        a, b = ro[k], rn[k]
        diff_cols = [c for c in cols if norm(a.get(c)) != norm(b.get(c))]
        if diff_cols:
            row_changed += 1
            for c in diff_cols:
                changed[c] += 1
                if len(per_col_examples[c]) < args.sample:
                    per_col_examples[c].append((k, a.get(c), b.get(c)))
    print("")
    print("=== 变动统计（共同 contig %d）===" % len([k for k in ro if k in rn]))
    print("  有任意字段变动的行: %d (%.1f%%)" % (row_changed, 100.0 * row_changed / max(len(ro), 1)))
    for c in cols:
        if changed[c]:
            print("    %-16s %6d" % (c, changed[c]))

    print("")
    for c in ["Family", "Genus", "Species", "Order", "Class"]:
        if not changed.get(c):
            continue
        print("=== %s 变动样本 ===" % c)
        for k, a, b in per_col_examples[c]:
            print("  %-52s %-24s -> %-24s" % (k[:52], (a or "-")[:24], (b or "-")[:24]))
        print("")

    if args.out:
        with open(args.out, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["contig_id", "column", "old", "new"])
            for k in ro:
                if k not in rn:
                    continue
                a, b = ro[k], rn[k]
                for c in cols:
                    if norm(a.get(c)) != norm(b.get(c)):
                        w.writerow([k, c, a.get(c), b.get(c)])
        print("变更明细已写: %s" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
