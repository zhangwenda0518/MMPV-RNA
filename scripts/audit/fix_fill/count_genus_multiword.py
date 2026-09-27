#!/usr/bin/env python3
"""统计 Genus 列含空格(物种级名/多词名混入 Genus)的行数, 并给出典型样例。

用法: python3 count_genus_multiword.py <product.tsv> [label]
"""
import csv, sys, collections

p = sys.argv[1]
label = sys.argv[2] if len(sys.argv) > 2 else p
NAV = {"", "na", "n/a", "-"}
n = 0
rows_total = 0
ex = collections.Counter()
with open(p, newline="", encoding="utf-8", errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        rows_total += 1
        g = (r.get("Genus") or "").strip()
        if g.lower() in NAV:
            continue
        if " " in g:
            n += 1
            ex[g] += 1
print("%-12s 行数=%d  Genus含空格=%d (%.2f%%)" % (label, rows_total, n, 100.0 * n / max(rows_total, 1)))
for v, c in ex.most_common(5):
    print("      %-52s x%d" % (v[:52], c))
