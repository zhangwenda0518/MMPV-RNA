#!/usr/bin/env python3
"""统计 09b rescue_evidence_scored.tsv 的 verdict 分布与 hmm_rescue 生效数。
用法: count_verdict.py <scored.tsv> [label]
列位置按表头定位，不写死索引。"""
import csv
import sys
import collections

path = sys.argv[1]
label = sys.argv[2] if len(sys.argv) > 2 else path

with open(path, newline="", encoding="utf-8", errors="replace") as fh:
    rd = csv.reader(fh, delimiter="\t")
    header = next(rd)
    idx_v = None
    idx_h = None
    for i, h in enumerate(header):
        hl = h.strip().lower()
        if hl == "verdict":
            idx_v = i
        if hl == "hmm_rescue":
            idx_h = i
    cnt = collections.Counter()
    hmm_on = 0
    n = 0
    for row in rd:
        n += 1
        if idx_v is not None and idx_v < len(row):
            cnt[row[idx_v].strip() or "EMPTY"] += 1
        if idx_h is not None and idx_h < len(row):
            v = row[idx_h].strip().lower()
            if v and v not in ("", "na", "nan", "false", "0"):
                hmm_on += 1

print(f"== {label}  行数={n}  列数={len(header)}  hmm_rescue列={'有' if idx_h is not None else '无'}")
print("   verdict: " + ", ".join(f"{k}={v}" for k, v in cnt.most_common()))
print(f"   hmm_rescue 生效={hmm_on}")
keep = cnt.get("KEEP", 0)
rev = cnt.get("REVIEW", 0)
print(f"   KEEP+REVIEW={keep + rev}")
