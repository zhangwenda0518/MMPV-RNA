#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用真实 rescue_report.tsv 验证修复后的统计口径"""
import collections

TSV = "/home/zhangwenda/data-test/out/08_Rescue/Plant/rescue_report.tsv"

total = 0
branch = collections.Counter()
notes = collections.Counter()

with open(TSV) as f:
    header = f.readline()
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 10:
            continue
        total += 1
        branch[parts[3]] += 1
        if parts[3] == 'fail':
            notes[parts[9]] += 1

failed = sum(notes.values())

print("=== 修复后口径 ===")
print(f"total (原始 centroids)     = {total}")
print(f"failed (fail_reasons 累加) = {failed}")
print(f"fail 占比                  = {failed/total*100:.1f}%")
print("")
print("=== 未拯救原因分布 (分母 = failed) ===")
print(f"{'原因':<32} {'数量':>8} {'占比':>8}")
for k, v in notes.most_common():
    pct = v / failed * 100
    flag = "  <-- 超100%!" if pct > 100 else ""
    print(f"{k:<32} {v:>8} {pct:>7.1f}%{flag}")
print(f"{'合计':<32} {failed:>8} {failed/failed*100:>7.1f}%")

print("")
print("=== 校验 ===")
print(f"原因加总 == failed: {sum(notes.values()) == failed}")
print(f"所有占比 <= 100%:   {all(v/failed*100 <= 100 for v in notes.values())}")
