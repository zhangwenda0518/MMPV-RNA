#!/usr/bin/env python3
"""统计 NCBI rankedlineage.dmp 里出现的 realm 取值与条数, 用于核对 KNOWN_REALMS 白名单覆盖面。

用法: python3 probe_realm_inventory.py [rankedlineage.dmp]
"""
import sys, collections

p = sys.argv[1] if len(sys.argv) > 1 else "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
cnt = collections.Counter()
with open(p, encoding="utf-8", errors="replace") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t|\t")
        if len(parts) < 10:
            continue
        r = parts[9].strip()
        if r:
            cnt[r] += 1
print("realm 取值共 %d 种:" % len(cnt))
for v, c in cnt.most_common():
    print("   %-20s %d" % (v, c))
