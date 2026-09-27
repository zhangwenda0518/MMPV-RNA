#!/usr/bin/env python3
"""stage_breakdown.py — 从 stdin 读 VSI run.log, 输出各阶段累计耗时与占比
用法: python stage_breakdown.py < run.log"""
import re
import sys
import datetime

WRAP = {"extendOneScaffold", "growScaffoldWithAssembly", "checkCoverage"}
pat_stage = re.compile(r"^([A-Za-z][A-Za-z0-9]*):$")
pat_time = re.compile(r"^(Start|End) time: (\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)")

st = None
t0 = {}
dur = {}
cnt = {}
for line in sys.stdin:
    line = line.strip()
    m = pat_stage.match(line)
    if m:
        st = m.group(1)
        continue
    m = pat_time.match(line)
    if m and st:
        ts = datetime.datetime.strptime(m.group(2), "%Y-%m-%d %H:%M:%S").timestamp()
        if m.group(1) == "Start":
            t0[st] = ts
        elif st in t0:
            dur[st] = dur.get(st, 0.0) + ts - t0[st]
            cnt[st] = cnt.get(st, 0) + 1
            del t0[st]

total = sum(v for k, v in dur.items() if k not in WRAP)
print("%-30s %10s %7s %7s" % ("stage", "minutes", "share", "calls"))
for k, v in sorted(dur.items(), key=lambda kv: -kv[1]):
    share = 100.0 * v / total if k not in WRAP else 0.0
    print("%-30s %10.1f %6.1f%% %7d" % (k, v / 60, share, cnt[k]))
print("%-30s %10.1f" % ("TOTAL(leaf stages)", total / 60))
