#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""10_Reports / rescue 报告层变化清单（独立复核用，只读）。

用法: python3 audit_reports.py <数据集目录>
"""
import filecmp
import hashlib
import os
import sys
from datetime import datetime

D = sys.argv[1]
BAK = os.path.join(D, "10_Reports.bak_v66f_20260916")
NEW = os.path.join(D, "10_Reports")


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def ts(p):
    return datetime.fromtimestamp(os.path.getmtime(p)).strftime("%Y-%m-%d %H:%M:%S")


allb = sorted(os.listdir(BAK))
alln = sorted(os.listdir(NEW))
only_b = [f for f in allb if f not in set(alln)]
only_n = [f for f in alln if f not in set(allb)]
print("bak 文件数 %d / new 文件数 %d" % (len(allb), len(alln)))
print("only_in_bak: %s" % (only_b or "-"))
print("only_in_new: %s" % (only_n or "-"))
print()

changed = []
for f in allb:
    pb, pn = os.path.join(BAK, f), os.path.join(NEW, f)
    if not os.path.isfile(pb) or not os.path.isfile(pn):
        continue
    if not filecmp.cmp(pb, pn, shallow=False):
        changed.append(f)

print("内容变化文件 %d 个:" % len(changed))
for f in changed:
    pb, pn = os.path.join(BAK, f), os.path.join(NEW, f)
    print("  %-42s bak[mtime=%s size=%7d md5=%s]  new[mtime=%s size=%7d md5=%s]" %
          (f, ts(pb), os.path.getsize(pb), md5(pb)[:8],
           ts(pn), os.path.getsize(pn), md5(pn)[:8]))
print()
print("内容不变文件 %d 个:" % (len(allb) - len(changed) - len(only_b)))
for f in allb:
    if f not in changed and f not in only_b:
        print("  %s" % f)

# rescue
print()
print("=== rescue 报告层 ===")
res = [
    (os.path.join(D, "08_Rescue", "Plant", "rescue_report.tsv"), "08_Rescue/Plant/rescue_report.tsv"),
    (os.path.join(D, "08_Rescue", "Plant", "rescue_summary.md"), "08_Rescue/Plant/rescue_summary.md"),
    (os.path.join(D, "08_Rescue", "Plant", "rescue_report.tsv.bak_v661"), "08_Rescue/Plant/rescue_report.tsv.bak_v661"),
    (os.path.join(D, "08_Rescue", "Plant", "rescue_summary.md.bak_v661"), "08_Rescue/Plant/rescue_summary.md.bak_v661"),
    (os.path.join(NEW, "rescue_report.tsv"), "10_Reports/rescue_report.tsv"),
]
for p, label in res:
    if os.path.exists(p):
        print("  %-46s mtime=%s size=%6d md5=%s" % (label, ts(p), os.path.getsize(p), md5(p)))
    else:
        print("  %-46s MISSING" % label)
