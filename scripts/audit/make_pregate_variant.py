#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从 v6.4 补丁生成「闸门关闭」变体，用于观察 cascade 定音后的裸共识状态。

只注释掉三道后置闸门的调用，不动共识引擎本身：
    harmonize_family_genus / enforce_rank_containment / enforce_species_containment
目的：量出「逐级淘汰结束后、相容性闸门介入前」还残留多少行科属不配套。
"""
import hashlib
import os
import sys

SRC = sys.argv[1] if len(sys.argv) > 1 else "/tmp/virus_classifier_analysis.R.cascade_v64"
DST = sys.argv[2] if len(sys.argv) > 2 else "/tmp/virus_classifier_analysis.R.pregate"

CALLS = [
    "  wide <- harmonize_family_genus(wide, stacked, tool_weights)\n",
    "  wide <- enforce_rank_containment(wide, GENUS_FAMILY_REF)\n",
    "  wide <- enforce_species_containment(wide, GENUS_FAMILY_REF, out_dir)\n",
]


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


src = open(SRC, encoding="utf-8").read()
out = src
for c in CALLS:
    n = out.count(c)
    if n != 1:
        sys.exit("调用点计数异常 %d: %r" % (n, c))
    out = out.replace(c, "  # [PREGATE-DIAG] " + c.lstrip())

open(DST, "w", encoding="utf-8").write(out)
print("源  md5 %s  %d 字节" % (md5(SRC), os.path.getsize(SRC)))
print("产物 %s  md5 %s  %d 字节" % (DST, md5(DST), os.path.getsize(DST)))
print("已注释调用点 %d 个" % len(CALLS))
