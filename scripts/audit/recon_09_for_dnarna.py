#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""09 层 DNA/RNA 补列：只读侦察
列出 ~/MMPV-paper 下所有 09* 目录里的 .tsv，打印行数/列数/表头，判断是否有 Genus / Species 列可 join VMR。
不写任何文件。
"""
import csv, os, glob

ROOT = os.path.expanduser("~/MMPV-paper")
cands = []
for proj in sorted(os.listdir(ROOT)):
    pdir = os.path.join(ROOT, proj)
    if not os.path.isdir(pdir):
        continue
    for d in sorted(os.listdir(pdir)):
        if not (d.startswith("09") or d.startswith("10_Reports")):
            continue
        base = os.path.join(pdir, d)
        if not os.path.isdir(base):
            continue
        for dp, _, fns in os.walk(base):
            for fn in fns:
                if fn.lower().endswith((".tsv", ".txt")) and ("vir" in fn.lower() or "virus" in fn.lower() or "plant" in fn.lower()):
                    cands.append(os.path.join(dp, fn))

print("候选文件数 = %d" % len(cands))
for p in sorted(cands):
    try:
        with open(p, newline="", encoding="utf-8", errors="replace") as f:
            rd = csv.reader(f, delimiter="\t")
            hdr = next(rd)
            n = sum(1 for _ in rd)
    except Exception as e:
        print("!! %s : %s" % (p, e))
        continue
    cols = [h.strip().strip('"') for h in hdr]
    low = [c.lower() for c in cols]
    has_g = any("genus" in c for c in low)
    has_s = any(("species" in c) or ("virus_name" in c) or c == "virus" for c in low)
    print("\n%s" % p.replace(os.path.expanduser("~"), "~"))
    print("   行=%d 列=%d | Genus列=%s Species列=%s | size=%d" % (n, len(cols), has_g, has_s, os.path.getsize(p)))
    print("   header: %s" % " | ".join(cols[:26]))
