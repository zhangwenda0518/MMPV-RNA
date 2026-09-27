#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""诊断：24 个目标文件的换行风格（读备份）+ 当前文件是否被错误追加在 \\r 之后。"""
import glob, os

ROOT = os.path.expanduser("~/MMPV-paper")
TARGETS = [
    "09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv",
    "09_Virome_Analysis/HQ_analysis/HQ_plant_viruses_info.tsv",
    "09_Virome_Analysis/all_plant_analysis/all_plant_viruses_genus_summary.tsv",
]
files = []
for rel in TARGETS:
    for pat in (os.path.join(ROOT, "*", "02_novel_virus", "*_out", rel), os.path.join(ROOT, "*", "onekp-virus", rel)):
        files += sorted(glob.glob(pat))

print("%-58s %6s %6s %6s %6s %6s" % ("文件", "CRLF", "LF", "CR", "原列", "状态"))
bad = []
for p in files:
    bak = p + ".bak_dnarna_20260830"
    src = bak if os.path.exists(bak) else p
    raw = open(src, "rb").read()
    crlf = raw.count(b"\r\n"); lf = raw.count(b"\n") - crlf; cr = raw.count(b"\r") - crlf
    cur = open(p, "rb").read()
    # 当前文件里是否出现 “\r\t...Genome_Type” 这种错误插入（说明追加落在 \r 之后）
    mangled = b"\r\tGenome_" in cur or b"\r\tDNA" in cur or b"\r\tRNA" in cur or b"\r\t" in cur.split(b"\n")[0]
    ncol = len(cur.split(b"\n")[0].split(b"\t"))
    tag = "OK" if not mangled else "!!追加落在\\r后"
    if mangled:
        bad.append(p)
    print("%-58s %6d %6d %6d %6d %6s" % (p.split("MMPV-paper/")[-1][-58:], crlf, lf, cr, ncol, tag))

print("\n受影响文件数=%d" % len(bad))
for p in bad:
    print("  " + p)
