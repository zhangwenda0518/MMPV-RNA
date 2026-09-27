#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""1) 列出 05_Taxonomy/Votus.integrated 目录里的证据文件
   2) 抽 3 条冲突 contig，在 05/09/10_Reports 各层里找它的行，并在 Votus 证据文件里搜它的记录
"""
import os, glob, csv

ROOT = os.path.expanduser("~/MMPV-paper")
BASE = os.path.join(ROOT, "goji-virome/02_novel_virus/RNA-Lycium_barbarum_out")
CONTIGS = ["CRR527041_clean_NODE_1593_length_1853_cov_4.065169", "contig_1940", "contig_312"]

print("=== 1) 05_Taxonomy/Votus.integrated 目录 ===")
for d in sorted(glob.glob(os.path.join(BASE, "05_Taxonomy", "*"))):
    print("  [dir] %s" % d.split("05_Taxonomy/")[-1])
for d in sorted(glob.glob(os.path.join(BASE, "05_Taxonomy", "*"))):
    if os.path.isdir(d):
        for f in sorted(os.listdir(d)):
            p = os.path.join(d, f)
            if os.path.isfile(p):
                sz = os.path.getsize(p)
                try:
                    with open(p, "r", encoding="utf-8", errors="replace") as fh:
                        h = fh.readline().rstrip("\n")
                        n = sum(1 for _ in fh)
                    print("    %-58s %9d B  行=%-7d %s" % (f, sz, n, h[:110]))
                except Exception as e:
                    print("    %-58s %9d B  (读失败 %s)" % (f, sz, e))

print("\n=== 2) 冲突 contig 在各层的行 ===")
targets = [
    ("05 final_integrated_classification", "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"),
    ("09 All_plant.viruses_info", "09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv"),
    ("10_Reports final_integrated", "10_Reports/final_integrated_classification.tsv"),
    ("10_Reports All_plant.viruses_info", "10_Reports/All_plant.viruses_info.tsv"),
]
for label, rel in targets:
    q = os.path.join(BASE, rel)
    if not os.path.exists(q):
        print("  %-34s 不存在" % label); continue
    with open(q, newline="", encoding="utf-8", errors="surrogateescape") as f:
        rd = csv.reader(f, delimiter="\t"); hdr = next(rd)
        hits = [r for r in rd if r and any(c.startswith(CONTIGS[0][:20]) for c in r[:2]) or (r and r[0] in CONTIGS[1:])]
    print("  %-34s 列=%d  命中=%d" % (label, len(hdr), len(hits)))
    for h in hits:
        print("      " + " | ".join(h[:13] + h[-6:]))

print("\n=== 3) 在 05_Taxonomy 全部文件里搜这 3 条 contig ===")
for f in sorted(glob.glob(os.path.join(BASE, "05_Taxonomy", "**", "*"), recursive=True)):
    if not os.path.isfile(f) or os.path.getsize(f) > 200 * 1024 * 1024:
        continue
    try:
        txt = open(f, "rb").read()
    except Exception:
        continue
    for c in CONTIGS:
        if c.encode() in txt:
            print("  命中 %-52s <- %s" % (c, f.split("MMPV-paper/")[-1]))
