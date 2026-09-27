#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""口径 A 的植物病毒成果影响面：54 行置空具体是谁。"""
import csv
import os
from collections import Counter

BASE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
CAL = os.path.join(BASE, "05_Taxonomy", "Votus.integrated", "calibration_20260914")


def p(*a):
    print(*a, flush=True)


def rd(path, sep="\t"):
    with open(path, encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter=sep))


act = {}
for d in rd(os.path.join(CAL, "final_integrated_classification.calibrated.tsv")):
    act[(d.get("contig_id") or "").strip().strip('"')] = (d.get("calib_action") or "").strip().strip('"')

p("=== 置空行在矫正前的 科/属/种 构成（全表口径）===")
rows = rd(os.path.join(BASE, "05_Taxonomy", "Votus.integrated", "final_integrated_classification.tsv"))
blank = [r for r in rows if act.get((r.get("contig_id") or "").strip().strip('"')) == "blank"]
p("置空 %d 行" % len(blank))
c = Counter()
for r in blank:
    c["%s | %s | %s" % ((r.get("Family") or "").strip().strip('"'),
                        (r.get("Genus") or "").strip().strip('"'),
                        (r.get("Species") or "").strip().strip('"'))] += 1
for k, v in c.most_common(18):
    p("   %-72s %4d" % (k, v))

p("\n=== 这 54 行在植物病毒成果表里的原样 ===")
pf = os.path.join(BASE, "10_Reports", "All_plant.viruses_info.tsv")
prows = rd(pf)
hdr = list(prows[0].keys())
p("列: %s" % hdr)
hit = [r for r in prows if act.get((r.get("contig_id") or "").strip()) == "blank"]
p("命中 %d / %d 行" % (len(hit), len(prows)))
cg = Counter()
for r in hit:
    cg["%s | %s | %s" % ((r.get("Family") or "").strip(), (r.get("Genus") or "").strip(),
                         (r.get("Species") or "").strip())] += 1
for k, v in cg.most_common(30):
    p("   %-72s %4d" % (k, v))

p("\n-- 逐行明细（前 25，含 primary_tool / confidence）--")
for r in hit[:25]:
    p("   %-56s %-22s %-18s %-6s %s" % (r["contig_id"][:56], (r.get("Family") or "")[:22],
                                        (r.get("Genus") or "")[:18],
                                        (r.get("primary_tool") or ""), (r.get("confidence") or "")))

p("\n=== 属汇总表会怎么变 ===")
gs = os.path.join(BASE, "10_Reports", "all_plant_viruses_genus_summary.tsv")
if os.path.isfile(gs):
    for d in rd(gs):
        if (d.get("Genus") or "").strip() in {k.split(" | ")[1] for k in cg}:
            p("   %-18s n_contigs=%-5s n_rescued=%-5s" % (d.get("Genus"), d.get("n_contigs"), d.get("n_rescued")))
