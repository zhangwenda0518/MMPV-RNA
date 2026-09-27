#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""465 行与植物病毒成果表的交集 + 小汇总表内容直读。"""
import csv
import os

BASE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
CAL = os.path.join(BASE, "05_Taxonomy", "Votus.integrated", "calibration_20260914")


def p(*a):
    print(*a, flush=True)


def rd(path, sep="\t"):
    with open(path, encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter=sep))


b465 = set(l.strip() for l in open(os.path.join(CAL, "biavirus_arbitration", "selected_ids.txt"), encoding="utf-8") if l.strip())
b1246 = set()
for d in rd(os.path.join(CAL, "final_integrated_classification.calibrated.tsv")):
    if (d.get("calib_action") or "").strip().strip('"') == "blank":
        b1246.add((d.get("contig_id") or "").strip().strip('"'))
p("465 集 %d 条 / 置空集 %d 条 / 交集 %d 条" % (len(b465), len(b1246), len(b465 & b1246)))

p("\n=== 植物病毒成果表里的交集 ===")
for rel in ["10_Reports/All_plant.viruses_info.tsv",
            "10_Reports/plant_final_taxonomy.tsv",
            "09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv",
            "10_Reports/plant_virus_summary.tsv"]:
    fp = os.path.join(BASE, rel)
    if not os.path.isfile(fp):
        p("   [缺] %s" % rel)
        continue
    rows = rd(fp)
    key = "contig_id" if rows and "contig_id" in rows[0] else None
    if key is None:
        p("   [无 contig_id] %s  列=%s" % (rel, list(rows[0].keys())[:12] if rows else "-"))
        continue
    i465 = [r for r in rows if (r.get("contig_id") or "").strip() == "" or (r["contig_id"].strip() in b465)]
    i465 = [r for r in rows if r["contig_id"].strip() in b465]
    i1246 = [r for r in rows if r["contig_id"].strip() in b1246]
    p("   %-52s 总 %5d | 465集 %3d | 置空集 %3d" % (rel, len(rows), len(i465), len(i1246)))
    for r in i465[:8]:
        p("        %-52s Fam=%-22s Gen=%-14s Sp=%s" % (r["contig_id"][:52],
          (r.get("Family") or "").strip(), (r.get("Genus") or "").strip(), (r.get("Species") or "").strip()))

p("\n=== 小汇总表全文 ===")
for rel in ["10_Reports/HQ_plant_viruses_info.tsv",
            "10_Reports/all_plant_viruses_genus_summary.tsv"]:
    fp = os.path.join(BASE, rel)
    if os.path.isfile(fp):
        p("\n-- %s --" % rel)
        with open(fp, encoding="utf-8", errors="replace") as fh:
            for i, line in enumerate(fh):
                if i > 14:
                    p("   ...")
                    break
                p("   " + line.rstrip()[:220])

p("\n=== 465 集里落在 06_HostPrediction 的宿主预测 ===")
hp = os.path.join(BASE, "06_HostPrediction", "ensemble_host_summary.tsv")
if os.path.isfile(hp):
    rows = [r for r in rd(hp) if (r.get("contig_id") or "").strip() in b465]
    from collections import Counter
    p("   命中 %d 条；宿主列 = %s" % (len(rows), list(rows[0].keys())[:14] if rows else "-"))
    if rows:
        for c in list(rows[0].keys()):
            vals = Counter((r.get(c) or "").strip() for r in rows)
            if 1 < len(vals) <= 6:
                p("      %-22s %s" % (c, dict(vals.most_common(6))))
