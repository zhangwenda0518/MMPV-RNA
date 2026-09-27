#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""1) 冲突行到底在哪些层存在（05 / 09 / 10_Reports）2) 冲突按亚型分类计数"""
import os, csv, glob
from collections import Counter

ROOT = os.path.expanduser("~/MMPV-paper")
BASE = os.path.join(ROOT, "goji-virome/02_novel_virus/RNA-Lycium_barbarum_out")
VMR = os.path.expanduser("~/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")
CONTIGS = ["CRR527041_clean_NODE_1593_length_1853_cov_4.065169",
           "SRR23215176_clean_NODE_605_length_3594_cov_9.685885",
           "CRR527041_clean_NODE_1995_length_1716_cov_23.029215"]

def rows(p):
    with open(p, newline="", encoding="utf-8", errors="surrogateescape") as f:
        return list(csv.reader(f, delimiter="\t"))

print("=== 1) 这 3 条 contig 在 4 个文件里的样子 ===")
files = [
    ("05_Taxonomy/Votus.integrated/final_integrated_classification.tsv", "05 整合表(20列)"),
    ("09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv", "09 All_plant(23列)"),
    ("10_Reports/final_integrated_classification.tsv", "10_Reports 整合表"),
    ("10_Reports/All_plant.viruses_info.tsv", "10_Reports All_plant"),
]
for rel, label in files:
    p = os.path.join(BASE, rel)
    if not os.path.exists(p):
        print("  %-22s 不存在" % label); continue
    r = rows(p); idx = {h: i for i, h in enumerate(r[0])}
    for cid in CONTIGS:
        hit = [x for x in r[1:] if x and x[0] == cid]
        if hit:
            x = hit[0]
            print("  %-22s %-46s Family=%-16s Genus=%-12s Species=%s" % (
                label, cid[:46], x[idx["Family"]], x[idx["Genus"]], x[idx["Species"]]))
        else:
            print("  %-22s %-46s (无)" % (label, cid[:46]))

# ---- 亚型分类 ----
sp2fam, fam2realm = {}, {}
with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t"); next(rd)
    for r in rd:
        if len(r) < 27: continue
        sp, fam = r[17].strip().strip('"'), r[13].strip().strip('"')
        if sp: sp2fam.setdefault(sp, fam)
        if fam: fam2realm.setdefault(fam, r[3].strip().strip('"'))

DNA_FAM = set()
with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t"); next(rd)
    for r in rd:
        if len(r) < 27: continue
        fam, g = r[13].strip().strip('"'), r[25].strip().strip('"')
        if fam and "DNA" in g.upper(): DNA_FAM.add(fam)

sub = Counter(); rows_dna_fam = 0
for pat in ("*/*/*_out/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv",
            "*/onekp-virus/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv"):
    for q in sorted(glob.glob(os.path.join(ROOT, pat))):
        r = rows(q); idx = {h: i for i, h in enumerate(r[0])}
        for x in r[1:]:
            if not x or len(x) < 13: continue
            sp, fam = x[idx["Species"]].strip(), x[idx["Family"]].strip()
            if sp not in sp2fam or not fam or fam == "NA": continue
            vf = sp2fam[sp]
            if fam == vf: continue
            if fam in DNA_FAM:
                sub["A 高等级阶元指向 DNA 病毒科，属种是 RNA 病毒"] += 1
            elif "viridae" in fam and "viridae" in vf:
                sub["B 同为病毒科但分属不同类群(RT 类为主)"] += 1
            else:
                sub["C 其他"] += 1
            rm, kd = x[idx["Realm"]].strip(), x[idx["Kingdom"]].strip()
            if rm and kd and rm not in ("NA",) and kd not in ("NA",):
                if (rm == "Riboviria" and kd != "Orthornavirae" and kd != "Pararnavirae") or \
                   (rm in ("Varidnaviria", "Monodnaviria") and kd in ("Orthornavirae", "Pararnavirae")):
                    sub["  其中 Realm 与 Kingdom 不同域(彻底的自相矛盾)"] += 1
print("\n=== 2) 625 行冲突的亚型拆解 ===")
for k, v in sub.most_common():
    print("  %-46s %d" % (k, v))
