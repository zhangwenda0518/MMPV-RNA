#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quantify concentration of R1-only contradictions by concrete label (new products only).

Usage: python3 freq_r1only.py
"""
import csv
import sys
from collections import Counter, defaultdict

LEVELS = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm"]
DMP_COL = {"Realm": "superkingdom", "Kingdom": "kingdom", "Phylum": "phylum",
           "Class": "class", "Order": "order", "Family": "family",
           "Genus": "genus", "Species": "species"}
COLS = ["tax_id", "tax_name", "species", "genus", "family", "order",
        "class", "phylum", "kingdom", "superkingdom"]
PAIRS = [("Realm", "Kingdom"), ("Kingdom", "Phylum"), ("Phylum", "Class"),
         ("Class", "Order"), ("Order", "Family"), ("Family", "Genus"), ("Genus", "Species")]


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None
    return v


DMP = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
PRODUCTS = [
    ("Alternaria", "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Alternaria_alternata_out/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"),
    ("amarum", "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_amarum_out/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"),
    ("Aphis", "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Aphis_gossypii_out/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"),
    ("barbarum", "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"),
    ("chinense", "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_chinense_out/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"),
    ("Fusarium", "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Fusarium_nematophilum_out/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"),
    ("onekp", "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"),
    ("ruthenicum", "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_ruthenicum_out/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"),
]

r1 = {}
r2 = defaultdict(set)
with open(DMP, encoding="utf-8", errors="ignore") as f:
    for line in f:
        parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
        if len(parts) < 10:
            continue
        rec = dict(zip(COLS, parts))
        nm = rec["tax_name"].lower()
        if nm and nm not in r1:
            r1[nm] = {L: rec[DMP_COL[L]] for L in LEVELS}
        for i in range(len(LEVELS) - 1):
            low, high = LEVELS[i], LEVELS[i + 1]
            lv, hv = rec[DMP_COL[low]], rec[DMP_COL[high]]
            if lv and hv:
                r2[(low, lv.lower())].add(hv)

for lab, path in PRODUCTS:
    cnt = Counter()
    tot_only = 0
    with open(path, encoding="utf-8", errors="ignore", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t", quotechar='"')
        for raw in rd:
            row = {k.strip(): norm(v) for k, v in raw.items()}
            for high, low in PAIRS:
                hv, lv = row.get(high), row.get(low)
                if not hv or not lv:
                    continue
                lvv = lv.lower()
                e = r1.get(lvv)
                pv = (e.get(high) or "").strip() if e is not None else ""
                c1 = bool(pv) and pv != hv
                ps = r2.get((low, lvv))
                c2 = bool(ps) and hv not in ps
                if c1 and not c2:
                    tot_only += 1
                    cnt[("%s-%s" % (high, low), low, lv, high, hv, "ref:" + (pv or "-"))] += 1
    print("==== %s  R1only_pair_occurrences=%d  distinct_labels=%d" % (lab, tot_only, len(cnt)))
    for k, v in cnt.most_common(12):
        print("   %5d  %s  %s=[%s]  %s=[%s]  %s" % (v, k[0], k[1], k[2], k[3], k[4], k[5]))
    print()
