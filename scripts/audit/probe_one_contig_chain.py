#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""单条 contig 的完整证据链：8 个工具各判了什么 + 最终整合表怎么拼出来的"""
import os, csv, glob

BASE = os.path.expanduser("~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out")
V = os.path.join(BASE, "05_Taxonomy/Votus.integrated")
CONTIGS = ["CRR527041_clean_NODE_1593_length_1853_cov_4.065169",
           "SRR23215176_clean_NODE_605_length_3594_cov_9.685885",
           "CRR527041_clean_NODE_1995_length_1716_cov_23.029215"]

def rows(p, delim="\t"):
    with open(p, newline="", encoding="utf-8", errors="surrogateescape") as f:
        return list(csv.reader(f, delimiter=delim))

for cid in CONTIGS:
    print("=" * 130)
    print("contig: %s" % cid)
    for label, fn in [("最终整合表", "final_integrated_classification.tsv")]:
        p = os.path.join(V, fn)
        r = rows(p); hdr = r[0]
        for x in r[1:]:
            if x and x[0] == cid:
                print("  --- %s ---" % label)
                for h, v in zip(hdr, x):
                    print("      %-28s = %s" % (h, v))
    print("  --- 8 个工具各自的判定 ---")
    for p in sorted(glob.glob(os.path.join(V, "standardized_*.tsv"))):
        r = rows(p)
        hit = [x for x in r[1:] if x and x[0] == cid]
        if hit:
            x = hit[0]
            print("      %-16s Realm=%s | Kingdom=%s | Family=%s | Genus=%s | Species=%s" % (
                os.path.basename(p).replace("standardized_", "").replace(".tsv", ""),
                x[1], x[2], x[6], x[7], x[8]))
        else:
            print("      %-16s (无记录)" % os.path.basename(p))
    print("  --- comparison_* 表里的取值分布 ---")
    for p in sorted(glob.glob(os.path.join(V, "comparison_*.tsv"))):
        r = rows(p)
        hit = [x for x in r[1:] if x and x[0] == cid]
        if hit:
            x = hit[0]
            d = dict(zip(r[0], x))
            print("      %-22s all_vals=%s | consensus=%s | is_consistent=%s" % (
                os.path.basename(p).replace("comparison_", "").replace(".tsv", ""),
                d.get("all_vals_str", "")[:60], d.get("consensus_taxon", ""), d.get("is_consistent", "")))
