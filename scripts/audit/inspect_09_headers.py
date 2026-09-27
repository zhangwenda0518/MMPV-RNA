#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""09 层候选表结构对比（只读）：找出可 join VMR 的键列，以及是否已有 DNA/RNA 类字段。"""
import csv, os, glob

ROOT = os.path.expanduser("~/MMPV-paper")
PATTERNS = [
    ("*_out", "09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv"),
    ("*_out", "09_Virome_Analysis/HQ_analysis/HQ_plant_viruses_info.tsv"),
    ("*_out", "09_Virome_Analysis/all_plant_analysis/all_plant_viruses_genus_summary.tsv"),
    ("*_out", "09_Virome_Analysis/viroid_analysis/Viroid.all_info.tsv"),
    ("*_out", "09_Virome_Analysis/rescue_detection/summary/all_viruses.summary.tsv"),
    ("*_out", "09_Virome_Analysis/rescue_detection/rescue_ref_info.tsv"),
]

def head(p, n=1):
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        rd = csv.reader(f, delimiter="\t")
        try:
            hdr = [h.strip().strip('"') for h in next(rd)]
        except StopIteration:
            return None, []
        rows = [r for _, r in zip(range(n), rd)]
    return hdr, rows

for base, rel in PATTERNS:
    hits = sorted(glob.glob(os.path.join(ROOT, "*", "*", base, rel))) + \
           sorted(glob.glob(os.path.join(ROOT, "*", base, rel)))
    if not hits:
        print("### %s  -> 无命中" % rel)
        continue
    print("### %s  -> %d 个" % (rel, len(hits)))
    hdr, rows = head(hits[0], 1)
    if hdr is None:
        print("    (空文件)")
        continue
    print("    列数=%d" % len(hdr))
    print("    表头: %s" % " | ".join(hdr))
    dn = [h for h in hdr if any(k in h.lower() for k in ("genome", "dna", "rna", "type", "struc"))]
    print("    可能已含 DNA/RNA 类字段: %s" % (dn if dn else "无"))
    if rows:
        print("    样例: %s" % " | ".join(str(x)[:28] for x in rows[0]))
    # 检查各项目表头是否一致
    sigs = {}
    for h in hits:
        hh, _ = head(h, 0)
        sigs.setdefault(tuple(hh) if hh else None, []).append(h)
    print("    表头一致性: %d 种不同表头" % len(sigs))
    for k, v in sigs.items():
        if k is None:
            print("      (空) x%d" % len(v))
        else:
            print("      x%d 列数=%d : %s" % (len(v), len(k), " | ".join(k[:8])))
    print()
