#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""落地后验收：1) CRLF 文件是否正常 2) 05/10_Reports 层未被触碰 3) 各项目 DNA/RNA 计数"""
import csv, os, glob, hashlib
from collections import Counter

ROOT = os.path.expanduser("~/MMPV-paper")

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:12]

def rows(p):
    with open(p, newline="", encoding="utf-8", errors="surrogateescape") as f:
        return list(csv.reader(f, delimiter="\t"))

print("=== 1) CRLF 文件检查 ===")
for p in ["onekp-virome/onekp-virus/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv",
          "onekp-virome/onekp-virus/09_Virome_Analysis/all_plant_analysis/all_plant_viruses_genus_summary.tsv",
          "goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv"]:
    q = os.path.join(ROOT, p)
    raw = open(q, "rb").read()
    eol = "CRLF" if b"\r\n" in raw else ("CR" if b"\r" in raw else "LF")
    r = rows(q)
    print("  %-72s eol=%-4s 行=%d 列=%d 错插(\\r\\t)=%d 末3列=%s" % (
        p.split("/")[-1] + " [" + p.split("/")[0] + "]", eol, len(r), len(r[0]),
        raw.count(b"\r\t"), r[0][-3:]))

print("\n=== 2) 其他层是否被触碰 ===")
checks = [
    ("05_Taxonomy/Votus.integrated/final_integrated_classification.tsv", "goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"),
    ("10_Reports/final_integrated_classification.tsv", "goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"),
    ("10_Reports/All_plant.viruses_info.tsv", "goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"),
]
for rel, base in checks:
    q = os.path.join(ROOT, base, rel)
    if os.path.exists(q):
        r = rows(q)
        print("  %-58s 存在 列=%d 首列=%s 含Genome_Type=%s md5=%s" % (
            rel, len(r[0]), r[0][0], "Genome_Type" in r[0], md5(q)))
    else:
        print("  %-58s 不存在" % rel)

print("\n=== 3) 各项目 DNA/RNA 计数（All_plant / HQ） ===")
pats = ["*/*/*_out/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv",
        "*/*/onekp-virus/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv",
        "*/*/*_out/09_Virome_Analysis/HQ_analysis/HQ_plant_viruses_info.tsv",
        "*/*/onekp-virus/09_Virome_Analysis/HQ_analysis/HQ_plant_viruses_info.tsv"]
files = []
for pat in pats:
    files += sorted(glob.glob(os.path.join(ROOT, pat)))
tot = Counter()
for q in files:
    r = rows(q)
    idx = {h: i for i, h in enumerate(r[0])}
    c = Counter(x[idx["Genome_Type"]] or "未定" for x in r[1:] if x and x[0])
    sp = Counter(x[idx["Genome_Source"]] for x in r[1:] if x and x[0])
    tag = "All_plant" if "HQ_" not in q else "HQ"
    proj = q.split("/02_novel_virus/")[0].split("/")[-1] if "02_novel_virus" in q else q.split("/onekp-virome/")[0].split("/")[-1]
    name = q.split("/")[-2] if "HQ" not in q else q.split("/")[-2]
    print("  %-10s %-16s %-24s DNA=%-5d RNA=%-6d 未定=%-4d | species=%d genus=%d tax=%d" % (
        tag, proj, q.split("02_novel_virus/")[-1].split("/")[0] if "02_novel_virus" in q else "onekp-virus",
        c["DNA"], c["RNA"], c["未定"], sp["species"], sp["genus"], sp["tax"]))
    tot.update(c)
print("  合计 DNA=%d RNA=%d 未定=%d" % (tot["DNA"], tot["RNA"], tot["未定"]))
