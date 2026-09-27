#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chk_verdict.py — 判定 Biavirus 那批 contig 双方投票的原始证据
ACVirus side : ACVirus_results/Votus.acvirus
RVDB side    : diamond_output/Votus_diamond_lca_raw.tsv + mmseqs_results/Votus_lca.tsv
"""
import os, sys, csv
from collections import Counter, defaultdict

DS = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
CL = os.path.join(DS, "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")
CLASSED = os.path.join(DS, "05_Taxonomy/Votus.classed")
PLANT = os.path.join(DS, "09_Virome_Analysis/All_plant.viruses.fasta")


def head(path, n=2):
    print("=" * 120)
    print("FILE:", path)
    if not os.path.exists(path):
        print("  [缺失]")
        return []
    print("  size = %.1f MB" % (os.path.getsize(path) / 1e6))
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            if i >= n:
                break
            rows.append(line.rstrip("\n"))
            print("  L%d: %s" % (i + 1, line.rstrip("\n")[:400]))
    return rows


def fa_ids(path):
    ids = set()
    if not os.path.exists(path):
        return ids
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


# ---------- 1. 取目标 contig ----------
targets, plant_ids = [], fa_ids(PLANT)
with open(CL, "r", encoding="utf-8", errors="replace") as f:
    rd = csv.DictReader(f, delimiter="\t")
    HDR = rd.fieldnames
    for r in rd:
        if (r.get("Genus") or "").strip() == "Biavirus":
            targets.append(r["contig_id"])
targets = sorted(set(targets))
print("=" * 120)
print("A. 产出表里 Genus == Biavirus 的 contig：%d 个" % len(targets))
print("   其中落在 All_plant 子集：%d 个" % len([t for t in targets if t in plant_ids]))
for t in targets:
    print("     %-55s plant=%s" % (t, t in plant_ids))

T = set(targets)

# ---------- 2. ACVirus 自己的输出 ----------
print()
head(os.path.join(CLASSED, "ACVirus_results/Votus.acvirus"), 3)

# ---------- 3. diamond_lca raw (RVDB) ----------
raw = os.path.join(CLASSED, "diamond_output/Votus_diamond_lca_raw.tsv")
head(raw, 3)
print()
print("=" * 120)
print("C. RVDB 侧命中明细（diamond_lca raw）")
subj = Counter()
subj_lineage = {}
qn = 0
if os.path.exists(raw):
    with open(raw, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            q = p[0].split()[0]
            if q not in T:
                continue
            s = p[1].split()[0]
            subj[s] += 1
            subj_lineage[s] = p[-1]
            print("  %-52s -> %-20s pident=%-6s len=%-6s ev=%-10s bits=%-8s %s"
                  % (q, s, p[2] if len(p) > 2 else "?", p[3] if len(p) > 3 else "?",
                     p[10] if len(p) > 10 else "?", p[11] if len(p) > 11 else "?",
                     (p[-1] or "")[:90]))
print()
print("  命中过的 subject 汇总（去重 %d 个）:" % len(subj))
for s, c in subj.most_common():
    print("    %-20s x%-4d %s" % (s, c, (subj_lineage.get(s) or "")[:120]))

# ---------- 4. mmseqs lca ----------
lca = os.path.join(CLASSED, "mmseqs_results/Votus_lca.tsv")
head(lca, 2)
print()
print("=" * 120)
print("D. mmseqs lca 里这 13 个 contig 的行")
if os.path.exists(lca):
    with open(lca, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            q = line.split("\t")[0].split()[0] if line.strip() else ""
            if q in T:
                print("  " + line.rstrip("\n")[:400])
