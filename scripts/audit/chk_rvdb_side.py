#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chk_rvdb_side.py — CAT / diamond_lca / mmseqs 三家（同用 RVDB-v31）对这 13 个 contig 的原始命中
"""
import os
from collections import Counter, defaultdict

B = ("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/"
     "RNA-Lycium_barbarum_out/05_Taxonomy/Votus.classed")

T = [
 "CRR1440136_clean_NODE_2764", "CRR1440137_clean_NODE_812",
 "CRR2703946_clean_NODE_2830", "CRR2703948_clean_NODE_611",
 "CRR527041_clean_NODE_304", "CRR527042_clean_NODE_477",
 "CRR527046_clean_NODE_441", "CRR527046_clean_NODE_560",
 "SRR22742702_clean_NODE_373", "SRR23107136_clean_NODE_15609",
 "SRR23215179_clean_NODE_2046", "SRR24305226_clean_NODE_2281",
 "SRR33536650_clean_NODE_3026",
]
S = set(T)


def base_of(q):
    q = q.split()[0]
    for t in T:
        if q.startswith(t):
            return t
    return None


# ---------- 1. diamond_lca raw ----------
raw = os.path.join(B, "diamond_output/Votus_diamond_lca_raw.tsv")
print("=" * 120)
print("A. diamond_lca raw: %s" % raw)
print("   size = %.1f MB" % (os.path.getsize(raw) / 1e6))
with open(raw, "r", encoding="utf-8", errors="replace") as f:
    n = 0
    for line in f:
        if n >= 3:
            break
        print("   L%d: %s" % (n + 1, line.rstrip()[:600]))
        n += 1

print("\n" + "=" * 120)
print("B. 这 13 个 contig 在 RVDB（diamond_lca raw）里的全部命中")
per = defaultdict(list)
subj_all = Counter()
lineage_of = {}
with open(raw, "r", encoding="utf-8", errors="replace") as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        b = base_of(p[0])
        if b is None:
            continue
        per[b].append(p)
        s = p[1].split()[0]
        subj_all[s] += 1
        lineage_of[s] = p[-1]

for t in T:
    rows = per.get(t, [])
    print("\n  [%s]  %d 条命中" % (t, len(rows)))
    for r in rows[:10]:
        print("      " + " | ".join(x[:75] for x in r))

print("\n" + "=" * 120)
print("C. 命中过的 RVDB subject 汇总（去重 %d 个）" % len(subj_all))
for s, c in subj_all.most_common(40):
    print("    %-22s x%-4d %s" % (s, c, (lineage_of.get(s) or "")[:110]))

# ---------- 2. mmseqs lca ----------
lca = os.path.join(B, "mmseqs_results/Votus_lca.tsv")
print("\n" + "=" * 120)
print("D. mmseqs Votus_lca.tsv")
if os.path.exists(lca):
    with open(lca, "r", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            if i < 2:
                print("   L%d: %s" % (i + 1, line.rstrip()[:400]))
            if base_of(line) if line.strip() else False:
                print("   * " + line.rstrip()[:400])

# ---------- 3. mmseqs tophit_aln ----------
aln = os.path.join(B, "mmseqs_results/Votus_tophit_aln")
print("\n" + "=" * 120)
print("E. mmseqs tophit_aln 前 2 行")
if os.path.exists(aln):
    with open(aln, "r", encoding="utf-8", errors="replace") as f:
        for i in range(2):
            l = f.readline()
            if l:
                print("   L%d: %s" % (i + 1, l.rstrip()[:400]))
