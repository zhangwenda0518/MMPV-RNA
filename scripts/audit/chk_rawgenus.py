#!/usr/bin/env python3
"""根因核查：CAT/diamond_lca/metabuli 是库没属级注释，还是解析吃了属"""
import csv
from collections import Counter, defaultdict
from pathlib import Path

P = Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/"
         "05_Taxonomy/Votus.classed/Votus_combined_taxonomy.tsv")
TGT = "CRR1126135_clean_NODE_68_length_2879"
RANKS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
BAD = {"", "-", "NA", "na", "N/A", "no rank", "undefined", "unknown", "null", "default",
       "Unclassified"}

with open(P, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.DictReader(f, delimiter="\t")
    print("列名:", rd.fieldnames)
    rows = list(rd)
print("总行数:", len(rows))

idcol = "seq_name" if "seq_name" in (rows[0] if rows else {}) else "contig_id"
print("用的 id 列:", idcol)

print()
print("=== 目标 contig 在 combined 表里的原始记录 ===")
hit = False
for r in rows:
    if (r.get(idcol) or "").startswith(TGT):
        hit = True
        print("  tool=%s" % r.get("tool"))
        for k in RANKS:
            if k in r:
                print("      %-14s %r" % (k, r[k]))
if not hit:
    print("  未命中, 试前缀匹配")
    for r in rows[:2]:
        print("  样例:", {k: r[k] for k in list(r)[:12]})

print()
print("=== 各工具各层级的非空率（看属级是不是整体缺） ===")
cnt = defaultdict(Counter)
tot = Counter()
present = set()
for r in rows:
    t = r.get("tool") or "?"
    present.add(t)
    tot[t] += 1
    for k in RANKS:
        v = (r.get(k) or "").strip()
        if v not in BAD:
            cnt[t][k] += 1
RANKS2 = [k for k in RANKS if k in (rows[0] if rows else {})] or RANKS
print("  %-14s %6s %s" % ("tool", "rows", " ".join("%9s" % k[:9] for k in RANKS2)))
for t in sorted(tot):
    print("  %-14s %6d %s" % (t, tot[t],
          " ".join("%9d" % cnt[t][k] for k in RANKS2)))
