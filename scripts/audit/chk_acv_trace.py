#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chk_acv_trace.py — 追踪 ACVirus 对 13 个 contig 判成 Biavirus 的完整链路
产物目录：05_Taxonomy/Votus.classed/ACVirus_results/Votus.acvirus/
"""
import os, csv
from collections import Counter, defaultdict

B = ("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/"
     "RNA-Lycium_barbarum_out/05_Taxonomy/Votus.classed/ACVirus_results/Votus.acvirus")

T = [
 "CRR1440136_clean_NODE_2764", "CRR1440137_clean_NODE_812",
 "CRR2703946_clean_NODE_2830", "CRR2703948_clean_NODE_611",
 "CRR527041_clean_NODE_304", "CRR527042_clean_NODE_477",
 "CRR527046_clean_NODE_441", "CRR527046_clean_NODE_560",
 "SRR22742702_clean_NODE_373", "SRR23107136_clean_NODE_15609",
 "SRR23215179_clean_NODE_2046", "SRR24305226_clean_NODE_2281",
 "SRR33536650_clean_NODE_3026",
]

FILES = ["diamond_result.tsv", "diamond_pre_deal.tsv", "network.tsv",
         "result_pre.txt", "medium_result.csv", "final_result.tsv",
         "final_result_with_confidence.tsv", "final_node.csv", "final_network.csv"]

print("=" * 120)
print("A. 各产物表头")
for fn in FILES:
    p = os.path.join(B, fn)
    if not os.path.exists(p):
        print("  [缺失] %s" % fn); continue
    with open(p, "r", encoding="utf-8", errors="replace") as f:
        lines = [next(f, "") for _ in range(2)]
    print("\n--- %s  (%.1f MB) ---" % (fn, os.path.getsize(p) / 1e6))
    for i, l in enumerate(lines):
        print("   L%d: %s" % (i + 1, l.rstrip()[:500]))

def scan(fn, ncol_key=0, label=""):
    p = os.path.join(B, fn)
    if not os.path.exists(p): return
    print("\n" + "=" * 120)
    print("B. %s — %s" % (fn, label))
    per = defaultdict(list)
    with open(p, "r", encoding="utf-8", errors="replace") as f:
        rd = csv.reader(f, delimiter="\t" if fn.endswith(("tsv", "txt")) else ",")
        for row in rd:
            if not row: continue
            q = row[0].split()[0]
            base = q
            for t in T:
                if q.startswith(t):
                    base = t; break
            if base in T:
                per[base].append(row)
    for t in T:
        rows = per.get(t, [])
        print("\n  [%s] %d 行" % (t, len(rows)))
        for r in rows[:14]:
            print("      " + " | ".join(x[:70] for x in r[:14]))
        if len(rows) > 14:
            print("      ... 余 %d 行" % (len(rows) - 14))

scan("diamond_result.tsv", label="ORF 级 DIAMOND 命中")
scan("final_result_with_confidence.tsv", label="最终带置信度判定")
scan("result_pre.txt", label="聚合后的查询命中表")
