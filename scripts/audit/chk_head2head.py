#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chk_head2head.py — 取 mmseqs 最佳比对数值（RVDB 侧），与 ACVirus 侧做同口径对比
"""
import os
from collections import defaultdict

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


def base_of(q):
    q = q.split()[0]
    for t in T:
        if q.startswith(t):
            return t
    return None


aln = os.path.join(B, "mmseqs_results/Votus_tophit_aln")
print("=" * 130)
print("A. mmseqs tophit_aln（RVDB 侧最佳比对）：query | target | fident | alnLen | qlen* | evalue | bits")
print("=" * 130)
rows = defaultdict(list)
with open(aln, "r", encoding="utf-8", errors="replace") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        b = base_of(p[0])
        if b:
            rows[b].append(p)

for t in T:
    rs = rows.get(t, [])
    if not rs:
        print("  [%-38s] 无比对" % t[:38]); continue
    for r in rs:
        print("  %-38s -> %-14s fid=%-7s alnLen=%-6s qspan=%s-%s/%s ev=%-11s bits=%s"
              % (t[:38], r[1], r[2], r[3], r[6], r[7], r[5], r[10], r[11]))

# 对照：同一文件里所有 contig 的 fident 分布 + 身份最高的前 15 个
print("\n" + "=" * 130)
print("B. 全数据集 mmseqs 最佳比对 fident 分布（判断这些 Biavirus contig 是否系统性偏低）")
allf = []
with open(aln, "r", encoding="utf-8", errors="replace") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        try:
            allf.append((float(p[2]), p[0], p[1], p[3], p[10], p[11]))
        except Exception:
            pass
allf.sort(reverse=True)
print("  总比对 %d 条" % len(allf))
n = len(allf)
for q in [0, 25, 50, 75, 90, 99, 100]:
    pass
vals = sorted(x[0] for x in allf)
def pct(p):
    return vals[min(n - 1, int(n * p / 100))]
print("  fident 分位：p0=%.3f p10=%.3f p25=%.3f p50=%.3f p75=%.3f p90=%.3f p100=%.3f"
      % (pct(0), pct(10), pct(25), pct(50), pct(75), pct(90), pct(100)))
print("\n  这 13 个 contig 的 fident 在全数据集里的排名（越高越前面）")
idx = {x[1]: i for i, x in enumerate(allf)}
for pf, q, tg, al, ev, bt in allf:
    b = base_of(q)
    if b:
        print("     %-38s fident=%.3f rank=%d/%d  target=%s ev=%s" % (b[:38], pf, idx[q] + 1, n, tg, ev))
        T.remove(b)

print("\n" + "=" * 130)
print("C. fident 最高的前 15 条（头部是什么）")
for pf, q, tg, al, ev, bt in allf[:15]:
    print("     %-46s fident=%.3f alnLen=%-5s %-14s ev=%s" % (q[:46], pf, al, tg, ev))
