#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Biavirus 裁决第三层：数据库覆盖对照 + ACVirus 自身证据链。

要回答的核心问题：465 行的 Genus=Biavirus 是「独立证据支持」还是「ACVirus 单工具
把一个参考基因组当兜底」。RVDB v31 不含 Biavirus 蛋白，所以必须看 ACVirus 自己的库。
"""
import csv
import os
from collections import Counter, defaultdict

DMP = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
P2T = "/home/zhangwenda/database/virus-db/RVDB-v31/protid2taxid.map"
ACDB = "/home/zhangwenda/database/virus-db/acvirus_db"
BK = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy/Votus.classed/ACVirus_results/Votus.acvirus"


def p(*a):
    print(*a, flush=True)


# dmp
dmp = {}
with open(DMP, encoding="utf-8", errors="replace") as fh:
    for line in fh:
        if not line or line[0] == "#":
            continue
        f = [x.strip().strip("|").strip() for x in line.rstrip("\n").split("\t|\t")]
        if len(f) >= 10:
            dmp[f[0]] = tuple(f[1:10])

# RVDB 蛋白计数
cnt = Counter()
with open(P2T, encoding="utf-8", errors="replace") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            cnt[f[1]] += 1
p("RVDB v31 map: 蛋白条目 %d，涉及 taxid %d" % (sum(cnt.values()), len(cnt)))

p("\n=== A. family == Schizomimiviridae 的全部 taxid（name/species/genus/taxid/RVDB蛋白数）===")
tot = 0
for tx, r in dmp.items():
    if r[3] == "Schizomimiviridae":
        p("   %-38s sp=%-34s gen=%-14s taxid=%-9s nprot=%d"
          % (r[0], r[1] or "-", r[2] or "-", tx, cnt.get(tx, 0)))
        tot += cnt.get(tx, 0)
p("   该科在 RVDB v31 的蛋白总数 = %d" % tot)

p("\n=== B. genus=Biavirus / Kratosvirus 的 taxid 及其 RVDB 蛋白数 ===")
for g in ("Biavirus", "Kratosvirus"):
    for tx, r in dmp.items():
        if r[2] == g:
            p("   %-20s %-38s fam=%-22s taxid=%-9s nprot=%d"
              % (g, r[0], r[3] or "-", tx, cnt.get(tx, 0)))

p("\n=== C. 我们 hit 里出现最多的 taxid 及其归属（top20）===")
hc = Counter()
for path in ("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated/calibration_20260914/biavirus_arbitration/hits_cat.tsv",):
    a2t = {}
    with open(P2T, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                a2t[f[0]] = f[1]
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1] in a2t:
                hc[a2t[f[1]]] += 1
for tx, n in hc.most_common(20):
    r = dmp.get(tx)
    if r:
        p("   %6d  %-40s gen=%-16s fam=%-22s" % (n, r[0], r[2] or "-", r[3] or "-"))
    else:
        p("   %6d  taxid=%s (dmp 无记录)" % (n, tx))

p("\n=== D. Dishui Lake large algae virus 1 ===")
for tx, r in dmp.items():
    if "Dishui" in r[0]:
        p("   %-40s sp=%-30s gen=%-14s fam=%-22s taxid=%-9s nprot=%d"
          % (r[0], r[1] or "-", r[2] or "-", r[3] or "-", tx, cnt.get(tx, 0)))

p("\n=== E. ACVirus 库里有没有 HG999358 / Aureococcus / Kratosvirus ===")
for fn in ("vmr.tsv", "fixed_vmr_b.tsv", "species.tsv", "taxa.txt",
           "processed_accessions_b.tsv", "processed_accessions_b.fa_names.tsv",
           "taxon_min_coverage.csv"):
    fp = os.path.join(ACDB, fn)
    if not os.path.exists(fp):
        continue
    hits = {"HG999358": 0, "Aureococcus": 0, "Kratosvirus": 0, "Biavirus": 0,
            "Prymnesium": 0, "quantuckense": 0}
    with open(fp, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            for k in hits:
                if k in line:
                    hits[k] += 1
    p("   %-42s %s" % (fn, {k: v for k, v in hits.items() if v}))

p("\n=== F. ACVirus 中间产物 ===")
for fn in sorted(os.listdir(BK)):
    fp = os.path.join(BK, fn)
    if os.path.isfile(fp):
        p("   %-28s %10d B" % (fn, os.path.getsize(fp)))
for fn in ("diamond_result.tsv", "diamond_pre_deal.tsv", "result_pre.txt",
           "network.tsv", "medium_result.csv"):
    fp = os.path.join(BK, fn)
    if os.path.exists(fp):
        p("\n   --- head %s ---" % fn)
        with open(fp, encoding="utf-8", errors="replace") as fh:
            for i, line in enumerate(fh):
                if i >= 3:
                    break
                p("   " + line.rstrip()[:300])

p("\n=== G. taxon_min_coverage.csv 里 Schizomimiviridae/Kratosvirus/Biavirus 行 ===")
fp = os.path.join(ACDB, "taxon_min_coverage.csv")
if os.path.exists(fp):
    with open(fp, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if any(k in line for k in ("Schizomimiviridae", "Kratosvirus", "Biavirus")):
                p("   " + line.rstrip()[:200])
