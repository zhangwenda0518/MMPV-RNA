#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Mimiviridae 这类科能不能拿来判"植物"？先看 VMR 自己怎么标宿主，再拿 625 行冲突行验证两侧的宿主标注。"""
import os, csv, glob
from collections import Counter, defaultdict

VMR = os.path.expanduser("~/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")
ROOT = os.path.expanduser("~/MMPV-paper")

fam_host = defaultdict(Counter)
fam_genome = defaultdict(Counter)
sp_fam, sp_host, sp_genome, sp_name = {}, {}, {}, {}
gen_fam = {}

with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t")
    hdr = next(rd)
    for r in rd:
        if len(r) < 27:
            continue
        fam = r[13].strip().strip('"'); gen = r[15].strip().strip('"')
        sp = r[17].strip().strip('"'); host = r[26].strip().strip('"')
        gnm = r[25].strip().strip('"')
        if fam:
            fam_host[fam][host] += 1
            fam_genome[fam][gnm] += 1
        if sp:
            sp_fam.setdefault(sp, fam); sp_host.setdefault(sp, host); sp_genome.setdefault(sp, gnm)
        if gen:
            gen_fam.setdefault(gen, fam)

print("=== 1) VMR 里关键科自己的宿主标注（Host source 列原样） ===")
for fam in ["Mimiviridae", "Schizomimiviridae", "Phycodnaviridae", "Nimaviridae", "Potyviridae",
            "Caulimoviridae", "Pseudoviridae", "Metaviridae", "Tymoviridae", "Secoviridae",
            "Iridoviridae", "Poxviridae", "Marseilleviridae", "Pithoviridae", "Alphaflexiviridae"]:
    if fam in fam_host:
        print("  %-20s n=%-5d host=%s | genome=%s" % (
            fam, sum(fam_host[fam].values()),
            dict(fam_host[fam].most_common(3)), dict(fam_genome[fam].most_common(2))))
    else:
        print("  %-20s (VMR 无此科)" % fam)

print("\n=== 2) VMR 全库：科的宿主标注 top 25（看 Host source 列的取值体系） ===")
for fam, c in sorted(fam_host.items(), key=lambda kv: -sum(kv[1].values()))[:25]:
    print("  %-24s n=%-6d %s" % (fam, sum(c.values()), dict(c.most_common(3))))

print("\n=== 2b) Host source 列全部取值形态 top 30 ===")
allh = Counter()
for c in fam_host.values():
    allh.update(c)
for h, n in allh.most_common(30):
    print("  %-40s %d" % (repr(h), n))


def rows(p):
    with open(p, newline="", encoding="utf-8", errors="surrogateescape") as f:
        return list(csv.reader(f, delimiter="\t"))


def is_plant(h):
    h = (h or "").lower()
    return any(k in h for k in ("plant", "host", "angiosperm", "dicot", "monocot"))


conf = 0
fam_side, sp_side = Counter(), Counter()
pair = Counter()
fam_plant = sp_plant = both_plant = neither = 0
for pat in ("*/*/*_out/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv",
            "*/onekp-virus/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv"):
    for q in sorted(glob.glob(os.path.join(ROOT, pat))):
        r = rows(q)
        idx = {h: i for i, h in enumerate(r[0])}
        for x in r[1:]:
            if not x or len(x) < 13:
                continue
            sp = x[idx["Species"]].strip(); fam = x[idx["Family"]].strip()
            if sp not in sp_fam or not fam or fam == "NA":
                continue
            if fam == sp_fam[sp]:
                continue
            conf += 1
            fam_side[fam] += 1
            sp_side[sp_fam[sp]] += 1
            pair[(fam, sp_fam[sp])] += 1
            fh = fam_host.get(fam, Counter()).most_common(1)
            sh = fam_host.get(sp_fam[sp], Counter()).most_common(1)
            fh = fh[0][0] if fh else ""
            sh = sh[0][0] if sh else ""
            fp, spt = is_plant(fh), is_plant(sh)
            fam_plant += fp; sp_plant += spt
            if fp and spt: both_plant += 1
            if not fp and not spt: neither += 1

print("\n=== 3) 625 行冲突行：科侧与种侧的 VMR 宿主标注 ===")
print("  冲突行数 %d" % conf)
print("  [科侧] top12")
for fam, n in fam_side.most_common(12):
    hh = dict(fam_host.get(fam, Counter()).most_common(2))
    print("    %-24s n=%-4d host=%s" % (fam, n, hh))
print("  [种侧/VMR 该种所属科] top12")
for fam, n in sp_side.most_common(12):
    hh = dict(fam_host.get(fam, Counter()).most_common(2))
    print("    %-24s n=%-4d host=%s" % (fam, n, hh))
print("\n  科侧首要宿主标为植物(含 plant/host 字样)的行: %d" % fam_plant)
print("  种侧首要宿主标为植物的行:                 %d" % sp_plant)
print("  两侧都标植物: %d | 两侧都非植物: %d" % (both_plant, neither))
print("\n  [冲突组合 top 15] 科 -> 该种真实科")
for (a, b), n in pair.most_common(15):
    print("    %-24s -> %-24s %d" % (a, b, n))
