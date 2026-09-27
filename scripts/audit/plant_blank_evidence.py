#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""对被置空、且落在植物病毒成果表里的 contig，查 7 工具逐票（Family/Genus）。
判据：置空后残留的 Family 有几票支持？能否支撑发表。
"""
import csv
import os
from collections import Counter, defaultdict

BASE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
RAW = os.path.join(BASE, "05_Taxonomy", "Votus.classed", "Votus_combined_taxonomy.tsv")
PLANT = os.path.join(BASE, "10_Reports", "All_plant.viruses_info.tsv")
CAL = os.path.join(BASE, "05_Taxonomy", "Votus.integrated", "calibration_20260914")


def p(*a):
    print(*a, flush=True)


def rd(path, sep="\t"):
    with open(path, encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter=sep))


act = {}
for d in rd(os.path.join(CAL, "final_integrated_classification.calibrated.tsv")):
    act[(d.get("contig_id") or "").strip().strip('"')] = (d.get("calib_action") or "").strip().strip('"')

plant = rd(PLANT)
hdr = list(plant[0].keys())
idc = "contig_id"
targets = [r for r in plant if act.get((r.get(idc) or "").strip()) == "blank"]
tset = set((r[idc] or "").strip() for r in targets)
p("植物成果表里被置空的 contig = %d" % len(tset))

raw = rd(RAW)
p("原始逐工具表列: %s" % list(raw[0].keys()))
byid = defaultdict(list)
for r in raw:
    k = (r.get("seq_name") or r.get(idc) or "").strip()
    if k in tset:
        byid[k].append(r)
p("命中 %d 个 contig（共 %d 票行）\n" % (len(byid), sum(len(v) for v in byid.values())))

# 找工具名与 Family/Genus 列
tcol = next((c for c in hdr if "tool" in c.lower()), None)
fcol = "Family"
gcol = "Genus"
p("工具列 = %s / 科列 = %s / 属列 = %s" % (tcol, fcol, gcol))

nfam_support = Counter()
detail = []
for cid in sorted(byid):
    votes = byid[cid]
    fams = Counter((v.get(fcol) or "").strip() for v in votes if (v.get(fcol) or "").strip())
    gens = Counter((v.get(gcol) or "").strip() for v in votes if (v.get(gcol) or "").strip())
    integ = next((r for r in targets if (r[idc] or "").strip() == cid), None)
    ifam = (integ.get(fcol) or "").strip() if integ else ""
    sup = fams.get(ifam, 0)
    nfam_support[sup] += 1
    detail.append((cid, ifam, sup, len(fams), fams.most_common(3), gens.most_common(2),
                   (integ.get("primary_tool") or "") if integ else ""))

p("\n=== 残留 Family 的工具支持票数分布 ===")
for k in sorted(nfam_support):
    p("   %d 票支持: %d 行" % (k, nfam_support[k]))

p("\n=== 逐 contig（支持票升序）===")
for cid, ifam, sup, nuniq, fams, gens, pt in sorted(detail, key=lambda x: x[2])[:30]:
    p("   %-52s 表内科=%-20s 支持%2d票 | 科票=%s" % (cid[:52], ifam[:20], sup, fams))
    p("       属票=%s  primary_tool=%s" % (gens, pt))
