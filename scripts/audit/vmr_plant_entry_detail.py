#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 VMR 宿主列含 plant 的争议条目逐个摊开: 谁、什么种、宿主列原文、所在科属。
只读。用法: python3 /tmp/vmr_plant_entry_detail.py
"""
import csv
import re

VMR = "/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv"
GEN_INTEREST = ["Miraophiovirus", "Oleurovirus", "Alphambiguivirus", "Betambiguivirus",
                "Rimosavirus", "Sylvanvirus", "Biavirus"]
FAM_INTEREST = ["Ouroboviridae", "Discoviridae", "Ambiguiviridae", "Partitiviridae",
                "Caulimoviridae", "Steitzviridae", "Fiersviridae", "Circoviridae",
                "Draupnirviridae"]
PLANT_PAT = re.compile(r"plant|viridiplant", re.I)

rows = list(csv.reader(open(VMR, encoding="utf-8", errors="replace"), delimiter="\t"))
hdr = rows[0]
print("列: ", {i: h for i, h in enumerate(hdr) if i in (0, 1, 13, 15, 17, 20, 26)})
ID, SPEC, FAM, GEN, HOST = 0, 1, 13, 15, 26

print()
print("=" * 90)
print("A. 属级争议条目: 含 plant 的具体种与宿主列原文")
print("=" * 90)
seen_gen = {}
for r in rows[1:]:
    if len(r) <= HOST:
        continue
    g = r[GEN].strip()
    if g not in GEN_INTEREST:
        continue
    h = r[HOST].strip()
    seen_gen.setdefault(g, []).append((r[FAM].strip(), r[SPEC].strip(), r[20].strip(), h))
for g in GEN_INTEREST:
    lst = seen_gen.get(g)
    print()
    print("  [%s]  VMR 记录 %d 条" % (g, len(lst) if lst else 0))
    if not lst:
        continue
    fams = sorted({x[0] for x in lst})
    hosts = {}
    for _, _, _, h in lst:
        hosts[h] = hosts.get(h, 0) + 1
    print("     所在科: %s" % ", ".join(fams))
    print("     宿主列分布: %s" % ", ".join("%s(%d)" % (k, v) for k, v in
                                          sorted(hosts.items(), key=lambda x: -x[1])))
    for fam, sp, vname, h in lst:
        if PLANT_PAT.search(h):
            print("     * %-34s %-30s host=%s" % (fam, sp[:34], h))

print()
print("=" * 90)
print("B. 科级争议条目: 含 plant 的具体种与宿主列原文")
print("=" * 90)
seen_fam = {}
for r in rows[1:]:
    if len(r) <= HOST:
        continue
    f = r[FAM].strip()
    if f not in FAM_INTEREST:
        continue
    seen_fam.setdefault(f, []).append((r[GEN].strip(), r[SPEC].strip(), r[HOST].strip()))
for f in FAM_INTEREST:
    lst = seen_fam.get(f)
    print()
    print("  [%s]  VMR 记录 %d 条" % (f, len(lst) if lst else 0))
    if not lst:
        continue
    hosts = {}
    for _, _, h in lst:
        hosts[h] = hosts.get(h, 0) + 1
    print("     宿主列分布: %s" % ", ".join("%s(%d)" % (k, v) for k, v in
                                          sorted(hosts.items(), key=lambda x: -x[1])))
    for gen, sp, h in lst:
        if PLANT_PAT.search(h) and "(S)" not in h:
            print("     * %-22s %-38s host=%s" % (gen[:22], sp[:38], h))
print()
print("DONE")
