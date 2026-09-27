#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""只读：Discoviridae / Ouroboviridae 两科的全部 VMR 记录展开（定"有没有真植物病毒"）"""
import csv
from collections import Counter
VMR = "/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv"
FAMS = ["Discoviridae", "Ouroboviridae"]
with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t")
    next(rd)
    rows = [r for r in rd if len(r) >= 27]

for fam in FAMS:
    sub = [r for r in rows if r[13].strip().strip('"') == fam]
    print("===== %s : %d 行 =====" % (fam, len(sub)))
    print("  Genome:", dict(Counter(r[25].strip().strip('"') for r in sub)))
    print("  Host source:", dict(Counter(r[26].strip().strip('"') for r in sub)))
    for r in sub:
        print("   genus=%-18s species=%-46s host=%-18s genome=%s" % (
            r[15].strip().strip('"'), r[17].strip().strip('"'),
            r[26].strip().strip('"'), r[25].strip().strip('"')))
    print()

print("===== 这三条'植物记录'的宿主展开：查 Virus name(s) 与 accession =====")
KEYS = ["rice dwarf", "Trichosanthes", "Forsythia"]
for r in rows:
    vn = r[20].strip().strip('"')
    if any(k.lower() in vn.lower() for k in KEYS):
        print("  family=%-15s genus=%-18s species=%-46s" % (r[13].strip().strip('"'), r[15].strip().strip('"'), r[17].strip().strip('"')))
        print("      virus_name=%s" % vn[:120])
        print("      host_source=%-18s genome=%-14s genbank=%s" % (r[26].strip().strip('"'), r[25].strip().strip('"'), r[23].strip().strip('"')[:60]))
