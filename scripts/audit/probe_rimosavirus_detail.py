#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""只读：Rimosavirus / Sylvanvirus 相关行的完整展开（属级证据收口）"""
import csv
VMR = "/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv"
TARGETS = ["Rimosavirus"]
with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t")
    hdr = next(rd)
    rows = [r for r in rd if len(r) >= 27]

for t in TARGETS:
    print("===== %s 全部行 =====" % t)
    for r in rows:
        if r[15].strip().strip('"') == t:
            print("  species=%-34s | family=%-14s | genome=%-10s | host_source=%-20s | exemplar=%s | isolate=%s" % (
                r[17].strip().strip('"'), r[13].strip().strip('"'), r[25].strip().strip('"'),
                r[26].strip().strip('"'), r[19].strip().strip('"')[:30], r[22].strip().strip('"')[:30]))

print("\n===== VMR 里 属名 或 种名 含 Rimosavirus 的物种汇总 =====")
sp = {}
for r in rows:
    if "Rimosavirus" in r[17] or "Rimosavirus" in r[15]:
        sp[(r[15].strip().strip('"'), r[17].strip().strip('"'))] = r[26].strip().strip('"')
for k, v in sorted(sp.items()):
    print("  genus=%-14s species=%-34s host_source=%s" % (k[0], k[1], v))

print("\n===== 对照：Tombusviridae 各属的 Host source 形态（该科是否本就是植物病毒科） =====")
fam = {}
for r in rows:
    if r[13].strip().strip('"') == "Tombusviridae":
        fam.setdefault(r[15].strip().strip('"'), set()).add(r[26].strip().strip('"'))
print("  Tombusviridae 属数 = %d" % len(fam))
for g, v in sorted(fam.items()):
    print("   %-18s %s" % (g, sorted(v)))
