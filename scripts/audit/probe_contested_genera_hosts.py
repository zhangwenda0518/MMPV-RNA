#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""只读：6 个争议属名的宿主证据盘点。
输出每属在 VMR / PlantVirusDB 植物库 / ICTV 属表 / 管线参考库 / NCBI 名字表 的命中情况。"""
import csv, os, subprocess

VMR = "/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv"
PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"
ICTVG = "/home/zhangwenda/MMPV-RNA/database/ICTV_MSL41_Genus.tsv"
REFLIB = "/home/zhangwenda/MMPV-RNA/database/final.cluster.ref_info.tsv"
NAMES = "/home/zhangwenda/database/taxonomy/names.dmp"

TARGETS = ["Sylvanvirus", "Crucivirus", "Rimosavirus",
           "Klosneuvirus", "Hokovirus", "Indivirus"]
# 参考：已撤回的两个真植物属
EXTRA = ["Miraophiovirus", "Oleurovirus"]


def head(p, n=1):
    out = []
    if not os.path.isfile(p):
        return ["MISSING " + p]
    with open(p, encoding="utf-8", errors="replace") as f:
        for _ in range(n):
            out.append(f.readline().rstrip("\n"))
    return out


print("=== 参照文件表头 ===")
for p in (PLANT, ICTVG, REFLIB):
    print("---", p)
    for l in head(p):
        cols = l.split("\t")
        print("   列数 %d:" % len(cols), " | ".join("%d=%s" % (i, c.strip('"')) for i, c in enumerate(cols, 1)))

# ---- VMR 建索引 ----
vmr_by_genus = {}
vmr_all_genus = set()
if os.path.isfile(VMR):
    with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = next(rd)
        for row in rd:
            if len(row) < 27:
                continue
            g = row[15].strip().strip('"')
            vmr_all_genus.add(g)
            vmr_by_genus.setdefault(g, []).append(row)
print("\nVMR 属数 =", len(vmr_all_genus), " 行数 =", sum(len(v) for v in vmr_by_genus.values()))

# ---- Plant.tsv 建索引 ----
plant_genus, plant_species_blob = {}, []
plant_hdr = None
if os.path.isfile(PLANT):
    with open(PLANT, newline="", encoding="utf-8", errors="replace") as f:
        rd = csv.reader(f, delimiter="\t")
        plant_hdr = [c.strip().strip('"') for c in next(rd)]
        gi = None
        for cand in ("Genus", "genus", "Virus Genus", "VMR_Genus"):
            if cand in plant_hdr:
                gi = plant_hdr.index(cand)
                break
        n = 0
        for row in rd:
            n += 1
            if gi is not None and len(row) > gi:
                plant_genus.setdefault(row[gi].strip().strip('"'), 0)
                plant_genus[row[gi].strip().strip('"')] += 1
            plant_species_blob.append(row)
    print("Plant.tsv 行数 =", n, " 属列 =", plant_hdr[gi] if gi is not None else "未识别")
    print("Plant.tsv 属数 =", len(plant_genus))

for name in TARGETS + EXTRA:
    tag = " [已撤回参照]" if name in EXTRA else ""
    print("\n######## %s%s" % (name, tag))
    # VMR
    rows = vmr_by_genus.get(name, [])
    if rows:
        print("  VMR: %d 行" % len(rows))
        from collections import Counter
        hostsrc = Counter(r[26].strip().strip('"') for r in rows)
        genome = Counter(r[25].strip().strip('"') for r in rows)
        print("    Genome:", dict(genome))
        print("    Host source:", dict(hostsrc.most_common()))
        for r in rows[:6]:
            print("    - species=%s | family=%s | genome=%s | host_source=%s" %
                  (r[17].strip().strip('"'), r[13].strip().strip('"'), r[25].strip().strip('"'), r[26].strip().strip('"')))
    else:
        print("  VMR: 0 命中")
    # Plant library
    pg = plant_genus.get(name, 0)
    hits = [r for r in plant_species_blob if any(name.lower() in (c or "").lower() for c in r[:6])]
    print("  PlantVirusDB classified_clean/Plant.tsv: 属列命中 %d ; 前 6 列含该名字的行 %d" % (pg, len(hits)))
    for r in hits[:3]:
        print("    -", " | ".join((c or "")[:40] for c in r[:6]))
    # ICTV 属表
    for p, label in ((ICTVG, "ICTV_MSL41_Genus.tsv"), (REFLIB, "final.cluster.ref_info.tsv")):
        if os.path.isfile(p):
            r = subprocess.run(["grep", "-cw", name, p], capture_output=True, text=True)
            print("  %s: grep -cw = %s" % (label, r.stdout.strip()))
    # NCBI names.dmp
    if os.path.isfile(NAMES):
        r = subprocess.run(["grep", "-cw", name, NAMES], capture_output=True, text=True)
        print("  NCBI names.dmp: grep -cw = %s" % r.stdout.strip())
