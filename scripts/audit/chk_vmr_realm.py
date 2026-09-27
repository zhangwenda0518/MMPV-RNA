#!/usr/bin/env python3
"""核对 VMR_MSL41 各层级列的实际填充率, 解释 Realm 为何弃票"""
import csv
from collections import Counter
from pathlib import Path

VMR = Path("/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")
RANKS = ["Realm", "Subrealm", "Kingdom", "Subkingdom", "Phylum", "Subphylum",
         "Class", "Subclass", "Order", "Suborder", "Family", "Subfamily",
         "Genus", "Subgenus", "Species", "Genome"]

rows = list(csv.DictReader(open(VMR, newline="", encoding="utf-8", errors="replace"), delimiter="\t"))
print("文件: %s" % VMR)
print("行数: %d   列数: %d" % (len(rows), len(rows[0])))
print("\n%-12s %8s %8s %8s" % ("层级", "非空", "空/NA", "唯一值"))
for r in RANKS:
    if r not in rows[0]:
        print("%-12s  <列不存在>" % r)
        continue
    vals = [(x.get(r) or "").strip() for x in rows]
    nb = sum(1 for v in vals if v and v.upper() != "NA")
    uv = len({v for v in vals if v and v.upper() != "NA"})
    print("%-12s %8d %8d %8d" % (r, nb, len(vals) - nb, uv))

print("\nRealm 列取值分布 top10:")
for k, v in Counter((x.get("Realm") or "").strip() for x in rows).most_common(10):
    print("   %-30s %5d" % (k if k else "(空)", v))
print("\nKingdom 列取值分布 top10:")
for k, v in Counter((x.get("Kingdom") or "").strip() for x in rows).most_common(10):
    print("   %-30s %5d" % (k if k else "(空)", v))
print("\nGenome 列取值分布 top12:")
for k, v in Counter((x.get("Genome") or "").strip() for x in rows).most_common(12):
    print("   %-30s %5d" % (k if k else "(空)", v))
print("\n表头前 24 列: %s" % list(rows[0].keys())[:24])
