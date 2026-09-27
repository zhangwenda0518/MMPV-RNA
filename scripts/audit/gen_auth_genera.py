#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从权威库 Plant.tsv 提取合法植物病毒属清单 (每行一个属)。"""
import csv

PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"
OUT = "/tmp/auth_plant_genera.txt"

auth_gen = set()
auth_fam = set()
with open(PLANT, errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        for p in (row.get("Virus_lineage", "") or "").split(";"):
            p = p.strip()
            if not p:
                continue
            if p.endswith("viridae"):
                auth_fam.add(p)
            elif p.endswith("virus"):
                auth_gen.add(p)

with open(OUT, "w") as f:
    for g in sorted(auth_gen):
        f.write(g + "\n")

print(f"权威库合法植物属: {len(auth_gen)} 个 → {OUT}")
print(f"权威库合法植物科: {len(auth_fam)} 个")
print("样例:", sorted(auth_gen)[:10])
# 验证关键属是否在库
for g in ["Deltapartitivirus", "Fabavirus", "Nepovirus", "Ipomovirus",
          "Biavirus", "Sylvanvirus", "Potyvirus"]:
    print(f"  {g}: {'✓ 在库' if g in auth_gen else '✗ 不在库'}")
