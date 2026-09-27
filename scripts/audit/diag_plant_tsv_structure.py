#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plant.tsv 结构诊断: 位置口径(index5/6) 与 科锚定口径 的差异, 关键科在两套口径/概率表/黑名单中的状态。"""
import importlib.util
import os
import re
from collections import Counter, defaultdict

PLANT = os.path.expanduser("~/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv")
print("Plant.tsv 存在: %s" % os.path.isfile(PLANT))

pos_fam, anc_fam = Counter(), Counter()
pos_gen, anc_gen = Counter(), Counter()
len_dist = Counter()
cat = Counter()
key_pos, key_anc = defaultdict(list), defaultdict(list)
KEY = ["Mitoviridae", "Mitovirus", "Ourmiaviridae", "Botourmiaviridae", "Ourmiavirus",
       "Tomosaviridae", "Virtovirus", "Tombusviridae", "Partitiviridae", "Endornaviridae",
       "Aspiviridae", "Kanorauviridae", "Shabtayvirus", "Pseudoviridae", "Pseudovirus",
       "Genomoviridae", "Gemycircularvirus", "Spiciviridae", "Spicivirus",
       "Chrysoviridae", "Alphachrysovirus", "Artoviridae", "Pestiviridae", "Orthopestivirus"]

n = 0
with open(PLANT, encoding="utf-8", errors="replace") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    icat = hdr.index("Host_Category")
    for line in fh:
        fs = line.rstrip("\n").split("\t")
        if len(fs) < 4:
            continue
        n += 1
        lin = fs[3]
        parts = lin.split(";")
        len_dist[len(parts)] += 1
        cat[fs[icat] if len(fs) > icat else "?"] += 1
        pf = parts[5].strip() if len(parts) > 5 else ""
        pg = parts[6].strip() if len(parts) > 6 else ""
        af = ag = ""
        for i, t in enumerate(parts):
            t = t.strip()
            if t.endswith("viridae"):
                af = t
                ag = parts[i + 1].strip() if len(parts) > i + 1 else ""
                break
        if pf:
            pos_fam[pf] += 1
        if pg:
            pos_gen[pg] += 1
        if af:
            anc_fam[af] += 1
        if ag:
            anc_gen[ag] += 1
        for k in KEY:
            if k in lin:
                if pf or pg:
                    key_pos[k].append((pf, pg, lin))
                if af:
                    key_anc[k].append((af, ag))

print("行数 %d ; 字段数分布 %s" % (n, len_dist.most_common()))
print("Host_Category %s" % cat.most_common())
print("\n位置口径: 科 %d / 属 %d" % (len(pos_fam), len(pos_gen)))
print("锚定口径: 科 %d / 属 %d" % (len(anc_fam), len(anc_gen)))
print("\n差异科 (锚定有而位置无, 前 25): %s" % sorted(set(anc_fam) - set(pos_fam))[:25])
print("差异科 (位置有而锚定无, 前 25): %s" % sorted(set(pos_fam) - set(anc_fam))[:25])

print("\n关键名检查:")
for k in KEY:
    p = len(key_pos.get(k, []))
    a = len(key_anc.get(k, []))
    ex = (key_pos.get(k) or key_anc.get(k) or [("", "", "")])[0]
    print("  %-20s 位置口径行=%-6d 锚定口径行=%-6d 例: 位置(%s|%s) 锚定(%s|%s)"
          % (k, p, a, ex[0] or "-", ex[1] or "-",
             (key_anc.get(k) or [("", "")])[0][0] if key_anc.get(k) else "-",
             (key_anc.get(k) or [("", "")])[0][1] if key_anc.get(k) else "-"))

# 概率表 + 黑名单
HP = "/home/zhangwenda/MMPV-RNA/database/cross_analysis"
if not os.path.isdir(HP):
    HP = os.path.expanduser("~/plant_virus_db/1.virus-host_db/C-host_classify/cross_analysis")
print("\n概率表目录: %s" % HP)
import pandas as pd
famhp = pd.read_csv(os.path.join(HP, "family_host_probability.tsv"), sep="\t", low_memory=False)
genhp = pd.read_csv(os.path.join(HP, "genus_host_probability.tsv"), sep="\t", low_memory=False)
print("科级表 %d 行 ; 属级表 %d 行" % (len(famhp), len(genhp)))
fi, gi = famhp.set_index("Family"), genhp.set_index("Genus")
for k in ["Mitoviridae", "Botourmiaviridae", "Ourmiaviridae", "Tomosaviridae", "Tombusviridae",
          "Kanorauviridae", "Pseudoviridae", "Genomoviridae", "Spiciviridae", "Artoviridae",
          "Pestiviridae", "Chrysoviridae", "Virgaviridae", "Solemoviridae", "Partitiviridae"]:
    if k in fi.index:
        r = fi.loc[k]
        print("  [科表] %-18s %-10s Plant=%-6s 总=%-7s %s"
              % (k, r["Predicted_Host"], r["Plant_Records"], r["Total_Records"], r["Confidence_Level"]))
    else:
        print("  [科表] %-18s 无条目" % k)
for k in ["Mitovirus", "Duamitovirus", "Dagavirus", "Ourmiavirus", "Virtovirus", "Shabtayvirus",
          "Pseudovirus", "Sirevirus", "Hemivirus", "Gemycircularvirus", "Spicivirus",
          "Orthopestivirus", "Alphachrysovirus", "Tobamovirus", "Polerovirus", "Sobemovirus",
          "Enamovirus", "Polemovirus", "Tombusvirus", "Alphanecrovirus", "Betanecrovirus",
          "Machlomovirus", "Panicovirus", "Carmovirus", "Alphacarmovirus", "Betacarmovirus",
          "Gammacarmovirus", "Dianthovirus", "Aureusvirus", "Pelarspovirus", "Alohovirus",
          "Albetovirus", "Avenavirus", "Zeavirus", "Rosadnavirus", "Gallantivirus", "Umbravirus",
          "Polemovirus", "Goravirus", "Furovirus", "Hordeivirus", "Pecluvirus", "Pomovirus",
          "Tobravirus", "Biavirus", "Miraophiovirus", "Efemunavirus", "Tungrovirus", "Caulimovirus"]:
    if k in gi.index:
        r = gi.loc[k]
        print("  [属表] %-18s %-10s Plant=%-6s 总=%-7s %s"
              % (k, r["Predicted_Host"], r["Plant_Records"], r["Total_Records"], r["Confidence_Level"]))
    else:
        print("  [属表] %-18s 无条目" % k)

P = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
FAMS = set(m.NON_PLANT_FAMILIES_FALLBACK)
GENS = set(m.NON_PLANT_GENERA)
print("\n现装黑名单: 科 %d / 属 %d ; FAMILY_FIRST_VETO=%s" % (len(FAMS), len(GENS), m.FAMILY_FIRST_VETO))
for k in ["Mitoviridae", "Botourmiaviridae", "Ourmiaviridae", "Tomosaviridae", "Tombusviridae",
          "Kanorauviridae", "Virgaviridae", "Solemoviridae", "Partitiviridae", "Aspiviridae"]:
    print("  科在名单: %-18s %s" % (k, k in FAMS))
for k in ["Mitovirus", "Duamitovirus", "Dagavirus", "Ourmiavirus", "Virtovirus", "Shabtayvirus",
          "Pseudovirus", "Sirevirus", "Hemivirus", "Tobamovirus", "Polerovirus", "Sobemovirus",
          "Enamovirus", "Biavirus", "Efemunavirus", "Furovirus", "Pomovirus", "Tobravirus"]:
    print("  属在名单: %-18s %s" % (k, k in GENS))
