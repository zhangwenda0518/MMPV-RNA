#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""12 个 ICTV 宿主含植物科的属级细查 v2 (Plant.tsv 位置口径为准) + 黑名单/白名单交叉校验 + 影响实测。

口径统一: Virus_lineage 按 ';' 分 9 段, index5=Family, index6=Genus (与管线 3_host_probability.py 一致)
证据层级:
  E1 Plant.tsv 有该属记录                      → 植物属 (不得拉黑)
  E2 属级概率表 Plant_Records>0                → 植物属 (不得拉黑)
  E3 属级概率表存在且 Plant_Records=0 且非 Plant → 非植物属 (可拉黑, 需看是否真在结果里)
  E4 两张表都无该属                            → 证据不足 (不动)
  E5 Plant.tsv 无该属但科 ICTV 宿主含植物      → DB 缺口 (标注, 不拉黑)
"""
import importlib.util
import os
import sys
from collections import Counter, defaultdict

import pandas as pd

PLANT = os.path.expanduser("~/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv")
HP = "/home/zhangwenda/MMPV-RNA/database/cross_analysis"
TAX = {
    "goji": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/"
            "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv",
    "onekp": "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus/"
             "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv",
}
TARGET = ["Artoviridae", "Chrysoviridae", "Genomoviridae", "Kanorauviridae", "Mitoviridae",
          "Ourmiaviridae", "Pestiviridae", "Pseudoviridae", "Spiciviridae", "Tomosaviridae",
          "Virgaviridae", "Solemoviridae"]
PREV = {"Artoviridae", "Chrysoviridae", "Genomoviridae", "Kanorauviridae", "Mitoviridae",
        "Ourmiaviridae", "Pestiviridae", "Pseudoviridae", "Spiciviridae"}
ICTV_HOST = {"Artoviridae": "plants, invertebrates", "Chrysoviridae": "fungi, plants, invertebrates",
             "Genomoviridae": "fungi, plants, invertebrates, vertebrates",
             "Kanorauviridae": "plants, invertebrates, vertebrates", "Mitoviridae": "fungi, plants",
             "Ourmiaviridae": "invertebrates, plants",
             "Pestiviridae": "vertebrates, invertebrates, plants",
             "Pseudoviridae": "protists, fungi, plants, invertebrates",
             "Spiciviridae": "plants, invertebrates", "Tomosaviridae": "plants",
             "Virgaviridae": "plants", "Solemoviridae": "plants"}

# ---------- Plant.tsv (位置口径) ----------
pl_gen, pl_fam, fam2gen = Counter(), Counter(), defaultdict(Counter)
with open(PLANT, encoding="utf-8", errors="replace") as fh:
    fh.readline()
    for line in fh:
        fs = line.rstrip("\n").split("\t")
        if len(fs) < 8:
            continue
        p = fs[3].split(";")
        f = p[5].strip() if len(p) > 5 else ""
        g = p[6].strip() if len(p) > 6 else ""
        if f:
            pl_fam[f] += 1
        if g:
            pl_gen[g] += 1
            fam2gen[f][g] += 1
print("Plant.tsv: 科 %d / 属 %d" % (len(pl_fam), len(pl_gen)))

genhp = pd.read_csv(os.path.join(HP, "genus_host_probability.tsv"), sep="\t", low_memory=False).set_index("Genus")
famhp = pd.read_csv(os.path.join(HP, "family_host_probability.tsv"), sep="\t", low_memory=False).set_index("Family")
print("属级概率表 %d 属 ; 科级概率表 %d 科" % (len(genhp), len(famhp)))

P = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
FAMS = set(m.NON_PLANT_FAMILIES_FALLBACK)
GENS = set(m.NON_PLANT_GENERA)
print("现装黑名单: 科 %d / 属 %d" % (len(FAMS), len(GENS)))

# ---------- 交叉校验: 黑名单里的植物属/科 ----------
bad_g = sorted(g for g in GENS if pl_gen.get(g, 0) > 0 or
               (g in genhp.index and genhp.loc[g, "Plant_Records"] > 0))
bad_f = sorted(f for f in FAMS if pl_fam.get(f, 0) > 0 or
               (f in famhp.index and famhp.loc[f, "Plant_Records"] > 0))
print("\n[校验1] 黑名单中带植物证据的属: %s" % (bad_g or "无"))
print("[校验2] 黑名单中带植物证据的科: %s" % (bad_f or "无"))

# ---------- 各树属级情况 ----------
tax = {}
for tag, path in TAX.items():
    tax[tag] = pd.read_csv(path, sep="\t", usecols=["Family", "Genus"], low_memory=False) if os.path.isfile(path) else None
    print("05 表 %s: %s" % (tag, "缺失" if tax[tag] is None else "%d 行" % len(tax[tag])))

rows, cand = [], []
for fam in TARGET:
    print("\n" + "=" * 92)
    print("【%s】%s | ICTV Host=%s | Plant.tsv 行=%d 属=%d | 科表=%s"
          % (fam, ("旧名单" if fam in PREV else "新出现"), ICTV_HOST[fam], pl_fam.get(fam, 0),
             len(fam2gen.get(fam, {})),
             (famhp.loc[fam, "Predicted_Host"] + "/Plant=" + str(famhp.loc[fam, "Plant_Records"]))
             if fam in famhp.index else "无条目"))
    seen = {}
    for tag in tax:
        if tax[tag] is None:
            continue
        sub = tax[tag][tax[tag]["Family"] == fam]
        for g, n in Counter(sub["Genus"].fillna("NA").astype(str)).items():
            if g in ("NA", "nan", ""):
                continue
            d = seen.setdefault(g, {"rows": {}, "plant": pl_gen.get(g, 0),
                                    "fam_plant_genus": g in fam2gen.get(fam, {})})
            d["rows"][tag] = n
    for g, d in sorted(seen.items(), key=lambda kv: -sum(kv[1]["rows"].values())):
        gh = genhp.loc[g] if g in genhp.index else None
        if d["plant"] > 0 or (gh is not None and gh["Plant_Records"] > 0):
            v = "E1/E2 植物属-白名单"
        elif gh is not None and gh["Predicted_Host"] != "Plant":
            v = "E3 非植物属-候选黑名单"
        elif gh is None:
            v = "E4 证据不足-不动"
        else:
            v = "E5 DB缺口-不拉黑"
        print("   %-22s 05行=%s 植物库记录=%-5d 属表=%-12s 总=%-6s %s | %s"
              % (g, d["rows"], d["plant"],
                 (gh["Predicted_Host"] if gh is not None else "无条目"),
                 (int(gh["Total_Records"]) if gh is not None else "-"), v,
                 ("已在黑名单" if g in GENS else "")))
        if v == "E3 非植物属-候选黑名单":
            cand.append({"Family": fam, "Genus": g, "Rows_goji": d["rows"].get("goji", 0),
                         "Rows_onekp": d["rows"].get("onekp", 0),
                         "GHP_Predicted_Host": gh["Predicted_Host"],
                         "GHP_Total_Records": int(gh["Total_Records"]),
                         "GHP_Confidence": gh["Confidence_Level"],
                         "Verdict": v})
        rows.append({"Family": fam, "Genus": g, "Rows_goji": d["rows"].get("goji", 0),
                     "Rows_onekp": d["rows"].get("onekp", 0), "Plant_tsv_records": d["plant"],
                     "Genus_HP_Predicted_Host": (gh["Predicted_Host"] if gh is not None else "无条目"),
                     "Genus_HP_Plant_Records": (int(gh["Plant_Records"]) if gh is not None else -1),
                     "Genus_HP_Total_Records": (int(gh["Total_Records"]) if gh is not None else -1),
                     "Genus_HP_Confidence": (gh["Confidence_Level"] if gh is not None else "NA"),
                     "Current_Blacklist": g in GENS, "Verdict": v})

out = pd.DataFrame(rows)
out.to_csv("/tmp/PLANT_FAMILY_GENUS_TRIAGE.tsv", sep="\t", index=False)
print("\n写出 /tmp/PLANT_FAMILY_GENUS_TRIAGE.tsv (%d 行)" % len(out))
print(out.groupby("Verdict")["Genus"].count().to_string())
print("\n候选黑名单 (出现在 12 科结果里的非植物属):")
cd = pd.DataFrame(cand)
if len(cd):
    print(cd.to_string(index=False))
