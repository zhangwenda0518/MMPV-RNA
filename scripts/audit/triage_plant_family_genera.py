#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""12 个 ICTV 官方"宿主含植物"科的属级细查 (结合 Plant.tsv + 管线自带的属级宿主概率表)。

关键口径 (与管线 3_host_probability.py 一致, 避免自造解析):
  Virus_lineage 补齐为 9 段后按位置取:  index5=Family, index6=Genus, index7=Species
  即 parts = (lineage + ';;;;;;;;').split(';')

判据:
  A. 属在 Plant.tsv 有记录 或 genus_host_probability.tsv 的 Plant_Records>0 → 植物属, 不得进黑名单
  B. 属无植物证据且属级宿主预测为非植物 → 可进非植物属黑名单
  C. 属在任何表中都查不到 → 证据不足, 不动 (避免误否决)
"""
import os
from collections import Counter, defaultdict

import pandas as pd

PLANT = os.path.expanduser("~/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv")
HP = os.path.expanduser("~/plant_virus_db/1.virus-host_db/C-host_classify/cross_analysis") \
    if os.path.isdir(os.path.expanduser("~/plant_virus_db/1.virus-host_db/C-host_classify/cross_analysis")) \
    else "/home/zhangwenda/MMPV-RNA/database/cross_analysis"
GENHP = os.path.join(HP, "genus_host_probability.tsv")
FAMHP = os.path.join(HP, "family_host_probability.tsv")
TAXES = {
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

print("概率表目录: %s" % HP)


def split3(lin):
    parts = (str(lin) + ";;;;;;;;").split(";")
    return parts[5].strip(), parts[6].strip(), parts[7].strip()


pl = pd.read_csv(PLANT, sep="\t", low_memory=False, usecols=["Virus_lineage", "Host_Category"])
fam_r, gen_r, fam2gen = Counter(), Counter(), defaultdict(Counter)
for lin in pl["Virus_lineage"]:
    f, g, _ = split3(lin)
    if f:
        fam_r[f] += 1
    if g:
        gen_r[g] += 1
        fam2gen[f][g] += 1
print("Plant.tsv %d 行 → 科 %d / 属 %d (位置口径)" % (len(pl), len(fam_r), len(gen_r)))

genhp = pd.read_csv(GENHP, sep="\t", low_memory=False)
famhp = pd.read_csv(FAMHP, sep="\t", low_memory=False)
genhp_idx = genhp.set_index("Genus")
famhp_idx = famhp.set_index("Family")
print("属级概率表 %d 行 / 科级概率表 %d 行" % (len(genhp), len(famhp)))
print("属级表中 Plant_Records>0 的属数: %d"
      % int((genhp["Plant_Records"] > 0).sum()))
inter = sorted(set(genhp_idx.index) & set(gen_r))
print("Plant.tsv 属 ∩ 属级概率表属 = %d ; Plant.tsv 属中不在概率表的: %s"
      % (len(inter), sorted(set(gen_r) - set(genhp_idx.index))[:12]))

tax_gens = {}
for tag, path in TAXES.items():
    if not os.path.isfile(path):
        print("  [跳过] %s 05 表不存在" % tag)
        continue
    d = pd.read_csv(path, sep="\t", usecols=["Family", "Genus"], low_memory=False)
    tax_gens[tag] = d

rows = []
for fam in TARGET:
    tag_prev = "旧" if fam in PREV else "新"
    fr = fam_r.get(fam, 0)
    fhp = famhp_idx.loc[fam] if fam in famhp_idx.index else None
    print("\n" + "=" * 96)
    print("【%s】%s   Plant.tsv 行=%d ; 植物库属 %d 个"
          % (fam, tag_prev, fr, len(fam2gen.get(fam, {}))))
    if fhp is not None:
        print("  科级概率表: Predicted_Host=%s Plant_Records=%s 置信=%s"
              % (fhp["Predicted_Host"], fhp["Plant_Records"], fhp["Confidence_Level"]))
    else:
        print("  科级概率表: 无该科条目")
    top = fam2gen.get(fam, Counter()).most_common(8)
    print("  植物库属 Top: %s" % (top if top else "无"))
    for tag, d in tax_gens.items():
        sub = d[d["Family"] == fam]
        gc = Counter(sub["Genus"].fillna("NA").astype(str))
        print("  %s 05 表: %d 行 ; 属分布 %s" % (tag, len(sub), gc.most_common(10)))
        for g, n in gc.items():
            if g in ("NA", "nan", ""):
                continue
            ghp = genhp_idx.loc[g] if g in genhp_idx.index else None
            rows.append({
                "Family": fam, "Family_Is_New": tag_prev, "Source": tag, "Genus": g,
                "Rows_in_05": n,
                "Plant_tsv_records": gen_r.get(g, 0),
                "In_Plant_tsv_family": fam in fam2gen and g in fam2gen[fam],
                "GHP_Predicted_Host": (ghp["Predicted_Host"] if ghp is not None else "无记录"),
                "GHP_Plant_Records": (int(ghp["Plant_Records"]) if ghp is not None else -1),
                "GHP_Total_Records": (int(ghp["Total_Records"]) if ghp is not None else -1),
                "GHP_Confidence": (ghp["Confidence_Level"] if ghp is not None else "NA"),
                "Verdict": (
                    "植物属-不得拉黑" if (gen_r.get(g, 0) > 0 or
                                     (ghp is not None and ghp["Plant_Records"] > 0))
                    else ("非植物属-可拉黑" if ghp is not None else "证据不足-不动")),
            })

out = pd.DataFrame(rows).drop_duplicates(subset=["Family", "Genus"])
outp = "/home/zhangwenda/MMPV-RNA/scripts/audit/PLANT_FAMILY_GENUS_TRIAGE.tsv"
out.to_csv(outp, sep="\t", index=False)
print("\n写出: %s (%d 行)" % (outp, len(out)))
print("\n判定汇总:")
print(out.groupby(["Verdict"])["Genus"].count().to_string())
print("\n【候选黑名单】(非植物属, 出现在这些科的 05 结果里):")
cand = out[(out["Verdict"] == "非植物属-可拉黑") & (out["Rows_in_05"] > 0)]
for r in cand.itertuples(index=False):
    print("   %-18s %-20s 05行=%-5d 属级预测=%-10s 植物records=%-4d 置信=%s"
          % (r.Family, r.Genus, r.Rows_in_05, r.GHP_Predicted_Host, r.GHP_Plant_Records,
             r.GHP_Confidence))
print("\n【证据不足】(不动): %s" % sorted(set(out[out["Verdict"] == "证据不足-不动"]["Genus"])))
