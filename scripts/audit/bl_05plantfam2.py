"""这几个科是不是植物科: 用 Plant.tsv 的 Virus_lineage(唯一能查科的字段) + 下游植物表实际落行."""
import csv
import os
from collections import Counter, defaultdict

FAMS7 = ["Hepaciviridae", "Pestiviridae", "Zimmerviridae", "Ambiguiviridae",
         "Autographiviridae", "Ourmiaviridae", "Parahypoviridae"]
ALSO = ["Botourmiaviridae", "Partitiviridae", "Mimiviridae", "Phycodnaviridae", "Marseilleviridae",
        "Pithoviridae", "Iridoviridae", "Schizomimiviridae", "Retroviridae", "Astroviridae",
        "Orpheoviridae", "Epsomviridae", "Mycoalphaviridae", "Hydriviridae", "Flaviviridae"]
ALL = FAMS7 + ALSO
PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"

# ── 1. Plant.tsv 的 Virus_lineage 里到底有没有这些科 ──────────────────
hit_virus = defaultdict(set)          # 科 -> 病毒名集合
hit_cat = defaultdict(Counter)        # 科 -> 宿主类别计数
hit_example = defaultdict(list)
n = 0
with open(PLANT, errors="replace") as f:
    rd = csv.DictReader(f, delimiter="\t")
    for r in rd:
        n += 1
        lin = (r.get("Virus_lineage") or "").lower()
        for k in ALL:
            if k.lower() in lin:
                hit_virus[k].add(r.get("Virus_Name") or r.get("Accession"))
                hit_cat[k][(r.get("Host_Category") or "").strip()] += 1
                if len(hit_example[k]) < 4:
                    hit_example[k].append(f"{r.get('Virus_Name')} <- {r.get('Host_Name')}")
print(f"Plant.tsv 行数 {n}")
print("\n=== 这些科在植物病毒库 Virus_lineage 里的出现 ===")
print(f"{'科名':<24}{'病毒数':>6}  宿主类别")
for k in ALL:
    vs = len(hit_virus.get(k, ()))
    cats = dict(hit_cat.get(k, {}))
    print(f"  {k:<22}{vs:>6}  {cats}")
for k in ALL:
    if hit_example.get(k):
        print(f"\n  [{k}] 例: " + " | ".join(hit_example[k]))

# ── 2. 下游植物表里落进这些科的实际行 ────────────────────────────────
DOWN = [
    ("枸杞 Lycium barbarum",
     "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/10_Reports/All_plant.viruses_info.tsv"),
    ("OneKP",
     "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv"),
]
for label, p in DOWN:
    print(f"\n\n########## {label} 下游植物表 ##########")
    with open(p, errors="replace") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    print(f"  总行数 {len(rows)}")
    for k in ALL:
        sub = [r for r in rows if (r.get("Family") or "").strip() == k]
        if not sub:
            continue
        print(f"\n  --- Family={k}: {len(sub)} 行 ---")
        for r in sub[:8]:
            print(f"    {r.get('contig_id')} | G={r.get('Genus')} | S={r.get('Species')} "
                  f"| conf={r.get('confidence')} | tool={r.get('primary_tool')} | src={r.get('source')}")
        if len(sub) > 8:
            print(f"    ... 余 {len(sub) - 8} 行")
