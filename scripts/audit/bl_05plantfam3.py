"""精确 token 判定这几个科的宿主归属 + 量化被置空/被漂移的属中有多少是植物病毒属.

子串陷阱: botourmiaviridae 含 ourmiaviridae; 必须按分号切分后精确比对.
"""
import csv
import os
import re
from collections import Counter, defaultdict

MSL = "/home/zhangwenda/database/virus-db/acvirus_db/taxa.txt"
REF = os.path.expanduser("~/database/taxonomy/genus_family_ref.tsv")
PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"
BASE = ("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
        "/05_Taxonomy/Votus.integrated")
CAL = f"{BASE}/calibration_20260914/final_integrated_classification.calibrated.tsv"

FAMS7 = ["Hepaciviridae", "Pestiviridae", "Zimmerviridae", "Ambiguiviridae",
         "Autographiviridae", "Ourmiaviridae", "Parahypoviridae"]
ALSO = ["Botourmiaviridae", "Partitiviridae", "Mimiviridae", "Phycodnaviridae", "Marseilleviridae",
        "Pithoviridae", "Iridoviridae", "Schizomimiviridae", "Retroviridae", "Astroviridae",
        "Orpheoviridae", "Epsomviridae", "Mycoalphaviridae", "Hydriviridae", "Flaviviridae"]
ALL = FAMS7 + ALSO

# ── 已知属名（用两套参照建全集） ─────────────────────────────────────
known_genus = set()
with open(MSL, errors="replace") as f:
    for r in csv.DictReader(f):
        g = (r.get("Genus") or "").strip().lower()
        if g:
            known_genus.add(g)
with open(REF, errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        g = (r.get("Genus") or "").strip().lower()
        if g:
            known_genus.add(g)
print(f"已知属名全集: {len(known_genus)}")

# ── Plant.tsv: 精确 token 抽科名/属名, 记宿主 ────────────────────────
fam_cat = defaultdict(Counter)
fam_virus = defaultdict(set)
plant_genus = set()
n = 0
with open(PLANT, errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        n += 1
        lin = (r.get("Virus_lineage") or "")
        toks = [t.strip().lower() for t in re.split(r"[;,]", lin) if t.strip()]
        cat = (r.get("Host_Category") or "").strip()
        for t in toks:
            if t in known_genus:
                plant_genus.add(t)
            for k in ALL:
                if t == k.lower():
                    fam_cat[k][cat] += 1
                    fam_virus[k].add(r.get("Virus_Name") or r.get("Accession"))
print(f"Plant.tsv {n} 行; 其中作为'属'出现的属名 {len(plant_genus)} 个")

print("\n=== 精确 token: 这些科在植物病毒库谱系里是否出现 ===")
for k in ALL:
    v, c = len(fam_virus.get(k, ())), dict(fam_cat.get(k, {}))
    mark = "植物" if c.get("Plant") else ("非植物" if c else "库中无")
    print(f"  {k:<22} 病毒 {v:>4} 宿主 {str(c):<22} -> {mark}")


def val(x):
    x = (x or "").strip().strip('"')
    return None if x in ("", "NA", "N/A", "-") else x


rows = list(csv.DictReader(open(CAL, errors="replace"), delimiter="\t"))
print(f"\n枸杞校准表 {len(rows)} 行")

# ── 被置空的属: 有多少是植物病毒属 ───────────────────────────────────
blank_plant = Counter()
blank_fam = Counter()
for r in rows:
    if val(r.get("calib_action")) != "blank":
        continue
    g = (val(r.get("calib_prev_genus")) or "").lower()
    if g in plant_genus:
        blank_plant["prev_genus 是植物病毒属"] += 1
        blank_fam[val(r.get("Family")) or "(空)"] += 1
    else:
        blank_plant["prev_genus 非植物病毒属/未知"] += 1
print("\n=== 1246 blank 行: 被置空的属是不是植物病毒属 ===")
for k, v in blank_plant.items():
    print(f"  {k:<30} {v} 行")
print("  这些植物属对应的行内科名 (前 10):")
for k, v in blank_fam.most_common(10):
    print(f"    {k:<26} {v} 行")

# ── 现存活行里: 科非植物科 但 属是植物病毒属 ────────────────────────
mismatch = Counter()
for r in rows:
    fam = val(r.get("Family"))
    g = (val(r.get("Genus")) or "").lower()
    if not fam or not g:
        continue
    if fam not in fam_cat or not fam_cat[fam].get("Plant"):
        if g in plant_genus:
            mismatch[(fam, g)] += 1
print(f"\n=== 科非植物科 但属是植物病毒属 的行: {sum(mismatch.values())} 行 ===")
for (a, b), v in mismatch.most_common(12):
    print(f"  科={a:<22} 属={b:<24} {v} 行")
