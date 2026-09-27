import csv
import collections
import glob
import os

base_dir = "/home/zhangwenda/data-test/out/06_HostPrediction/C9_ICTV_result"
PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"

# 权威库所有植物病毒名 (小写, 去空格)
auth_names = set()
auth_species = set()
with open(PLANT, errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        nm = (row.get("Virus_Name", "") or "").strip().lower()
        if nm:
            auth_names.add(nm)
        for p in (row.get("Virus_lineage", "") or "").split(";"):
            p = p.strip()
            if p:
                auth_species.add(p.lower())

print(f"权威库植物病毒名: {len(auth_names)}, 谱系 token: {len(auth_species)}")

# Mimiviridae 判 Plant 的 244 条: 它们的 Species 是否在权威库中
rows = []
for tsv in glob.glob(os.path.join(base_dir, "*.classified.tsv")):
    with open(tsv, errors="replace") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if (row.get("Predicted_Host", "") or "").strip() != "Plant":
                continue
            if (row.get("Family", "") or "").strip() != "Mimiviridae":
                continue
            rows.append(row)

in_auth = 0
gen_in_auth = 0
for row in rows:
    sp = (row.get("Species", "") or "").strip().lower()
    ge = (row.get("Genus", "") or "").strip()
    if sp and (sp in auth_names or sp in auth_species):
        in_auth += 1
    if ge and ge.lower() in auth_species:
        gen_in_auth += 1

print(f"\nMimiviridae 判 Plant 共 {len(rows)} 条:")
print(f"  Species 在权威库中: {in_auth} ({in_auth/len(rows)*100:.1f}%)")
print(f"  Genus 在权威库中:   {gen_in_auth} ({gen_in_auth/len(rows)*100:.1f}%)")

print("\n=== 结论验证: 这些 Species 的真科 ===")
# 用权威库反查这些 species 的真实分类
sp_to_fam = {}
with open(PLANT, errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        for p in (row.get("Virus_lineage", "") or "").split(";"):
            p = p.strip()
            if p:
                sp_to_fam[p.lower()] = row.get("Virus_lineage", "")

shown = 0
for row in rows:
    sp = (row.get("Species", "") or "").strip()
    if not sp:
        continue
    key = sp.lower()
    if key in sp_to_fam and shown < 5:
        lin = sp_to_fam[key]
        real_fam = [p for p in lin.split(";") if p.endswith("viridae")]
        print(f"  {sp:36} → 权威库真实科: {real_fam}")
        shown += 1
