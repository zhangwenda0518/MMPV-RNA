import csv
import collections

PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"

# 统计 Plant.tsv 中所有病毒科 + 属，作为"权威植物病毒集合"
fam_plant = collections.Counter()
genus_plant = collections.Counter()
with open(PLANT) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        lin = row.get("Virus_lineage", "") or ""
        parts = [p for p in lin.split(";") if p]
        # 找 Family / Genus (位置因谱系深度而异, 用后缀匹配)
        for i, p in enumerate(parts):
            if p.endswith("viridae") or p.endswith("virinae"):
                if p.endswith("viridae"):
                    fam_plant[p] += 1
            if p.endswith("virus") and i == len(parts) - 3:
                genus_plant[p] += 1

print(f"Plant.tsv 中出现的病毒科: {len(fam_plant)} 个")
print("top20:")
for k, v in fam_plant.most_common(20):
    print(f"  {k}: {v}")

print(f"\nPlant.tsv 中出现的病毒属: {len(genus_plant)} 个")
print("top20:")
for k, v in genus_plant.most_common(20):
    print(f"  {k}: {v}")

# 反查: 我们管线里判为 Plant 的科属, 有多少在权威库里不存在
print("\n=== 核对: 管线 HQ_plant_viruses 里的科是否在权威库中 ===")
FA = "/home/zhangwenda/data-test/out/09_Virome_Analysis/HQ_plant_viruses.fasta"
ids = set()
with open(FA) as f:
    for line in f:
        if line.startswith(">"):
            ids.add(line[1:].strip())
print(f"HQ_plant_viruses: {len(ids)} 条")
