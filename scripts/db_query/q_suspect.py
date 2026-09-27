import csv
import collections

PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"

# 权威库 Plant.tsv 的合法植物病毒属
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

print(f"权威库合法植物属: {len(auth_gen)}")
print(f"权威库合法植物科: {len(auth_fam)}")

# C9 判为 Plant 的属清单 (来自上一步的统计, 这里重新全量扫描精确版)
base_dir = "/home/zhangwenda/data-test/out/06_HostPrediction/C9_ICTV_result"
import glob, os

gen_plant = collections.Counter()
gen_plant_familylevel = collections.Counter()
for tsv in glob.glob(os.path.join(base_dir, "*.classified.tsv")):
    with open(tsv, errors="replace") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if (row.get("Predicted_Host", "") or "").strip() != "Plant":
                continue
            g = (row.get("Genus", "") or "").strip()
            gen_plant[g] += 1
            d = (row.get("Determination_Level", "") or "").split("(")[0].strip()
            if d in ("Family", "Order", "None"):
                gen_plant_familylevel[g] += 1

print(f"\nC9 判 Plant 的不同属总数: {len(gen_plant)}")

# 哪些属不在权威库中 (可疑黑名单)
suspect = []
for g, n in gen_plant.items():
    if g and g not in ("NA", "nan", "") and g not in auth_gen:
        suspect.append((g, n, gen_plant_familylevel.get(g, 0)))

suspect.sort(key=lambda x: -x[1])
print(f"\n=== 【不在权威库 Plant.tsv 中】的属: {len(suspect)} 个 (按条数排序) ===")
print(f"{'属':<28} {'C9总数':<8} {'科级判定数':<10}")
tot = 0
for g, n, fl in suspect:
    tot += n
    print(f"{g:<28} {n:<8} {fl:<10}")
print(f"\n合计: {tot} 条")
