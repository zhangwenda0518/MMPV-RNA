import glob
import os
import csv
import collections

base = "/home/zhangwenda/data-test/out/06_HostPrediction/C9_ICTV_result"

# 科级判定 + Plant 宿主 → 属于哪些科/属
fam_of_family_level = collections.Counter()
genus_of_family_level = collections.Counter()
# 所有 Plant 判定 (不分级别) 的科
fam_all_plant = collections.Counter()
genus_all_plant = collections.Counter()

for tsv in glob.glob(os.path.join(base, "*.classified.tsv")):
    with open(tsv, errors="replace") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if (row.get("Predicted_Host", "") or "").strip() != "Plant":
                continue
            fam = (row.get("Family", "") or "").strip()
            gen = (row.get("Genus", "") or "").strip()
            fam_all_plant[fam] += 1
            genus_all_plant[gen] += 1
            d = (row.get("Determination_Level", "") or "").split("(")[0].strip()
            if d in ("Family", "Order", "None"):
                fam_of_family_level[fam] += 1
                genus_of_family_level[gen] += 1

tot = sum(fam_all_plant.values())
print(f"C9 判为 Plant 的总数: {tot}")
print(f"\n=== Plant 判定中的科 top25 ===")
for k, v in fam_all_plant.most_common(25):
    print(f"  {k!r}: {v}")

print(f"\n=== 【科级/目级判定】的科 top25 (候选黑名单) ===")
for k, v in fam_of_family_level.most_common(25):
    print(f"  {k!r}: {v}")

print(f"\n=== 【科级/目级判定】的属 top25 (候选黑名单) ===")
for k, v in genus_of_family_level.most_common(25):
    print(f"  {k!r}: {v}")
