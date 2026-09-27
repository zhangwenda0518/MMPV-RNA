import glob
import os
import csv
import collections

base = "/home/zhangwenda/data-test/out/06_HostPrediction/C9_ICTV_result"

lvl = collections.Counter()
lvl_plant = collections.Counter()
genus_agree_1 = 0
total = 0
pass_conf_false = 0

for tsv in glob.glob(os.path.join(base, "*.classified.tsv")):
    with open(tsv) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            total += 1
            d = (row.get("Determination_Level", "") or "").split("(")[0].strip()
            lvl[d] += 1
            if (row.get("Predicted_Host", "") or "").strip() == "Plant":
                lvl_plant[d] += 1
            ga = row.get("Genus_agree", "") or ""
            if ga.startswith("1/"):
                genus_agree_1 += 1
            if (row.get("_pass_conf", "") or "").strip().lower() == "false":
                pass_conf_false += 1

print(f"C9 总记录: {total}")
print(f"\nDetermination_Level 分布:")
for k, v in lvl.most_common():
    print(f"  {k}: {v} ({v/total*100:.1f}%)")
print(f"\n其中 Predicted_Host==Plant 的 Determination_Level:")
tp = sum(lvl_plant.values())
for k, v in lvl_plant.most_common():
    print(f"  {k}: {v} ({v/tp*100:.1f}%)")
print(f"\nGenus_agree 仅 1 个工具支持: {genus_agree_1} ({genus_agree_1/total*100:.1f}%)")
print(f"_pass_conf == false: {pass_conf_false} ({pass_conf_false/total*100:.1f}%)")
