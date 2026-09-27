import csv
import collections
import glob
import os

base_dir = "/home/zhangwenda/data-test/out/06_HostPrediction/C9_ICTV_result"

# Mimiviridae 判 Plant 的明细
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

print(f"Mimiviridae 判 Plant: {len(rows)} 条")
print("\n=== 前 10 条 ===")
for row in rows[:10]:
    print(f"  Genus={row.get('Genus','')!r:22} Species={row.get('Species','')!r:34} "
          f"Level={row.get('Determination_Level','')!r:22} conf={row.get('Integrated_Confidence','')}")

# Level 分布
lvl = collections.Counter(r.get("Determination_Level", "") for r in rows)
print("\nLevel 分布:", dict(lvl))
gen = collections.Counter(r.get("Genus", "") for r in rows)
print("\nGenus 分布 top10:", dict(gen.most_common(10)))

# 权威库核对: Mimiviridae 的宿主
print("\n=== 对照: 权威库 Mimiviridae 宿主 ===")
FIN = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/host_extract/Final_Virus_Host_Lineage.tsv"
hosts = collections.Counter()
with open(FIN, errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        if "Mimiviridae" in (row.get("Virus_lineage", "") or ""):
            hosts[row.get("Host_Name", "")] += 1
            if sum(hosts.values()) > 5000:
                break
for k, v in hosts.most_common(6):
    print(f"  {k}: {v}")
