import csv
import collections

F = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/host_extract/Final_Virus_Host_Lineage.tsv"

# 先看表头
with open(F, errors="replace") as f:
    head = f.readline().rstrip("\n").split("\t")
print("列:", head)

fams = ["Mimiviridae", "Metaviridae", "Biavirus", "Partitiviridae"]
stats = {k: collections.Counter() for k in fams}
cnt = {k: 0 for k in fams}

with open(F, errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        vals = list(row.values())
        blob = "\t".join(str(v) for v in vals)
        for k in fams:
            if k in blob:
                cnt[k] += 1
                # 找 category 类列
                for col in row:
                    if "ategor" in col:
                        stats[k][row[col]] += 1
                break

for k in fams:
    print(f"\n=== {k}: 命中 {cnt[k]} 行 ===")
    for kk, vv in stats[k].most_common(10):
        print(f"    {kk!r}: {vv}")
