import csv
import collections

# 本项目 PlantVirusDB 参考库的真实宿主注释
PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"

fams = ["Mimiviridae", "Metaviridae", "Partitiviridae"]

# 每个科: accession -> virus_name / host / lineage
for fam in fams:
    n = 0
    hosts = collections.Counter()
    gens = collections.Counter()
    examples = []
    with open(PLANT) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            lin = row.get("Virus_lineage", "") or ""
            if fam not in lin:
                continue
            n += 1
            hosts[row.get("Host_Name", "")] += 1
            # 提取属
            parts = lin.split(";")
            if len(parts) >= 7:
                gens[parts[6]] += 1
            if len(examples) < 3:
                examples.append((row.get("Accession", ""), row.get("Virus_Name", ""),
                                 row.get("Host_Name", ""), row.get("Host_Category", "")))
    print(f"\n=== {fam} 在 Plant.tsv 中: {n} 条 ===")
    if n:
        print("  病毒属分布:")
        for k, v in gens.most_common(8):
            print(f"    {k}: {v}")
        print("  宿主 top8:")
        for k, v in hosts.most_common(8):
            print(f"    {k!r}: {v}")
        print("  示例:")
        for a, nm, h, hc in examples:
            print(f"    {a} | {nm[:40]} | host={h} | cat={hc}")
