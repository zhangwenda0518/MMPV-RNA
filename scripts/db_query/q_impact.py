import csv
import collections

# 9a 的植物病毒 fasta 有多少条
FA = "/home/zhangwenda/data-test/out/09_Virome_Analysis/HQ_plant_viruses.fasta"
ids = []
with open(FA) as f:
    for line in f:
        if line.startswith(">"):
            ids.append(line[1:].strip())
print(f"HQ_plant_viruses.fasta: {len(ids)} 条")

id_set = set(ids)

# 从 ensemble_host_summary 取这些 contig 的判定级别
# 注意: 判定级别来自 C9, 这里用 Final_Host + Host_ICTV 反查
F = "/home/zhangwenda/data-test/out/06_HostPrediction/ensemble_host_summary.tsv"
fin = collections.Counter()
plant_by_family = collections.Counter()
n = 0
with open(F) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        cid = row.get("contig_id", "")
        if cid not in id_set:
            continue
        n += 1
        fin[row.get("Final_Host", "")] += 1
        if row.get("Final_Host", "") == "Plant":
            plant_by_family[row.get("Family", "")] += 1

print(f"\nensemble_host_summary 中命中 {n} 条")
print("\nFinal_Host 分布:")
for k, v in fin.most_common(10):
    print(f"  {k!r}: {v}")

print("\n植物病毒里 Family 分布 top20:")
for k, v in plant_by_family.most_common(20):
    print(f"  {k!r}: {v}")
