import csv
import collections

F = "/home/zhangwenda/data-test/out/06_HostPrediction/ensemble_host_summary.tsv"
rows_bia_plant = []
method_cnt = collections.Counter()
ictv_cnt = collections.Counter()

with open(F) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        if (row.get("Genus", "") or "").strip() != "Biavirus":
            continue
        if row.get("Final_Host", "") != "Plant":
            continue
        rows_bia_plant.append(row)
        method_cnt[row.get("Decision_Method", "")] += 1
        ictv_cnt[row.get("Host_ICTV", "")] += 1

print(f"Biavirus + Final_Host==Plant 总数: {len(rows_bia_plant)}")
print("\nDecision_Method 分布:")
for k, v in method_cnt.most_common():
    print(f"  {k!r}: {v}")
print("\nHost_ICTV 分布:")
for k, v in ictv_cnt.most_common():
    print(f"  {k!r}: {v}")

print("\n前 6 条明细:")
for row in rows_bia_plant[:6]:
    keys = ["contig_id", "Family", "Genus", "Host_ICTV", "pred|L1", "pred|L2",
            "Host", "Host_NCBI_lineage", "Final_Host", "Decision_Method"]
    print("  " + " | ".join(f"{k}={row.get(k,'')[:40]}" for k in keys))
