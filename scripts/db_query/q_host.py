import csv
import collections

F = "/home/zhangwenda/data-test/out/06_HostPrediction/ensemble_host_summary.tsv"
targets = [
    "ERR2040118_clean_NODE_150_length_3726_cov_61.0",
    "ERR2040121_clean_NODE_646_length_3062_cov_",
]
hit = 0
sample_host = collections.Counter()
bia_host = collections.Counter()

with open(F) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        cid = row["contig_id"]
        if cid.startswith("ERR2040118_") or cid.startswith("ERR2040121_"):
            sample_host[row.get("Final_Host", "")] += 1
        # Biavirus 属的记录
        if (row.get("Genus", "") or "").strip() == "Biavirus":
            bia_host[row.get("Final_Host", "")] += 1
        for t in targets:
            if cid == t or cid.startswith(t):
                hit += 1
                if hit <= 5:
                    print(f"  {cid[:48]} | Host_ICTV={row.get('Host_ICTV','')} | Final_Host={row.get('Final_Host','')} | Method={row.get('Decision_Method','')}")

print(f"\n精确命中: {hit}")
print("\nERR2040118/ERR2040121 的 Final_Host 分布:")
for k, v in sample_host.most_common(8):
    print(f"  {k!r}: {v}")
print("\n全表中 Genus==Biavirus 的 Final_Host 分布:")
for k, v in bia_host.most_common(8):
    print(f"  {k!r}: {v}")
