import csv
import collections

PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"

auth_fam = set()
with open(PLANT, errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        for p in (row.get("Virus_lineage", "") or "").split(";"):
            p = p.strip()
            if p.endswith("viridae"):
                auth_fam.add(p)

print(f"权威库 Plant.tsv 合法植物科 ({len(auth_fam)}):")
for f_ in sorted(auth_fam):
    print(f"  {f_}")

# C9 中判为 Plant 的科全集
base_dir = "/home/zhangwenda/data-test/out/06_HostPrediction/C9_ICTV_result"
import glob, os
c9_fams = collections.Counter()
for tsv in glob.glob(os.path.join(base_dir, "*.classified.tsv")):
    with open(tsv, errors="replace") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if (row.get("Predicted_Host", "") or "").strip() != "Plant":
                continue
            c9_fams[(row.get("Family", "") or "").strip()] += 1

print(f"\n=== C9 判 Plant 但【权威库完全没有】的科 (整科黑名单候选) ===")
for fam, n in c9_fams.most_common():
    if fam and fam != "NA" and fam not in auth_fam:
        print(f"  {fam}: {n} 条")
