import glob
import os
import csv
import collections

base = "/home/zhangwenda/data-test/out/06_HostPrediction/C9_ICTV_result"

# Partitiviridae / Mimiviridae / Metaviridae 的判定级别分布
fam_targets = {"Partitiviridae", "Mimiviridae", "Metaviridae"}
stats = {f: collections.Counter() for f in fam_targets}
genus_stats = {f: collections.Counter() for f in fam_targets}

for tsv in glob.glob(os.path.join(base, "*.classified.tsv")):
    with open(tsv) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            fam = (row.get("Family", "") or "").strip()
            if fam not in fam_targets:
                continue
            d = (row.get("Determination_Level", "") or "").split("(")[0].strip()
            ph = (row.get("Predicted_Host", "") or "").strip()
            stats[fam][d] += 1
            if ph == "Plant":
                genus_stats[fam][(row.get("Genus", "") or "").strip()] += 1

for fam in fam_targets:
    tot = sum(stats[fam].values())
    print(f"\n=== {fam} (总 {tot}) ===")
    print("  Determination_Level:")
    for k, v in stats[fam].most_common():
        print(f"    {k}: {v} ({v/tot*100:.1f}%)")
    if genus_stats[fam]:
        print("  Plant 宿主里的属 top8:")
        for k, v in genus_stats[fam].most_common(8):
            print(f"    {k!r}: {v}")
