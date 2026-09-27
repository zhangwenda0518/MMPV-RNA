import glob
import os
import csv

base = "/home/zhangwenda/data-test/out/06_HostPrediction/C9_ICTV_result"
targets = [
    "ERR2040156_clean_NODE_74_length_3607_cov",
    "ERR2040149_clean_NODE_575_length_2136_co",
]

found = 0
for tsv in glob.glob(os.path.join(base, "*.classified.tsv")):
    if found >= 2:
        break
    with open(tsv) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            cid = row.get("contig_id", "") or row.get("\ufeffcontig_id", "")
            for t in targets:
                if cid.startswith(t) and found < 2:
                    found += 1
                    print(f"[{os.path.basename(tsv)}]")
                    for k, v in row.items():
                        print(f"    {k} = {str(v)[:70]}")
                    print()
print(f"命中 {found} 条")
print("\nC9 结果文件的列名样例:")
first = glob.glob(os.path.join(base, "*.classified.tsv"))[0]
with open(first) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        print("  ", list(row.keys()))
        break
