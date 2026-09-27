import csv
import collections

CLS = "/home/zhangwenda/data-test/out/09b_Analysis_Verify/acvirus_classify"
c = collections.Counter()
n = 0
with open(f"{CLS}/final_result_with_confidence.tsv", errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        n += 1
        fam = (row.get("Family", "") or "NA").strip()
        c[fam] += 1

print(f"总记录: {n}")
fams = [k for k in c if k != "NA"]
print(f"有科名的科数: {len(fams)}")
print("\n=== 科分布 top20 ===")
for k, v in c.most_common(20):
    print(f"  {k:<26} {v}")

# 与 43 条涉及的科对比
affected = ["Tombusviridae", "Partitiviridae", "Aspiviridae", "Caulimoviridae",
            "Metaviridae", "Geminiviridae", "Potyviridae", "Secoviridae", "Pithoviridae"]
print("\n=== 受影响科是否仍存在 ===")
for fam in affected:
    print(f"  {fam:<22} {c.get(fam, 0)}")
