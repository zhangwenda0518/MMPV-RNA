import csv
import collections
import glob
import os

base_dir = "/home/zhangwenda/data-test/out/06_HostPrediction/C9_ICTV_result"

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

print("=== Mimiviridae 判 Plant: 各列一致性 ===")
for row in rows[:8]:
    print(f"  {row.get('contig_id','')[:38]}")
    print(f"    primary_tool={row.get('primary_tool','')!r} conf={row.get('confidence','')}")
    print(f"    Realm={row.get('Realm','')!r} Kingdom={row.get('Kingdom','')!r}")
    print(f"    Phylum={row.get('Phylum','')!r} Class={row.get('Class','')!r} Order={row.get('Order','')!r}")
    print(f"    Family={row.get('Family','')!r} Genus={row.get('Genus','')!r} Species={row.get('Species','')!r}")
    print(f"    Realm_agree={row.get('Realm_agree','')!r}")
    print(f"    Family_agree={row.get('Family_agree','')!r}")
    print(f"    Genus_agree={row.get('Genus_agree','')!r}")
    print()

# 统计
print("=== 汇总 ===")
print("primary_tool:", dict(collections.Counter(r.get("primary_tool","") for r in rows).most_common(5)))
print("Realm:", dict(collections.Counter(r.get("Realm","") for r in rows).most_common(5)))
print("Kingdom:", dict(collections.Counter(r.get("Kingdom","") for r in rows).most_common(5)))
print("Family_agree top5:", dict(collections.Counter(r.get("Family_agree","") for r in rows).most_common(5)))
