"""勘明: 662 涉及的 27 个属在 KEEP 中各有几条, 判断 SDT 重算范围。"""
import csv
from collections import defaultdict

ACVDIR = "/home/zhangwenda/data-test/out/09b_Analysis_Verify"
CLS = f"{ACVDIR}/acvirus_classify/final_result_with_confidence.tsv"
KEEP = f"{ACVDIR}/class_KEEP.fasta"

# 662 涉及的 27 属
killed_genera = set()
with open("/tmp/blacklist_kill_detail.tsv") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        g = (row.get("genus") or "").strip()
        if g and g != "NA":
            killed_genera.add(g)

keep_ids = set()
with open(KEEP) as f:
    for line in f:
        if line.startswith(">"):
            keep_ids.add(line[1:].split()[0])

gen_of = defaultdict(list)
with open(CLS, errors="replace") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        cid = (row.get("Nucleotide", "") or "").strip()
        g = (row.get("Genus", "") or "").strip()
        if cid in keep_ids and g and g != "NA":
            gen_of[g].append(cid)

print("=== 27 个受影响属在当前 KEEP 中的条数 ===")
print(f"{'属':<24} {'KEEP条数':>8}  {'影响':<10}")
print("-" * 50)
need = []
for g in sorted(killed_genera):
    n = len(gen_of.get(g, []))
    hit = "受影响" if g in gen_of else "已清空"
    flag = ""
    if 2 <= n <= 80:
        flag = "→ 需重算 SDT"
        need.append(g)
    elif n > 80:
        flag = "→ 超80上限,跳过"
    elif n == 1:
        flag = "→ 仅1条,不满足"
    else:
        flag = "→ 无序列"
    print(f"{g:<24} {n:>8}  {hit:<10} {flag}")

print(f"\n=== 需重算 SDT 的属 ({len(need)}) ===")
print(", ".join(need) if need else "无")
