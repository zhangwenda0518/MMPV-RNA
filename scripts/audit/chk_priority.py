#!/usr/bin/env python3
"""量出「属优先 / 种优先 / 多数派」三种取法各差多少行，以及科属覆盖率

回答两个问题：
  1 科属一层能覆盖多少行（是否真的够用）
  2 Fadolivirus 那类「只有属是异类」的行有多少，排序改一下能救回多少
"""
import csv
import importlib.util
from collections import Counter
from pathlib import Path

PIPE = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline")
VMR = Path("/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")
SRC = Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/"
           "RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated/"
           "final_integrated_classification.tsv")

spec = importlib.util.spec_from_file_location("ann", PIPE / "utils" / "annotate_nucleic_acid.py")
ann = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ann)
idx = ann.load_vmr(VMR)
ORDER = ann.LOOKUP_ORDER

with open(SRC, newline="", encoding="utf-8", errors="replace") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
n = len(rows)
print("输入 %s  %d 行" % (SRC.name, n))
print("当前脚本 LOOKUP_ORDER = %s" % ORDER)

# ---- 1. 各层级覆盖率（该层能否唯一判出大类）----
cov = Counter()
only = Counter()
per_rows = []
for row in rows:
    per = {}
    for rank in ORDER:
        b = ann.resolve_rank(row.get(rank, ""), rank, idx)
        if b:
            per[rank] = b
    per_rows.append(per)
    for r in per:
        cov[r] += 1
    if len(per) == 1:
        only[list(per)[0]] += 1

print()
print("=== 各层级可判行数（该层名字在参考库中唯一确定大类）===")
for r in ORDER:
    print("  %-8s %6d  (%5.2f%%)" % (r, cov[r], 100.0 * cov[r] / n))

gf = sum(1 for p in per_rows if "Genus" in p or "Family" in p)
print()
print("=== 「科属够不够用」===")
print("  Genus 或 Family 至少一层可判 : %d (%5.2f%%)" % (gf, 100.0 * gf / n))
print("  只有 Genus 一层可判          : %d (%5.2f%%)" % (only["Genus"], 100.0 * only["Genus"] / n))
print("  只有 Family 一层可判         : %d (%5.2f%%)" % (only["Family"], 100.0 * only["Family"] / n))
print("  科属都判不出、靠更粗层级      : %d" % sum(1 for p in per_rows if "Genus" not in p and "Family" not in p and p))
print("  整条谱系为空（无从判）        : %d (%5.2f%%)" % (sum(1 for p in per_rows if not p),
                                                  100.0 * sum(1 for p in per_rows if not p) / n))

# ---- 2. 三种取法 ----
ORD_A = ORDER                                     # 现状: Genus 优先
ORD_B = ["Species", "Genus"] + [r for r in ORDER if r not in ("Species", "Genus")]


def pick(per, order):
    for r in order:
        if r in per:
            return per[r]
    return ""


def majority(per):
    if not per:
        return ""
    top = Counter(per.values()).most_common()
    if len(top) > 1 and top[0][1] == top[1][1]:
        return pick(per, ORDER)
    return top[0][0]


res = {"majority": [], "genus_first": [], "species_first": []}
for per in per_rows:
    res["majority"].append(majority(per))
    res["genus_first"].append(pick(per, ORD_A))
    res["species_first"].append(pick(per, ORD_B))

print()
print("=== 三种取法分布 ===")
for k, v in res.items():
    print("  %-14s %s" % (k, dict(Counter(v))))

print()
print("=== 两两分歧 ===")
keys = list(res)
for i in range(len(keys)):
    for j in range(i + 1, len(keys)):
        a, b = keys[i], keys[j]
        d = sum(1 for x, y in zip(res[a], res[b]) if x != y)
        dist = Counter("%s->%s" % (x, y) for x, y in zip(res[a], res[b]) if x != y)
        print("  %-14s vs %-14s : %4d 行 (%5.2f%%)  %s" % (
            a, b, d, 100.0 * d / n, dict(dist)))

print()
print("=== 矛盾行里「一层是异类」的形态（矛盾行总数）===")
contra = [(row, per) for row, per in zip(rows, per_rows) if len(set(per.values())) > 1]
print("  矛盾行 %d" % len(contra))
shape = Counter()
for row, per in contra:
    c = Counter(per.values()).most_common()
    if len(c) > 1 and c[0][1] == c[1][1]:
        shape["平票"] += 1
    else:
        odd = [r for r, b in per.items() if b != c[0][0]]
        shape["异类层数=%d (%s)" % (len(odd), ",".join(sorted(odd)))] += 1
for k, v in shape.most_common(10):
    print("    %-40s %4d" % (k, v))
