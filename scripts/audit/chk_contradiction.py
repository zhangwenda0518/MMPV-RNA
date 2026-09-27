#!/usr/bin/env python3
"""独立量化分类表的谱系自相矛盾比例，并给出 majority / finest 两口径的分歧

不依赖已部署脚本的输出，直接读原始分类表 + VMR 自行复算，用于交叉验证。
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
POS = {r: i for i, r in enumerate(ORDER)}

with open(SRC, newline="", encoding="utf-8", errors="replace") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
print("输入: %s  (%d 行)" % (SRC.name, len(rows)))

n_none = 0
unanimous = 0
contra = 0
ties = 0
maj_counts, fin_counts = Counter(), Counter()
flip = Counter()          # 分歧方向
contra_by_family = Counter()
example = {}

for row in rows:
    per = {}
    for rank in ORDER:
        b = ann.resolve_rank(row.get(rank, ""), rank, idx)
        if b:
            per[rank] = b
    if not per:
        n_none += 1
        continue
    fin_rank = min(per, key=lambda r: POS[r])
    fin = per[fin_rank]
    if len(set(per.values())) == 1:
        unanimous += 1
        maj = fin
    else:
        contra += 1
        contra_by_family[(row.get("Family", "") or "NA").strip()] += 1
        votes = Counter(per.values()).most_common()
        if len(votes) > 1 and votes[0][1] == votes[1][1]:
            ties += 1
            maj = fin
        else:
            maj = votes[0][0]
        if maj != fin and (fin, maj) not in example:
            example[(fin, maj)] = (row, per)
    maj_counts[maj] += 1
    fin_counts[fin] += 1
    if maj != fin:
        flip["%s->%s" % (fin, maj)] += 1

n_dec = len(rows) - n_none
print()
print("=== 谱系内部一致性（只看能唯一判定大类的层级）===")
print("  可判定行        : %d (%.2f%%)" % (n_dec, 100.0 * n_dec / len(rows)))
print("  谱系全空        : %d (%.2f%%)" % (n_none, 100.0 * n_none / len(rows)))
print("  内部一致行      : %d (%.2f/可判定 %.2f%%)" % (
    unanimous, 100.0 * unanimous / len(rows), 100.0 * unanimous / max(n_dec, 1)))
print("  内部矛盾行      : %d (%.2f/可判定 %.2f%%)" % (
    contra, 100.0 * contra / len(rows), 100.0 * contra / max(n_dec, 1)))
print("  其中平票        : %d" % ties)

print()
print("=== 两口径分布 ===")
print("  majority : %s" % dict(maj_counts))
print("  finest   : %s" % dict(fin_counts))
print("  分歧方向 : %s  合计 %d 行 (%.2f%%)" % (
    dict(flip), sum(flip.values()), 100.0 * sum(flip.values()) / len(rows)))

print()
print("=== 矛盾行最多的 Family (前 8) ===")
for fam, c in contra_by_family.most_common(8):
    print("  %-28s %5d" % (fam[:28] or "NA", c))

print()
print("=== 分歧样例（谱系逐层取值）===")
for (fin, maj), (row, per) in list(example.items())[:4]:
    print("  finest=%s / majority=%s" % (fin, maj))
    print("     %s" % " | ".join("%s=%s:%s" % (r, (row.get(r) or "-").strip(), per[r])
                                for r in ORDER if r in per))
