#!/usr/bin/env python3
"""用表里自带的 *_agree 票数做权重/门槛，对比几种取法

A majority        每层等权，多数派（当前实现）
B agree_weighted  每层的票权 = 认同该层该值的工具数，再比总票权
C gate2_majority  丢弃 k<2 的层，再多数派
D gate2_finest    丢弃 k<2 的层，再取最细层
"""
import csv
import importlib.util
import re
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

AGREE_RE = re.compile(r"^\s*(\d+)\s*/\s*(\d+)\s*:\s*(.*)$")


def parse_agree(s):
    if not s:
        return None
    m = AGREE_RE.match(s)
    if not m:
        return None
    return (int(m.group(1)), int(m.group(2)))


with open(SRC, newline="", encoding="utf-8", errors="replace") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
n = len(rows)

print("=== Genus_agree 无法解析的值长什么样 ===")
c = Counter(r.get("Genus_agree", "") for r in rows if not parse_agree(r.get("Genus_agree", "")))
for k, v in c.most_common(6):
    print("  %-30r %d" % (k[:28], v))

# 每行预先算好: rank -> (class, agree_k)
data = []
for r in rows:
    per = {}
    for rank in ORDER:
        b = ann.resolve_rank(r.get(rank, ""), rank, idx)
        if not b:
            continue
        p = parse_agree(r.get(rank + "_agree", ""))
        k = p[0] if p else 1          # 无法解析时保守当 1 票
        per[rank] = (b, k)
    data.append(per)


def maj(d, gate=None):
    use = {r: v for r, v in d.items() if gate is None or v[1] >= gate}
    if not use:
        return ""
    cnt = Counter(v[0] for v in use.values())
    top = cnt.most_common()
    if len(top) > 1 and top[0][1] == top[1][1]:
        for r in ORDER:
            if r in use:
                return use[r][0]
    return top[0][0]


def agw(d):
    if not d:
        return ""
    w = Counter()
    for cls, k in d.values():
        w[cls] += k
    top = w.most_common()
    if len(top) > 1 and top[0][1] == top[1][1]:
        for r in ORDER:
            if r in d:
                return d[r][0]
    return top[0][0]


def finest(d, gate=None):
    for r in ORDER:
        if r in d and (gate is None or d[r][1] >= gate):
            return d[r][0]
    return ""


variants = {
    "A_majority":       [maj(d) for d in data],
    "B_agree_weighted": [agw(d) for d in data],
    "C_gate2_majority": [maj(d, 2) for d in data],
    "D_gate2_finest":   [finest(d, 2) for d in data],
    "E_finest":         [finest(d) for d in data],
}

print()
print("=== 各取法分布 ===")
for k, v in variants.items():
    c = Counter(v)
    na = c.get("", 0)
    print("  %-18s DNA %6d  RNA %6d  未判 %5d (%5.2f%%)" % (
        k, c.get("DNA", 0), c.get("RNA", 0), na, 100.0 * na / n))

print()
print("=== 与 A_majority 的分歧 ===")
base = variants["A_majority"]
for k, v in variants.items():
    if k == "A_majority":
        continue
    d = sum(1 for x, y in zip(base, v) if x != y)
    dist = Counter("%s->%s" % (x or "NA", y or "NA") for x, y in zip(base, v) if x != y)
    print("  %-18s %5d 行 (%5.2f%%)  %s" % (k, d, 100.0 * d / n, dict(dist)))

# 在 306 打架行上看各取法
print()
print("=== 306 行科属打架行上，各取法判什么 ===")
clash_idx = [i for i, d in enumerate(data)
             if d.get("Family") and d.get("Genus") and d["Family"][0] != d["Genus"][0]]
print("  打架行 %d" % len(clash_idx))
for k, v in variants.items():
    c = Counter(v[i] for i in clash_idx)
    print("  %-18s DNA %4d  RNA %4d  未判 %3d" % (k, c.get("DNA", 0), c.get("RNA", 0), c.get("", 0)))
