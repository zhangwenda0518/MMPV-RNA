#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 v6.6f 成品里的残余级间矛盾按机制分桶。

对每一行、每一对发生矛盾的相邻阶元，用两个阶元的 *_agree 列判断支持工具集合：
  disjoint  -> 「跨阵营」：上下级由不同工具阵营支撑（引擎侧拼接）
  有交集    -> 「同阵营」：同一批工具自己报的上下级在参照库里就不相容（工具/参照侧）
  缺标注    -> 无法判定

用法: python3 decomp_residual.py <manifest.tsv>   # 每行: 标签 \t (忽略) \t v66f成品.tsv
只读。
"""
import csv
import sys
from collections import Counter

RANKED = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
RANKS = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm"]
REF_IDX = {"Species": 2, "Genus": 3, "Family": 4, "Order": 5,
           "Class": 6, "Phylum": 7, "Kingdom": 8, "Realm": 9}
PAIRS = [
    ("Realm", "Kingdom", "Realm"),
    ("Kingdom", "Phylum", "Kingdom"),
    ("Phylum", "Class", "Phylum"),
    ("Class", "Order", "Class"),
    ("Order", "Family", "Order"),
    ("Family", "Genus", "Family"),
    ("Genus", "Species", "Genus"),
]


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None
    return v


def parse_agree(v):
    v = norm(v)
    if not v or ":" not in v:
        return set()
    return {t.strip() for t in v.split(":", 1)[1].split(",") if t.strip()}


man = []
with open(sys.argv[1], encoding="utf-8") as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line.strip() or line.startswith("#"):
            continue
        p = line.split("\t")
        man.append((p[0], p[-1]))

rows_by = {}
need = set()
for label, path in man:
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rows = [{k.strip(): norm(v) for k, v in r.items()} for r in csv.DictReader(fh, delimiter="\t")]
    rows_by[label] = rows
    for r in rows:
        for L in RANKS:
            if r.get(L):
                need.add(r[L].lower())

ref = {}
with open(RANKED, encoding="utf-8", errors="replace") as fh:
    for line in fh:
        parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
        if len(parts) < 10:
            continue
        nm = parts[1].lower()
        if nm in need and nm not in ref:
            ref[nm] = {L: parts[REF_IDX[L]] for L in RANKS}

print("数据集\t行数\t矛盾行\t跨阵营行\t同阵营行\t缺标注行\t跨阵营占比\t主导跨阵营边界")
for label, _p in man:
    rows = rows_by[label]
    n_conf = 0
    b_cross = b_same = b_na = 0
    cross_pairs = Counter()
    for r in rows:
        kinds = []
        for h, l, refcol in PAIRS:
            hv, lv = r.get(h), r.get(l)
            if not hv or not lv:
                continue
            e = ref.get(lv.lower())
            if e is None:
                continue
            pv = (e.get(refcol) or "").strip()
            if not pv or pv == hv:
                continue
            sh = parse_agree(r.get(h + "_agree"))
            sl = parse_agree(r.get(l + "_agree"))
            if not sh or not sl:
                kinds.append("na")
            elif sh.isdisjoint(sl):
                kinds.append("cross")
                cross_pairs["%s-%s" % (h, l)] += 1
            else:
                kinds.append("same")
        if not kinds:
            continue
        n_conf += 1
        if "cross" in kinds:
            b_cross += 1
        elif "same" in kinds:
            b_same += 1
        else:
            b_na += 1
    tot_pairs = sum(cross_pairs.values())
    print("%s\t%d\t%d\t%d\t%d\t%d\t%.1f%%\t%s (%d/%d)" % (
        label, len(rows), n_conf, b_cross, b_same, b_na,
        100.0 * b_cross / n_conf if n_conf else 0.0,
        cross_pairs.most_common(1)[0][0] if cross_pairs else "-",
        cross_pairs.most_common(1)[0][1] if cross_pairs else 0, tot_pairs))
