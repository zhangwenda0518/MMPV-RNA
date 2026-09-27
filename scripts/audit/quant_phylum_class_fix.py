#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""量化 barbarum 成品里 Phylum<->Class 矛盾行的成因构成，并预估两种补丁的覆盖面。

用法: python3 quant_phylum_class_fix.py <combined.tsv> <product.tsv>
输出:
  1) 矛盾行数
  2) 其中「Class 阵营与 Phylum 阵营不相交」的跨阵营翻转行数
  3) 阵营不相交行里，Phylum 阵营工具自己报的 Class 是否与 Phylum 参照相容（可被平票淘汰/参照闸门修好的行数）
  4) Class 阵营工具自己报的 Phylum 是否与参照相容（反向证据：矛盾主要出在哪一侧）
只读。
"""
import csv
import sys
from collections import Counter, defaultdict

RANKED = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"


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


def read_tsv(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        return [{k.strip(): norm(v) for k, v in r.items()} for r in rd]


combined_p, product_p = sys.argv[1], sys.argv[2]
prod = read_tsv(product_p)
comb = read_tsv(combined_p)

need = set()
for r in prod:
    for k in ("Phylum", "Class"):
        if r.get(k):
            need.add(r[k].lower())
    for t in parse_agree(r.get("Phylum_agree")) | parse_agree(r.get("Class_agree")):
        pass
for r in comb:
    for k in ("Phylum", "Class"):
        if r.get(k):
            need.add(r[k].lower())

ref = {}
with open(RANKED, encoding="utf-8", errors="replace") as fh:
    for line in fh:
        parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
        if len(parts) < 10:
            continue
        nm = parts[1].lower()
        if nm in need and nm not in ref:
            ref[nm] = {"class": parts[6], "phylum": parts[7]}


def ref_phylum_of_class(c):
    e = ref.get((c or "").lower())
    return (e or {}).get("phylum") or None


def ref_class_of_phylum(p):
    """不唯一，返回启发式：这个 phylum 在参照库里对应的 class 集合（可能为空）"""
    return None


# 工具原始报告
by_contig = defaultdict(list)
for r in comb:
    c = r.get("seq_name") or r.get("contig_id")
    if c:
        by_contig[c].append(r)

bad = []
for r in prod:
    ph, cl = r.get("Phylum"), r.get("Class")
    if not (ph and cl):
        continue
    rp = ref_phylum_of_class(cl)
    if rp and rp != ph:
        bad.append(r)

print("成品行数: %d" % len(prod))
print("Phylum<->Class 参照矛盾行: %d" % len(bad))

cross = 0
fixable = 0
class_side_ok = 0
no_alt = 0
detail = Counter()
for r in bad:
    pt = parse_agree(r.get("Phylum_agree"))
    ct = parse_agree(r.get("Class_agree"))
    ph = r.get("Phylum")
    cl = r.get("Class")
    cid = r["contig_id"]
    tools = {t.get("tool"): t for t in by_contig.get(cid, [])}
    if pt and ct and pt.isdisjoint(ct):
        cross += 1
        detail["跨阵营(Class 阵营与 Phylum 阵营不相交)"] += 1
        # Phylum 阵营工具的 Class 是否与 Phylum 参照相容
        good = False
        for t in pt:
            tr = tools.get(t)
            if tr is None:
                continue
            c2 = tr.get("Class")
            if c2 and ref_phylum_of_class(c2) == ph:
                good = True
                break
        if good:
            fixable += 1
        else:
            no_alt += 1
    # Class 阵营工具的 Phylum 是否与 Class 的参照父级一致
    ok = False
    for t in ct:
        tr = tools.get(t)
        if tr is None:
            continue
        p2 = tr.get("Phylum")
        if p2 and p2 == ref_phylum_of_class(cl):
            ok = True
            break
    if ok:
        class_side_ok += 1

print("其中 跨阵营翻转（Phylum 与 Class 的支持工具完全不重叠）: %d" % cross)
print("  可得替代 Class（Phylum 阵营工具自己报的 Class 与 Phylum 参照相容）: %d" % fixable)
print("  无替代（Phylum 阵营工具在该 contig 没报 Class 或报的也不相容）: %d" % no_alt)
print("Class 阵营工具的 Phylum 与参照一致（说明 Class 侧自洽、错在粗阶元）: %d" % class_side_ok)
