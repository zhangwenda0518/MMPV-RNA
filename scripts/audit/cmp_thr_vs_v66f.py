#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""对比「门槛 0.5（v6.6f）」与「门槛 0.499（半数也淘汰）」两份 barbarum 产物。

输出：
  1) 各阶元取值发生变化的行数
  2) Order / Family 由有值变空（信息损失）的行数
  3) Phylum<->Class 参照相容性：v66f 矛盾 -> 新产物相容（修好）、反向（弄坏）、两边都矛盾
只读。
"""
import csv
import sys
from collections import Counter

RANKED = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
RANKS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None
    return v


def read_tsv(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return {norm(r.get("contig_id")) if False else (r.get("contig_id") or "").strip().strip('"'):
                {k.strip(): norm(v) for k, v in r.items()}
                for r in csv.DictReader(fh, delimiter="\t")}


old_p, new_p = sys.argv[1], sys.argv[2]
old = read_tsv(old_p)
new = read_tsv(new_p)
keys = [k for k in old if k in new]
print("共同 contig 数: %d (old %d / new %d)" % (len(keys), len(old), len(new)))

need = set()
for src in (old, new):
    for r in src.values():
        for L in ("Phylum", "Class"):
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
            ref[nm] = parts[7]

chg = Counter()
lost = Counter()
pc = Counter()
for k in keys:
    o, n = old[k], new[k]
    for L in RANKS:
        if (o.get(L) or "") != (n.get(L) or ""):
            chg[L] += 1
    for L in ("Order", "Family"):
        if o.get(L) and not n.get(L):
            lost[L] += 1
        elif not o.get(L) and n.get(L):
            lost[L + "_gain"] += 1

    def bad(r):
        ph, cl = r.get("Phylum"), r.get("Class")
        if not (ph and cl):
            return False
        rp = ref.get(cl.lower())
        return bool(rp) and rp != ph
    bo, bn = bad(o), bad(n)
    pc["both_ok" if not bo and not bn else
       "fixed" if bo and not bn else
       "broken" if bn and not bo else "both_bad"] += 1

print("\n各阶元取值变化行数:")
for L in RANKS:
    print("  %-9s %6d" % (L, chg[L]))
print("\nOrder/Family 有值->空: Order %d / Family %d" % (lost["Order"], lost["Family"]))
print("Order/Family 空->有值: Order %d / Family %d" % (lost["Order_gain"], lost["Family_gain"]))
print("\nPhylum<->Class 参照相容性:")
for k in ("both_ok", "fixed", "broken", "both_bad"):
    print("  %-10s %6d" % (k, pc[k]))
