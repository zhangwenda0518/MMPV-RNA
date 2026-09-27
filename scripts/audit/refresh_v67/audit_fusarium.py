#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 独立复核脚本（非转述）：行数/集合、逐级相容性矛盾率、阶元取值变化、Order/Family 有值变空。
# 判据严格对齐 /tmp/chimera_multi.py 与 /tmp/cmp_thr_vs_v66f.py。
import csv
import hashlib
import sys
from collections import Counter

RANKED = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
RANKS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
# dmp 列序: tax_id|tax_name|species|genus|family|order|class|phylum|kingdom|superkingdom
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

OLD = sys.argv[1]
NEW = sys.argv[2]


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None
    return v


def raw_id(v):
    return (v or "").strip().strip('"')


def read(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        rows = []
        for r in rd:
            rows.append((raw_id(r.get("contig_id")),
                         {k.strip(): norm(v) for k, v in r.items()}))
        return rows


old = read(OLD)
new = read(NEW)
old_ids = [i for i, _ in old]
new_ids = [i for i, _ in new]
old_set, new_set = set(old_ids), set(new_ids)

print("== 文件指纹 ==")
print("OLD %s md5=%s rows=%d uniq_id=%d" % (OLD, md5(OLD), len(old), len(old_set)))
print("NEW %s md5=%s rows=%d uniq_id=%d" % (NEW, md5(NEW), len(new), len(new_set)))

print("\n== 行数与集合 ==")
print("行数相等: %s (old %d / new %d)" % (len(old) == len(new), len(old), len(new)))
print("contig_id 集合完全一致: %s" % (old_set == new_set))
print("仅旧有 %d 个: %s" % (len(old_set - new_set), sorted(old_set - new_set)[:10]))
print("仅新有 %d 个: %s" % (len(new_set - old_set), sorted(new_set - old_set)[:10]))

# ---- 参照库 ----
need = set()
for rows in (old, new):
    for _, r in rows:
        for L in RANKS:
            v = r.get(L)
            if v:
                need.add(v.lower())
ref = {}
with open(RANKED, encoding="utf-8", errors="replace") as fh:
    for line in fh:
        parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
        if len(parts) < 10:
            continue
        nm = parts[1].lower()
        if nm in need and nm not in ref:
            ref[nm] = {L: parts[REF_IDX[L]] for L in RANKS}
print("\n== 参照库 ==")
print("需查 name %d, 命中 %d" % (len(need), len(ref)))


def conflicts(rows):
    bcount = Counter()
    n_conf = 0
    for _, r in rows:
        bad = 0
        for h, l, refcol in PAIRS:
            hv, lv = r.get(h), r.get(l)
            if not hv or not lv:
                continue
            e = ref.get(lv.lower())
            if e is None:
                continue
            pv = (e.get(refcol) or "").strip()
            if pv and pv != hv:
                bcount["%s-%s" % (h, l)] += 1
                bad = 1
        n_conf += bad
    return bcount, n_conf


print("\n== 逐级相容性矛盾率（我自己重算）==")
hdr = ["口径", "行数"] + ["%s-%s" % (h, l) for h, l, _ in PAIRS] + ["矛盾行", "矛盾率"]
print("\t".join(hdr))
for tag, rows in (("OLD", old), ("NEW", new)):
    bc, nc = conflicts(rows)
    cells = [tag, str(len(rows))] + [str(bc["%s-%s" % (h, l)]) for h, l, _ in PAIRS]
    cells += [str(nc), "%.2f%%" % (100.0 * nc / len(rows))]
    print("\t".join(cells))

# ---- 阶元取值变化 / Order-Family 有值变空 ----
omap = dict(old)
nmap = dict(new)
common = [k for k in old_ids if k in new_set]
print("\n== 阶元取值变化行数（共同 contig %d）==" % len(common))
for L in RANKS:
    c = sum(1 for k in common if (omap[k].get(L) or "") != (nmap[k].get(L) or ""))
    print("  %-9s %6d" % (L, c))
lost = Counter()
for k in common:
    for L in ("Order", "Family"):
        o = omap[k].get(L)
        n = nmap[k].get(L)
        if o and not n:
            lost[L] += 1
        elif not o and n:
            lost[L + "_gain"] += 1
print("Order/Family 有值->空: Order %d / Family %d" % (lost["Order"], lost["Family"]))
print("Order/Family 空->有值: Order %d / Family %d" % (lost["Order_gain"], lost["Family_gain"]))
