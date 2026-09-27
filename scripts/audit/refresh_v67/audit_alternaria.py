#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""独立复核：Alternaria 门槛 v6.7 重跑产物 vs v6.6f 旧产物。

作者独立实现，仅对齐 /tmp/chimera_multi.py 的判据语义，不复制其代码：
  - 相邻阶元对 (H 高, L 低)：用 rankedlineage.dmp 查 L 的定型父级(H 列)，与行内 H 比较。
  - 仅当 行内 H 有值 且 行内 L 有值 且 dmp 命中 L 且 该命中 H 列有值 时才可比；不等即矛盾。
  - 空/NA/引号/尾随 '*' 归一为缺失。
额外做一次「dmp 同名取最后一次」的敏感性对照，确认判据对同名歧义不敏感。
"""
import csv
import sys
from collections import Counter

RANKED = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
RANKS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
# dmp 列序: 0 tax_id, 1 tax_name, 2 species, 3 genus, 4 family, 5 order, 6 class, 7 phylum, 8 kingdom, 9 superkingdom
COL = {"Species": 2, "Genus": 3, "Family": 4, "Order": 5,
       "Class": 6, "Phylum": 7, "Kingdom": 8, "Realm": 9}
PAIRS = [("Realm", "Kingdom"), ("Kingdom", "Phylum"), ("Phylum", "Class"),
         ("Class", "Order"), ("Order", "Family"), ("Family", "Genus"),
         ("Genus", "Species")]

OLD = sys.argv[1]
NEW = sys.argv[2]


def norm(v):
    if v is None:
        return None
    s = v.strip().strip('"').strip()
    while s.endswith("*"):
        s = s[:-1].strip()
    if not s or s.upper() in ("NA", "N/A"):
        return None
    return s


def load(path):
    rows = {}
    order = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for rec in rd:
            clean = {}
            for k, v in rec.items():
                clean[k.strip()] = norm(v)
            cid = (rec.get("contig_id") or "").strip().strip('"')
            order.append(cid)
            rows[cid] = clean
    return order, rows


def build_ref(names):
    first, last = {}, {}
    with open(RANKED, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
            if len(parts) < 10:
                continue
            nm = parts[1].lower()
            if nm not in names:
                continue
            entry = {r: parts[COL[r]] for r in RANKS}
            if nm not in first:
                first[nm] = entry
            last[nm] = entry
    return first, last


def conflicts(rows, ref, tags):
    pair_cnt = Counter()
    line_cnt = 0
    compared = 0
    for cid, r in rows.items():
        bad = 0
        for h, l in PAIRS:
            hv, lv = r.get(h), r.get(l)
            if not hv or not lv:
                continue
            e = ref.get(lv.lower())
            if e is None:
                continue
            pv = (e.get(h) or "").strip()
            if not pv:
                continue
            compared += 1
            if pv != hv:
                pair_cnt[tags[h] + "-" + tags[l]] += 1
                bad = 1
        line_cnt += bad
    return pair_cnt, line_cnt, compared


old_order, old_rows = load(OLD)
new_order, new_rows = load(NEW)

print("=== 1) 行数与 contig_id 集合 ===")
print("old 数据行=%d  new 数据行=%d" % (len(old_order), len(new_order)))
so, sn = set(old_order), set(new_order)
print("old 唯一 id=%d  new 唯一 id=%d" % (len(so), len(sn)))
print("仅 old 有: %d %s" % (len(so - sn), sorted(so - sn)[:5]))
print("仅 new 有: %d %s" % (len(sn - so), sorted(sn - so)[:5]))
print("集合完全相同: %s" % (so == sn))

print()
print("=== 2) contig_id 顺序是否一致 ===")
print("顺序相同: %s" % (old_order == new_order))

print()
print("=== 3) 相邻阶元矛盾（独立实现） ===")
need = set()
for src in (old_rows, new_rows):
    for r in src.values():
        for L in RANKS:
            if r.get(L):
                need.add(r[L].lower())
ref_first, ref_last = build_ref(need)
print("需查参照 name=%d  dmp 命中=%d" % (len(need), len(ref_first)))

for tag, rows in (("v66f(旧)", old_rows), ("online(新)", new_rows)):
    pc, lines, cmp_n = conflicts(rows, ref_first, {"Realm": "Realm", "Kingdom": "Kingdom",
                                                  "Phylum": "Phylum", "Class": "Class",
                                                  "Order": "Order", "Family": "Family",
                                                  "Genus": "Genus", "Species": "Species"})
    total = sum(pc.values())
    rowcount = len(rows)
    print("[%s] 可比格子=%d 矛盾行=%d 矛盾率=%.2f%%" % (tag, cmp_n, lines, 100.0 * lines / rowcount))
    for h, l in PAIRS:
        print("    %-16s %d" % (h + "-" + l, pc[h + "-" + l]))
    print("    矛盾计数总和=%d" % total)

pc2f, lf, _ = conflicts(old_rows, ref_first, {r: r for r in RANKS})
pc2l, ll, _ = conflicts(old_rows, ref_last, {r: r for r in RANKS})
print("敏感性(dmp 同名取最后) 旧产物 矛盾行 first=%d last=%d 同否=%s" % (lf, ll, lf == ll))

print()
print("=== 4) 各阶元取值变化行数（old vs new） ===")
chg = Counter()
lost = Counter()
for cid in old_order:
    if cid not in new_rows:
        continue
    o, n = old_rows[cid], new_rows[cid]
    for L in RANKS:
        if (o.get(L) or "") != (n.get(L) or ""):
            chg[L] += 1
    for L in ("Order", "Family"):
        if o.get(L) and not n.get(L):
            lost[L] += 1
        if not o.get(L) and n.get(L):
            lost[L + "_gain"] += 1
for L in RANKS:
    print("  %-9s %6d" % (L, chg[L]))
print("  Order 有值->空=%d  空->有值=%d" % (lost["Order"], lost["Order_gain"]))
print("  Family 有值->空=%d  空->有值=%d" % (lost["Family"], lost["Family_gain"]))

print()
print("=== 5) Phylum<->Class 参照相容性 ===")
need2 = set()
for src in (old_rows, new_rows):
    for r in src.values():
        for L in ("Phylum", "Class"):
            if r.get(L):
                need2.add(r[L].lower())
rf2, _ = build_ref(need2)
st = Counter()
for cid in old_order:
    if cid not in new_rows:
        continue

    def bad(r):
        ph, cl = r.get("Phylum"), r.get("Class")
        if not (ph and cl):
            return False
        e = rf2.get(cl.lower())
        if e is None:
            return False
        pv = (e.get("Phylum") or "").strip()
        return bool(pv) and pv != ph
    bo, bn = bad(old_rows[cid]), bad(new_rows[cid])
    st["both_ok" if not bo and not bn else "fixed" if bo and not bn else
       "broken" if bn and not bo else "both_bad"] += 1
for k in ("both_ok", "fixed", "broken", "both_bad"):
    print("  %-10s %6d" % (k, st[k]))
