#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""拆 L1（跨阵营）与 L2（同阵营）两层残余矛盾的细节。

L1: 跨阵营矛盾对里，上层阶元的票首占比分布 —— 验证「未过 CASCADE_MIN_SHARE 门槛 -> 不淘汰」这条机制。
L2: 打印若干「同阵营」矛盾行，看是哪一对阶元、参照期望值是什么。

用法: python3 decomp_layers.py <manifest.tsv>
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
        return None, set()
    head, tools = v.split(":", 1)
    s = {t.strip() for t in tools.split(",") if t.strip()}
    try:
        num, den = head.strip().split("/")
        return (int(num), int(den)), s
    except Exception:
        return None, s


man = []
with open(sys.argv[1], encoding="utf-8") as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line.strip() or line.startswith("#"):
            continue
        p = line.split("\t")
        man.append((p[0], p[-1]))

TARGET = {"RNA-Lycium_barbarum_out", "onekp-virus"}
rows_by = {}
need = set()
for label, path in man:
    if label not in TARGET:
        continue
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

for label, rows in rows_by.items():
    frac = Counter()
    n_cross = 0
    same_rows = []
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
            (num, den), sh = parse_agree(r.get(h + "_agree"))
            sl = parse_agree(r.get(l + "_agree"))[1]
            if not sh or not sl:
                kinds.append("na")
                continue
            if sh.isdisjoint(sl):
                n_cross += 1
                frac["%d/%d" % (num, den) if den else "?"] += 1
                kinds.append("cross")
            else:
                kinds.append("same")
                same_rows.append((r, h, l, hv, lv, pv))
    print("=" * 90)
    print("%s：跨阵营矛盾对 %d 个" % (label, n_cross))
    print("  上层阶元票首占比分布（未过 0.5 门槛即不淘汰）:")
    for k, v in frac.most_common(10):
        print("    %-8s %6d  (%.1f%%)" % (k, v, 100.0 * v / n_cross if n_cross else 0))
    tot_ok = sum(v for k, v in frac.items() if k not in ("?",) and int(k.split("/")[0]) / max(1, int(k.split("/")[1])) > 0.5)
    print("  其中上层票首占比 > 0.5（本应触发淘汰）的矛盾对: %d (%.1f%%)" % (tot_ok, 100.0 * tot_ok / n_cross if n_cross else 0))
    print("  同阵营矛盾行数: %d" % len({id(x[0]) for x in same_rows}))
    for r, h, l, hv, lv, pv in same_rows[:6]:
        print("    %s | %s=%s vs %s=%s(参照父级=%s) | %s_agree=%s %s_agree=%s" % (
            r.get("contig_id"), h, hv, l, lv, pv,
            h, r.get(h + "_agree"), l, r.get(l + "_agree")))
