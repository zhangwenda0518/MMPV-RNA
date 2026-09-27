#!/usr/bin/env python3
"""回查成品表中 值->NA / NA->值 的单元格, 并打印对应 contig 在 combined(新旧) 的逐工具行。

用法:
  python3 probe_lost_cells.py --prod-old <bak>/final_integrated_classification.tsv \
      --prod-new <new>/final_integrated_classification.tsv \
      --comb-old <bak combined> --comb-new <new combined> [--focus Epsilonentomopoxvirus,Alphaorpheovirus]
"""
import csv, argparse


def load_prod(p):
    with open(p, newline="") as f:
        rd = csv.reader(f, delimiter="\t")
        cols = next(rd)
        rows = {}
        for r in rd:
            if r:
                rows[r[0].strip('"')] = r
    return cols, rows


def load_comb(p):
    with open(p, newline="") as f:
        rd = csv.reader(f, delimiter="\t")
        cols = next(rd)
        d = {}
        for r in rd:
            if r:
                d.setdefault(r[0].strip('"'), []).append(r)
    return cols, d


ap = argparse.ArgumentParser()
ap.add_argument("--prod-old", required=True)
ap.add_argument("--prod-new", required=True)
ap.add_argument("--comb-old", required=True)
ap.add_argument("--comb-new", required=True)
ap.add_argument("--focus", default="")
ap.add_argument("--max", type=int, default=6)
a = ap.parse_args()

NAV = {"NA", "", "-", "nan"}
co, ro = load_prod(a.prod_old)
cn, rn = load_prod(a.prod_new)
ci = {c: i for i, c in enumerate(co)}
lost, gained = [], []
for k, x in ro.items():
    y = rn.get(k)
    if not y:
        continue
    for c in co:
        if c not in cn:
            continue
        i, j = ci[c], cn.index(c)
        va = (x[i] if i < len(x) else "").strip()
        vb = (y[j] if j < len(y) else "").strip()
        if va != vb:
            if va not in NAV and vb in NAV:
                lost.append((k, c, va))
            elif va in NAV and vb not in NAV:
                gained.append((k, c, vb))

print("值->NA %d 格 / %d contig ; NA->值 %d 格 / %d contig"
      % (len(lost), len(set(t[0] for t in lost)), len(gained), len(set(t[0] for t in gained))))
for k, c, v in lost:
    print("   LOST   %s  %s: %s -> NA" % (k, c, v))
for k, c, v in gained[:20]:
    print("   GAIN   %s  %s: NA -> %s" % (k, c, v))

cbo, dbo = load_comb(a.comb_old)
cbn, dbn = load_comb(a.comb_new)
print("combined 列: %s" % " | ".join(cbo))
focus = [t for t in a.focus.split(",") if t]

ids = []
for k, c, v in lost:
    if k not in ids:
        ids.append(k)
    if focus and any(f in v for f in focus):
        pass
ids = ids[: a.max]
for k in ids:
    print("=== %s" % k)
    for tag, D in (("OLD", dbo), ("NEW", dbn)):
        for r in D.get(k, []):
            print("   %s %s" % (tag, " | ".join(r)))
    if k in ro:
        print("   PROD_OLD %s" % " | ".join(ro[k][1:]))
    if k in rn:
        print("   PROD_NEW %s" % " | ".join(rn[k][1:]))
