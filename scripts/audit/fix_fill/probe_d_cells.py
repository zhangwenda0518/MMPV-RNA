#!/usr/bin/env python3
"""把 explain_lost_cells.py 的 D 桶(未解释)逐格摊开, 用 combined 真实票面判定归属。

用法:
  python3 probe_d_cells.py --comb C.tsv --prod-old O.tsv --prod-new N.tsv --ids a,b,c
每格输出: 该 contig 在 combined 各阶元的取值 x 计数, 以及旧/新成品该行取值。
"""
import csv, argparse

RANKS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
NAV = {"", "na", "n/a", "-", "no rank", "undefined", "unknown", "null", "default", "unclassified"}
PH = ("environmental", "uncultured", "unclassified", "unidentified", "unassigned",
      "not assigned", "no hit", "unknown", "undefined", "sp.")

ap = argparse.ArgumentParser()
ap.add_argument("--comb", required=True)
ap.add_argument("--prod-old", required=True)
ap.add_argument("--prod-new", required=True)
ap.add_argument("--ids", required=True)
a = ap.parse_args()
ids = [x.strip() for x in a.ids.split(",") if x.strip()]
idset = set(ids)


def norm(v):
    return (v or "").strip().strip('"')


votes = {k: {r: {} for r in RANKS} for k in idset}
with open(a.comb, newline="", encoding="utf-8", errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        k = norm(r.get("seq_name"))
        if k not in idset:
            continue
        for rk in RANKS:
            v = norm(r.get(rk))
            if v:
                votes[k][rk][v] = votes[k][rk].get(v, 0) + 1


def load(p):
    out = {}
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            k = norm(r.get("contig_id"))
            if k in idset:
                out[k] = r
    return out


old, new = load(a.prod_old), load(a.prod_new)
for k in ids:
    print("=" * 70)
    print("CONTIG %s" % k)
    ro, rn = old.get(k, {}), new.get(k, {})
    for rk in RANKS:
        vo, vn = norm(ro.get(rk)), norm(rn.get(rk))
        if vo == vn:
            continue
        vc = votes.get(k, {}).get(rk, {})
        top = sorted(vc.items(), key=lambda kv: -kv[1])[:5]
        valid = [(v, n) for v, n in top if v.lower() not in NAV]
        ph = [p for p, _ in valid if any(t in p.lower() for t in PH)]
        print("  %-8s old=%-34s new=%-14s | vote=%s | 占位候选=%s"
              % (rk, vo[:34] or "NA", vn[:14] or "NA",
                 "; ".join("%s x%d" % (v[:30], n) for v, n in top), ph or "-"))
