#!/usr/bin/env python3
"""列出成品表中与参照库(rankedlineage.dmp)相邻阶元不相容的行 (判据同 /tmp/chimera_multi.py)。

用法:
  python3 chimera_rows.py <prod.tsv> [--max 40] [--keys-only k1,k2]

输出: 每个矛盾行的 contig_id / primary_tool 与逐对矛盾明细(行内高阶元 vs 参照库给出的低阶元父级)。
"""
import csv, sys, os, argparse

RANKED = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
RANKS = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm"]
REF_IDX = {"Species": 2, "Genus": 3, "Family": 4, "Order": 5,
           "Class": 6, "Phylum": 7, "Kingdom": 8, "Realm": 9}
PAIRS = [("Realm", "Kingdom", "Realm"), ("Kingdom", "Phylum", "Kingdom"),
         ("Phylum", "Class", "Phylum"), ("Class", "Order", "Class"),
         ("Order", "Family", "Order"), ("Family", "Genus", "Family"),
         ("Genus", "Species", "Genus")]

ap = argparse.ArgumentParser()
ap.add_argument("prod")
ap.add_argument("--max", type=int, default=40)
ap.add_argument("--keys-only", default="")
a = ap.parse_args()


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None
    return v


with open(a.prod, newline="", encoding="utf-8", errors="replace") as fh:
    rows = [{k.strip(): norm(v) for k, v in r.items()} for r in csv.DictReader(fh, delimiter="\t")]

need = set()
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

only = set(t for t in a.keys_only.split(",") if t)
n_conf = 0
shown = 0
pairc = {}
for r in rows:
    cid = (r.get("contig_id") or "").strip('"')
    confs = []
    for h, l, refcol in PAIRS:
        hv, lv = r.get(h), r.get(l)
        if not hv or not lv:
            continue
        e = ref.get(lv.lower())
        if e is None:
            continue
        pv = (e.get(refcol) or "").strip()
        if pv and pv != hv:
            confs.append((h, l, hv, lv, pv))
            pairc["%s-%s" % (h, l)] = pairc.get("%s-%s" % (h, l), 0) + 1
    if confs:
        n_conf += 1
        if only and cid not in only:
            continue
        if shown < a.max:
            print("%s [%s]" % (cid or "?", r.get("primary_tool") or "?"))
            for h, l, hv, lv, pv in confs:
                print("    %s-%s: 行内 %s=%s vs 参照库 %s 父级 %s" % (h, l, h, hv, lv, pv))
            shown += 1

print("矛盾行 %d / %d 行 (%.2f%%)" % (n_conf, len(rows), 100.0 * n_conf / max(1, len(rows))))
print("逐对计数: " + ", ".join("%s=%d" % (k, v) for k, v in sorted(pairc.items())))
