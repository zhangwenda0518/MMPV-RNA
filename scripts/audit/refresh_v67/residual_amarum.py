#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Dump residual conflicts in one integrated table (v6.7 criteria). Read-only.
import csv
import sys

RANKED = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
RANK_ORDER = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
REF_IDX = {"Species": 2, "Genus": 3, "Family": 4, "Order": 5,
           "Class": 6, "Phylum": 7, "Kingdom": 8, "Realm": 9}
PAIRS = [("Realm", "Kingdom"), ("Kingdom", "Phylum"), ("Phylum", "Class"),
         ("Class", "Order"), ("Order", "Family"), ("Family", "Genus"), ("Genus", "Species")]


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None
    return v


path = sys.argv[1]
rows = []
need = set()
with open(path, newline="", encoding="utf-8", errors="replace") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        r = {k.strip(): norm(v) for k, v in r.items()}
        rows.append(r)
for r in rows:
    for L in RANK_ORDER:
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
            ref[nm] = {L: parts[REF_IDX[L]] for L in RANK_ORDER}
for r in rows:
    hits = []
    for h, l in PAIRS:
        hv, lv = r.get(h), r.get(l)
        if not hv or not lv:
            continue
        e = ref.get(lv.lower())
        if e is None:
            continue
        pv = (e.get(h) or "").strip()
        if pv and pv != hv:
            hits.append("%s=%s vs ref %s=%s" % (h, hv, h, pv))
    if hits:
        print("%s\t%s" % (r.get("contig_id"), " ; ".join(hits)))
