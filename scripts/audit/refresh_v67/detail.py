#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R1-vs-R2 divergence detail + sampler (my own independent script; imports no pipeline code).

Usage: python3 detail.py <rankedlineage.dmp> <tsv1> [tsv2 ...]

R1 = lineage of the first dmp line whose tax_name equals L (rank-agnostic name index)
R2 = set of H values over all dmp lines where L sits in the H-lower rank column
Prints per-file row/pair counts for both, plus classified samples of
R1-only and R2-only contradictions.
"""
import csv
import sys
from collections import Counter, defaultdict

LEVELS = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm"]
DMP_COL = {"Realm": "superkingdom", "Kingdom": "kingdom", "Phylum": "phylum",
           "Class": "class", "Order": "order", "Family": "family",
           "Genus": "genus", "Species": "species"}
COLS = ["tax_id", "tax_name", "species", "genus", "family", "order",
        "class", "phylum", "kingdom", "superkingdom"]


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None
    return v


DMP = sys.argv[1]
TSVS = sys.argv[2:]

r1 = {}
r2 = defaultdict(set)
r2any = set()
rankmap = defaultdict(set)

with open(DMP, encoding="utf-8", errors="ignore") as f:
    for line in f:
        parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
        if len(parts) < 10:
            continue
        rec = dict(zip(COLS, parts))
        nm = rec["tax_name"].lower()
        if nm and nm not in r1:
            r1[nm] = {L: rec[DMP_COL[L]] for L in LEVELS}
        for L in LEVELS:
            if rec[DMP_COL[L]].lower() == nm:
                rankmap[nm].add(L)
        for i in range(len(LEVELS) - 1):
            low, high = LEVELS[i], LEVELS[i + 1]
            lv, hv = rec[DMP_COL[low]], rec[DMP_COL[high]]
            if lv:
                lvv = lv.lower()
                r2any.add((low, lvv))
                if hv:
                    r2[(low, lvv)].add(hv)

pair_order = [(LEVELS[i + 1], LEVELS[i], "%s-%s" % (LEVELS[i + 1], LEVELS[i]))
              for i in range(len(LEVELS) - 1)][::-1]

print("REF r1_names=%d r2_keys=%d r2any=%d names_rankmap=%d" %
      (len(r1), len(r2), len(r2any), len(rankmap)), flush=True)


def reason_r1only(low, lvv, hv, pv, ps):
    if not ps:
        if (low, lvv) not in r2any:
            return "A: name never appears in dmp %s column (ranks=%s)" % (low, sorted(rankmap.get(lvv, [])))
        return "B: all dmp %s-column occurrences have empty %s parent" % (low, low)
    if hv in ps:
        return "C: multi-parent; %s-parent set=%s includes row value, first-occurrence lineage gives %s" % (
            low, sorted(ps)[:4], pv)
    return "D: other (ps=%s)" % sorted(ps)[:4]


def reason_r2only(high, low, lvv, hv, pv, ps):
    e = r1.get(lvv)
    if e is None:
        return "E: name absent from dmp entirely"
    if not pv:
        return "F: first-occurrence lineage has empty %s" % high
    if pv == hv:
        return "G: rank-ambiguous name; %s-column occurrences have %s-parents=%s (no %s)" % (
            low, high, sorted(ps)[:4], hv)
    return "H: other"


for TSV in TSVS:
    total = 0
    r1bad = set()
    r2bad = set()
    pr1 = Counter()
    pr2 = Counter()
    p1o = Counter()
    p2o = Counter()
    reason_ct = Counter()
    r1only_samples = []
    r2only_samples = []
    with open(TSV, encoding="utf-8", errors="ignore") as f:
        rd = csv.DictReader(f, delimiter="\t", quotechar='"')
        for idx, raw in enumerate(rd):
            total += 1
            row = {k.strip(): norm(v) for k, v in raw.items()}
            cid = row.get("contig_id")
            h1 = []
            h2 = []
            for high, low, name in pair_order:
                hv, lv = row.get(high), row.get(low)
                if not hv or not lv:
                    continue
                lvv = lv.lower()
                e = r1.get(lvv)
                pv = (e.get(high) or "").strip() if e is not None else ""
                ps = r2.get((low, lvv))
                c1 = bool(pv) and pv != hv
                c2 = bool(ps) and hv not in ps
                if c1:
                    pr1[name] += 1
                    h1.append((high, low, hv, lv, pv, ps))
                if c2:
                    pr2[name] += 1
                    h2.append((high, low, hv, lv, pv, ps))
                if c1 and not c2:
                    p1o[name] += 1
                    reason_ct["R1only_" + reason_r1only(low, lvv, hv, pv, ps)[0]] += 1
                if c2 and not c1:
                    p2o[name] += 1
                    reason_ct["R2only_" + reason_r2only(high, low, lvv, hv, pv, ps)[0]] += 1
            if h1:
                r1bad.add(idx)
            if h2:
                r2bad.add(idx)
            if h1 and not h2 and len(r1only_samples) < 8:
                r1only_samples.append((cid, h1))
            if h2 and not h1 and len(r2only_samples) < 5:
                r2only_samples.append((cid, h2))

    both = r1bad & r2bad
    r1o = r1bad - r2bad
    r2o = r2bad - r1bad
    print("=" * 70)
    print("FILE %s" % TSV)
    print("rows=%d r1bad=%d r2bad=%d both=%d r1only=%d r2only=%d" %
          (total, len(r1bad), len(r2bad), len(both), len(r1o), len(r2o)))
    print("PAIR  name r1 r2 r1only r2only")
    for high, low, name in pair_order:
        print("  %-15s r1=%6d r2=%6d r1only=%6d r2only=%6d" % (
            name, pr1[name], pr2[name], p1o[name], p2o[name]))
    print("REASONS " + repr(dict(sorted(reason_ct.items()))))
    print("-- R1ONLY SAMPLES (%d cap 8)" % len(r1only_samples))
    for cid, h1 in r1only_samples:
        for high, low, hv, lv, pv, ps in h1:
            print("  cid=%s pair=%s row[%s]=%s row[%s]=%s | r1_first_lineage[%s]=%s | r2_set=%s | %s" %
                  (cid, "%s-%s" % (high, low), high, hv, low, lv, high, pv,
                   sorted(ps)[:4] if ps else None,
                   reason_r1only(low, lv.lower(), hv, pv, ps)))
    print("-- R2ONLY SAMPLES (%d cap 5)" % len(r2only_samples))
    for cid, h2 in r2only_samples:
        for high, low, hv, lv, pv, ps in h2:
            print("  cid=%s pair=%s row[%s]=%s row[%s]=%s | r1_first_lineage[%s]=%s | r2_set=%s | %s" %
                  (cid, "%s-%s" % (high, low), high, hv, low, lv, high, pv,
                   sorted(ps)[:4] if ps else None,
                   reason_r2only(high, low, lv.lower(), hv, pv, ps)))
    sys.stdout.flush()
