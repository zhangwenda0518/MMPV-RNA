#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Independent re-check for the v6.7 tax-gate rerun (amarum).
# Criteria mirror /tmp/chimera_multi.py: for each adjacent rank pair (H high, L low),
# look up L in rankedlineage.dmp, take its canonical parent at rank H, compare to the row's H.
# Read-only.
import csv
import sys
from collections import Counter, OrderedDict

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


def read_tsv(path):
    rows = OrderedDict()
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            r = {k.strip(): norm(v) for k, v in r.items()}
            cid = r.get("contig_id")
            rows[cid] = r
    return rows


def load_ref(names):
    need = set(n.lower() for n in names if n)
    ref = {}
    with open(RANKED, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
            if len(parts) < 10:
                continue
            nm = parts[1].lower()
            if nm in need and nm not in ref:
                ref[nm] = {L: parts[REF_IDX[L]] for L in RANK_ORDER}
    return ref


def conflicts(rows, ref):
    bcount = Counter()
    n_conf = 0
    for r in rows.values():
        bad = 0
        for h, l in PAIRS:
            hv, lv = r.get(h), r.get(l)
            if not hv or not lv:
                continue
            e = ref.get(lv.lower())
            if e is None:
                continue
            pv = (e.get(h) or "").strip()
            if pv and pv != hv:
                bcount["%s-%s" % (h, l)] += 1
                bad = 1
        n_conf += bad
    return bcount, n_conf


def report(tag, rows, ref):
    bcount, n_conf = conflicts(rows, ref)
    n = len(rows)
    print("[%s] rows=%d" % (tag, n))
    for h, l in PAIRS:
        print("    %-16s %5d" % ("%s-%s" % (h, l), bcount["%s-%s" % (h, l)]))
    print("    TOTAL conflict_rows=%d rate=%.4f%%" % (n_conf, 100.0 * n_conf / n if n else 0.0))
    return n_conf


def main():
    old_p = sys.argv[1]
    new_p = sys.argv[2] if len(sys.argv) > 2 else None
    old = read_tsv(old_p)
    names = []
    for r in old.values():
        for L in RANK_ORDER:
            if r.get(L):
                names.append(r[L])
    new = None
    if new_p:
        new = read_tsv(new_p)
        for r in new.values():
            for L in RANK_ORDER:
                if r.get(L):
                    names.append(r[L])
    ref = load_ref(names)
    print("== ref entries loaded: %d ==" % len(ref))
    report("OLD " + old_p, old, ref)
    if new is not None:
        report("NEW " + new_p, new, ref)
        so, sn = set(old.keys()), set(new.keys())
        print("== set compare ==")
        print("    old_only=%d new_only=%d common=%d" % (len(so - sn), len(sn - so), len(so & sn)))
        print("    identical_set=%s" % (so == sn))
        chg = Counter()
        lost = Counter()
        for k in (so & sn):
            o, nw = old[k], new[k]
            for L in RANK_ORDER:
                if (o.get(L) or "") != (nw.get(L) or ""):
                    chg[L] += 1
            for L in ("Order", "Family"):
                if o.get(L) and not nw.get(L):
                    lost[L + "_lost"] += 1
                elif not o.get(L) and nw.get(L):
                    lost[L + "_gain"] += 1
        print("== per-rank changed rows (old vs new) ==")
        for L in RANK_ORDER:
            print("    %-9s %6d" % (L, chg[L]))
        print("== Order/Family value-presence change ==")
        print("    Order lost=%d gain=%d" % (lost["Order_lost"], lost["Order_gain"]))
        print("    Family lost=%d gain=%d" % (lost["Family_lost"], lost["Family_gain"]))


if __name__ == "__main__":
    main()
