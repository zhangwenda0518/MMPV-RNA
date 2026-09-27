#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Independent audit for the ruthenicum v6.7 gate rerun.

Written from scratch to cross-check (NOT transcribe) the executor's outputs.
Semantics are re-derived to match the stated judge:
  For each adjacent rank pair (H high, L low): look up L's typed parent from
  rankedlineage.dmp (columns tax_id|tax_name|species|genus|family|order|class|
  phylum|kingdom|superkingdom). Compare only when L has a definable parent and
  the row's H has a value; unequal => conflict row.

Usage: python3 audit_ruthenicum.py <old.tsv> <new.tsv>
Read-only.
"""
import csv
import sys
from collections import Counter, OrderedDict

RANKED = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"

# rank -> column index in rankedlineage.dmp holding that rank's own name
COL_OF = {
    "Species": 2, "Genus": 3, "Family": 4, "Order": 5,
    "Class": 6, "Phylum": 7, "Kingdom": 8, "Realm": 9,
}
RANKS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

# (high, low) adjacent pairs; parent column of `low` lives in COL_OF[high]
PAIRS = [
    ("Realm", "Kingdom"),
    ("Kingdom", "Phylum"),
    ("Phylum", "Class"),
    ("Class", "Order"),
    ("Order", "Family"),
    ("Family", "Genus"),
    ("Genus", "Species"),
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


def load(path):
    """Return list of normalized row dicts preserving order."""
    rows = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append({k.strip(): norm(v) for k, v in r.items()})
    return rows


def build_ref(names):
    """name(lower) -> dict(rank -> dmp field), first hit wins (like reference)."""
    names = {n.lower() for n in names}
    ref = {}
    with open(RANKED, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
            if len(parts) < 10:
                continue
            nm = parts[1].lower()
            if nm in names and nm not in ref:
                ref[nm] = {L: parts[i] for L, i in COL_OF.items()}
    return ref


def chimera(rows, ref):
    per_pair = Counter()
    conflict_rows = 0
    for r in rows:
        bad = False
        for h, l in PAIRS:
            hv, lv = r.get(h), r.get(l)
            if not hv or not lv:
                continue
            e = ref.get(lv.lower())
            if e is None:
                continue
            pv = (e.get(h) or "").strip()
            if pv and pv != hv:
                per_pair["%s-%s" % (h, l)] += 1
                bad = True
        if bad:
            conflict_rows += 1
    return per_pair, conflict_rows


def main():
    old_p, new_p = sys.argv[1], sys.argv[2]
    old = load(old_p)
    new = load(new_p)

    names = set()
    for rows in (old, new):
        for r in rows:
            for L in RANKS:
                if r.get(L):
                    names.add(r[L])
    ref = build_ref(names)
    print("REF matched names: %d / %d" % (len(ref), len(names)))

    # ---- row counts + contig set comparison ----
    old_ids = [r.get("contig_id") for r in old]
    new_ids = [r.get("contig_id") for r in new]
    so, sn = set(old_ids), set(new_ids)
    print("\n=== ROWS ===")
    print("old data rows : %d (unique ids %d)" % (len(old), len(so)))
    print("new data rows : %d (unique ids %d)" % (len(new), len(sn)))
    print("ids only in OLD: %d" % len(so - sn))
    print("ids only in NEW: %d" % len(sn - so))
    print("id sets identical: %s" % (so == sn))
    print("row ORDER identical: %s" % ([str(x) for x in old_ids] == [str(x) for x in new_ids]))

    # ---- chimera OLD / NEW ----
    print("\n=== CHIMERA: OLD (v6.6f) ===")
    po, co = chimera(old, ref)
    for h, l in PAIRS:
        print("  %-16s %6d" % ("%s-%s" % (h, l), po["%s-%s" % (h, l)]))
    print("  conflict rows %d / %d = %.4f%%" % (co, len(old), 100.0 * co / len(old)))

    print("\n=== CHIMERA: NEW (v6.7) ===")
    pn, cn = chimera(new, ref)
    for h, l in PAIRS:
        print("  %-16s %6d" % ("%s-%s" % (h, l), pn["%s-%s" % (h, l)]))
    print("  conflict rows %d / %d = %.4f%%" % (cn, len(new), 100.0 * cn / len(new)))

    # ---- per-rank value changes (by contig_id) ----
    omap = {str(r.get("contig_id")): r for r in old}
    nmap = {str(r.get("contig_id")): r for r in new}
    common = [k for k in omap if k in nmap]
    chg = Counter()
    filled_to_empty = Counter()
    empty_to_filled = Counter()
    for k in common:
        o, n = omap[k], nmap[k]
        for L in RANKS:
            ov, nv = o.get(L), n.get(L)
            if (ov or "") != (nv or ""):
                chg[L] += 1
            if ov and not nv:
                filled_to_empty[L] += 1
            elif not ov and nv:
                empty_to_filled[L] += 1
    print("\n=== PER-RANK CHANGES (common contigs %d) ===" % len(common))
    for L in RANKS:
        print("  %-9s changed %6d | filled->empty %6d | empty->filled %6d"
              % (L, chg[L], filled_to_empty[L], empty_to_filled[L]))

    # ---- gate-relevant: did the count of members per rank shift? ----
    print("\n=== NON-EMPTY COUNT PER RANK ===")
    for L in RANKS:
        o_c = sum(1 for r in old if r.get(L))
        n_c = sum(1 for r in new if r.get(L))
        print("  %-9s old %6d -> new %6d (%+d)" % (L, o_c, n_c, n_c - o_c))


if __name__ == "__main__":
    main()
