#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""独立复核脚本（大王的本轮审计，不复用 /tmp 参考实现）。

判据（与线上 chk_chimera.py / chimera_multi.py 对齐）：
  对每对相邻阶元 (H 高, L 低)，用 rankedlineage.dmp 查 L 的定型父级；
  仅当行内 H、L 都有值、且 L 在参照库中有该列的定型父级时才比较；
  pv != hv 即该行在该 pair 上矛盾（每行最多计一次矛盾行）。

dmp 列序: tax_id|tax_name|species|genus|family|order|class|phylum|kingdom|superkingdom
（Realm 用 index 9 的 superkingdom 作为上位参照，与参考实现一致）

用法:
  python3 chimera_audit.py stats <tsv> [<tsv> ...]   # 每个文件单算
  python3 chimera_audit.py cmp  <old.tsv> <new.tsv>  # 两文件对比
只读。
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
PAIRNAME = ["%s-%s" % (h, l) for h, l, _ in PAIRS]


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None
    return v


def read_rows(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return [({k.strip(): norm(v) for k, v in r.items()},
                 (r.get("contig_id") or "").strip().strip('"'))
                for r in csv.DictReader(fh, delimiter="\t")]


def load_ref(names):
    need = {n.lower() for n in names if n}
    ref = {}
    with open(RANKED, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
            if len(parts) < 10:
                continue
            nm = parts[1].lower()
            if nm in need and nm not in ref:
                ref[nm] = {L: parts[REF_IDX[L]] for L in RANKS}
    return ref


def conflict_counts(rows, ref):
    bc = Counter()
    n_rows_conf = 0
    for r, _cid in rows:
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
                bc[h + "-" + l] += 1
                bad = 1
        n_rows_conf += bad
    return bc, n_rows_conf


def cmd_stats(paths):
    rowsets = [read_rows(p) for p in paths]
    names = set()
    for rows in rowsets:
        for r, _ in rows:
            for L in RANKS:
                if r.get(L):
                    names.add(r[L])
    ref = load_ref(names)
    for p, rows in zip(paths, rowsets):
        bc, nc = conflict_counts(rows, ref)
        print("FILE\t%s" % p)
        print("  data_rows\t%d" % len(rows))
        for pn in PAIRNAME:
            print("  pair\t%s\t%d" % (pn, bc[pn]))
        print("  conflict_rows\t%d" % nc)
        print("  conflict_rate\t%.4f%%" % (100.0 * nc / len(rows)))


def cmd_cmp(op, np_):
    old = read_rows(op)
    new = read_rows(np_)
    names = set()
    for rows in (old, new):
        for r, _ in rows:
            for L in RANKS:
                if r.get(L):
                    names.add(r[L])
    ref = load_ref(names)
    for tag, rows in (("OLD", old), ("NEW", new)):
        bc, nc = conflict_counts(rows, ref)
        print("%s\trows=%d\tconflict_rows=%d\trate=%.4f%%" %
              (tag, len(rows), nc, 100.0 * nc / len(rows)))
        for pn in PAIRNAME:
            print("  %s\t%d" % (pn, bc[pn]))
    so = [c for _, c in old]
    sn = [c for _, c in new]
    print("set_old\t%d\tset_new\t%d" % (len(set(so)), len(set(sn))))
    print("only_in_old\t%d\tonly_in_new\t%d" % (len(set(so) - set(sn)), len(set(sn) - set(so))))
    o = {c: r for r, c in old}
    n = {c: r for r, c in new}
    keys = [k for k in o if k in n]
    print("common\t%d" % len(keys))
    chg = Counter()
    lost = Counter()
    for k in keys:
        for L in RANKS:
            if (o[k].get(L) or "") != (n[k].get(L) or ""):
                chg[L] += 1
        for L in ("Order", "Family"):
            if o[k].get(L) and not n[k].get(L):
                lost[L + "_lost"] += 1
            if not o[k].get(L) and n[k].get(L):
                lost[L + "_gain"] += 1
    for L in RANKS:
        print("changed\t%s\t%d" % (L, chg[L]))
    for k in ("Order_lost", "Family_lost", "Order_gain", "Family_gain"):
        print("%s\t%d" % (k, lost[k]))


if __name__ == "__main__":
    if sys.argv[1] == "stats":
        cmd_stats(sys.argv[2:])
    else:
        cmd_cmp(sys.argv[2], sys.argv[3])
