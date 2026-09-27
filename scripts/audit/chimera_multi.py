#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""8 数据集 × (线上成品 / v6.6f 产物) 的参照库逐级相容性对照。

判据与 /tmp/chk_chimera.py 一致：对每对相邻阶元 (H 高, L 低)，用参照库 rankedlineage.dmp
查 L 的定型父级，与行内 H 比较，不等即矛盾。

用法: python3 chimera_multi.py <manifest.tsv>
manifest 每行: 数据集标签 \t 线上成品.tsv \t v66f成品.tsv
只读，不改任何文件。
"""
import csv
import os
import sys
from collections import Counter, OrderedDict

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


def read_tsv(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        out = []
        for r in rd:
            out.append({k.strip(): norm(v) for k, v in r.items()})
        return out


def main():
    man = []
    with open(sys.argv[1], encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            man.append((parts[0], parts[1], parts[2] if len(parts) > 2 else None))

    # ---- 读全部成品，收集 name ----
    tables = OrderedDict()
    need = set()
    for label, p_old, p_new in man:
        for tag, p in (("online", p_old), ("v66f", p_new)):
            if not p or not os.path.exists(p):
                print("MISSING\t%s\t%s\t%s" % (label, tag, p))
                continue
            rows = read_tsv(p)
            tables[(label, tag)] = rows
            for r in rows:
                for L in RANKS:
                    v = r.get(L)
                    if v:
                        need.add(v.lower())
    print("装载: %d 张表, 需查参照 name %d 个" % (len(tables), len(need)), flush=True)

    ref = {}
    with open(RANKED, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
            if len(parts) < 10:
                continue
            nm = parts[1].lower()
            if nm in need and nm not in ref:
                ref[nm] = {L: parts[REF_IDX[L]] for L in RANKS}
    print("参照命中 name %d 个" % len(ref), flush=True)
    print(flush=True)

    hdr = ["数据集", "口径", "行数"] + ["%s-%s" % (h, l) for h, l, _ in PAIRS] + ["矛盾行", "矛盾率"]
    print("\t".join(hdr))
    for label, p_old, p_new in man:
        for tag in ("online", "v66f"):
            rows = tables.get((label, tag))
            if rows is None:
                continue
            bcount = Counter()
            n_conf = 0
            for r in rows:
                bad = 0
                for h, l, refcol in PAIRS:
                    hv, lv = r.get(h), r.get(l)
                    # 成对判据与 chk_chimera.py 一致：低阶元有参照父级、且行内高阶元有值时才比
                    if not hv or not lv:
                        continue
                    e = ref.get(lv.lower())
                    if e is None:
                        continue
                    pv = (e.get(refcol) or "").strip()
                    if pv and pv != hv:
                        bcount["%s-%s" % (h, l)] += 1
                        bad = 1
                n_conf += bad
            cells = [label, tag, str(len(rows))]
            for h, l, _ in PAIRS:
                cells.append(str(bcount["%s-%s" % (h, l)]))
            cells.append(str(n_conf))
            cells.append("%.2f%%" % (100.0 * n_conf / len(rows)))
            print("\t".join(cells), flush=True)


if __name__ == "__main__":
    main()
