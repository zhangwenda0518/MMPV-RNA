#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""评估「种级相容性闸门」在不同判据下的影响面。
对每行：Species 若是双名法，取首词当属名，用 genus_family_ref.tsv 查双参照定型科。
分类：
  A dual-agree-mismatch : 两参照一致且都不等于行内 Family（现属级闸门就治这类）
  B dual-disagree-neither: 两参照互相不一致，但行内 Family 与两边都不等（新判据目标）
  C dual-disagree-one   : 两参照不一致且行内 Family 等于其一（分类学版本漂移，应放过）
  D single-side         : 只有一边收录该属
  E 未覆盖              : 种名非双名法，或属名不在参照表
用法: python3 probe_sp_gate.py P1.tsv [P2.tsv ...]
"""
import csv, io, sys, os
from collections import defaultdict, Counter

REF = "/home/zhangwenda/database/taxonomy/genus_family_ref.tsv"
TAX = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
NA_TOKENS = {"", "na", "n/a", "nan", "null", "none", "undefined", "-"}


def norm(v):
    v = (v or "").strip().strip('"').replace("*", "").strip()
    return "" if v.lower() in NA_TOKENS else v


def load_ref():
    out = {}
    with open(REF, newline="", encoding="utf-8", errors="replace") as f:
        rd = csv.DictReader(f, delimiter="\t", quotechar='"')
        for r in rd:
            g = norm(r.get("Genus"))
            if not g:
                continue
            gl = g.lower()
            if gl in out:
                continue
            nf, vf = norm(r.get("NCBI_Family")), norm(r.get("VMR_Family"))
            out[gl] = (nf.lower() or None, vf.lower() or None)
    return out


def load(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        raw = f.read()
    lines = [ln for ln in raw.splitlines() if not ln.startswith("#")]
    return list(csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t"))


ref = load_ref()
print("属参照表条目: %d" % len(ref))

for path in sys.argv[1:]:
    if not os.path.isfile(path):
        print("缺失 %s" % path)
        continue
    rows = load(path)
    st = Counter()
    ex = defaultdict(list)
    for r in rows:
        fam, sp = norm(r.get("Family")), norm(r.get("Species"))
        if not (fam and sp):
            continue
        st["有科有种"] += 1
        if " " not in sp:
            st["E 种名非双名法"] += 1
            continue
        g = sp.split(" ", 1)[0].lower()
        if g not in ref:
            st["E 属名不在参照表"] += 1
            continue
        nf, vf = ref[g]
        fl = fam.lower()
        if nf and vf:
            if nf == vf:
                if fl != nf:
                    st["A dual-agree-mismatch"] += 1
                    ex["A"].append((fam, sp, nf))
                else:
                    st["A ok"] += 1
            else:
                if fl in (nf, vf):
                    st["C dual-disagree-mathes-one(放过)"] += 1
                else:
                    st["B dual-disagree-neither"] += 1
                    ex["B"].append((fam, sp, "%s|%s" % (nf, vf)))
        else:
            one = nf or vf
            if fl != one:
                st["D single-side-mismatch"] += 1
                ex["D"].append((fam, sp, one))
            else:
                st["D ok"] += 1
    print("\n=== %s" % path)
    for k in sorted(st):
        print("  %-34s %6d" % (k, st[k]))
    fixable = st["A dual-agree-mismatch"] + st["B dual-disagree-neither"]
    print("  → 现判据(仅A)可治 %d; 加严判据(A+B)可治 %d" % (st["A dual-agree-mismatch"], fixable))
    for tag in ("A", "B", "C", "D"):
        if ex[tag]:
            print("  [%s] 例:" % tag)
            for fam, sp, rf in ex[tag][:5]:
                print("      %-22s %-34s ref=%s" % (fam[:22], sp[:34], rf))
