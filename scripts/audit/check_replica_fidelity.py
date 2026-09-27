#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""基线保真度核验：本地复刻的「现行加权票首」OLD 与盘上产物逐格对比。

同时给出 B2（把留空阶元从 prev 中剔除，不用不可信值约束下游）的指标。
"""
import csv, os, re
from collections import defaultdict
exec(open("/tmp/rehearse_toptobottom.py").read().split("# ---------------- 指标")[0])

SHIPPED = os.path.join(BASE, "final_integrated_classification.tsv")
ship = {}
with open(SHIPPED, newline="", encoding="utf-8", errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        ship[r["contig_id"]] = {k: (r.get(k) or "").strip() or None for k in TAX}
print("产物行数 %d, 复刻 contig %d, 交集 %d"
      % (len(ship), len(OLD), len(set(ship) & set(OLD))))

print("\n=== 复刻 OLD vs 盘上产物 逐格一致率（只看双方都在的 contig） ===")
print("  %-9s %9s %9s %9s" % ("阶元", "双方有值", "一致", "一致率"))
for r in TAX:
    both = same = 0
    for cid in set(ship) & set(OLD):
        a, b = ship[cid].get(r), OLD[cid].get(r)
        if a and b:
            both += 1
            if a.lower() == b.lower():
                same += 1
    print("  %-9s %9d %9d %9s" % (r, both, same, "%.3f" % (same / both) if both else "n/a"))

print("\n=== 产物 vs 复刻 的科属不匹配（同一指标，验证基线） ===")
for tag, tab in (("产物", ship), ("复刻OLD", OLD)):
    f_sp = f_ge = n_sp = n_ge = 0
    for cid in set(ship) & set(OLD):
        row = tab.get(cid) or {}
        fam, sp, ge = row.get("Family"), row.get("Species"), row.get("Genus")
        if fam and sp:
            fs_ = vmr_family_of("Species", sp)
            if fs_:
                n_sp += 1
                if fam not in fs_:
                    f_sp += 1
        if fam and ge:
            fg_ = vmr_family_of("Genus", ge)
            if fg_:
                n_ge += 1
                if fam not in fg_:
                    f_ge += 1
    print("  %-8s Species侧 %4d/%6d   Genus侧 %4d/%6d" % (tag, f_sp, n_sp, f_ge, n_ge))

def summarize_g(pool):
    agg = defaultdict(float)
    fir = {}
    for i, (v, w, t) in enumerate(pool):
        agg[v] += w
        fir.setdefault(v, i)
    return sorted(agg.items(), key=lambda kv: (-kv[1], fir[kv[0]]))[0][0]


# ---- B2: 留空阶元不进 prev ----
B2 = {}
for cid in all_ids:
    B2[cid] = {}
    prev = {}
    for r in TAX:
        recs = cells[cid][r]
        if not recs:
            continue
        pool = [x for x in recs if compatible(r, x[0], prev)]
        if pool:
            v = summarize_g(pool)
        else:
            v = None
        if v:
            B2[cid][r] = v
            prev[r] = v
print("\n=== B2（留空不进 prev）指标 ===")
f_sp = f_ge = n_sp = n_ge = 0
for cid in all_ids:
    row = B2[cid]
    fam, sp, ge = row.get("Family"), row.get("Species"), row.get("Genus")
    if fam and sp:
        fs_ = vmr_family_of("Species", sp)
        if fs_:
            n_sp += 1
            if fam not in fs_:
                f_sp += 1
    if fam and ge:
        fg_ = vmr_family_of("Genus", ge)
        if fg_:
            n_ge += 1
            if fam not in fg_:
                f_ge += 1
print("  Species侧科不一致 %d/%d   Genus侧科属不一致 %d/%d" % (f_sp, n_sp, f_ge, n_ge))
print("  %-9s %8s %8s" % ("阶元", "改变", "变空"))
for r in TAX:
    ch = em = 0
    for cid in all_ids:
        o, n = OLD[cid].get(r), B2[cid].get(r)
        if (o or None) != (n or None):
            ch += 1
            if o and not n:
                em += 1
    print("  %-9s %8d %8d" % (r, ch, em))
