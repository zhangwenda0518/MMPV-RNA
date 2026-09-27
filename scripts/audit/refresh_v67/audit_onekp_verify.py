#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""onekp 门槛变更（v6.6f -> v6.7）独立复核脚本。

独立实现，不调用 /tmp 下任何辅助脚本。判据：
  对每对相邻阶元 (H 高, L 低)，用 rankedlineage.dmp 查 L 的定型父级（dmp 列序
  tax_id|tax_name|species|genus|family|order|class|phylum|kingdom|superkingdom），
  低阶元有参照父级、且行内高阶元有值时才比，不等即矛盾。
  （病毒谱系中 Realm 落在 superkingdom 列，King 落在 kingdom 列，已实证。）

用法: python3 audit_onekp_verify.py <old.tsv> <new.tsv>
只读，不改任何文件。
"""
import csv
import sys
import time
from collections import Counter

RANKED = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
RANKS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
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
IDX = {L: i for i, L in enumerate(RANKS)}


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None
    return v


def iter_rows(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            d = {}
            for k, v in r.items():
                if k is None:
                    continue
                d[k.strip()] = norm(v)
            yield d


def build_ref(needed):
    """name(lower) -> {rank: 参照值}；同名取文件首现（与既有口径一致），并统计重名。"""
    ref = {}
    dup = Counter()
    with open(RANKED, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
            if len(parts) < 10:
                continue
            nm = parts[1].lower()
            if nm not in needed:
                continue
            if nm in ref:
                dup[nm] += 1
                continue
            ref[nm] = {L: parts[REF_IDX[L]] for L in RANKS}
    return ref, dup


def conflict_stats(path, ref):
    """返回 (行数, 非空contig数, 七对分布Counter, 矛盾行数)。"""
    bcount = Counter()
    n_conf = 0
    n = 0
    nonempty = 0
    for d in iter_rows(path):
        n += 1
        if d.get("contig_id"):
            nonempty += 1
        bad = 0
        for h, l, refcol in PAIRS:
            hv, lv = d.get(h), d.get(l)
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
    return n, nonempty, bcount, n_conf


def main():
    old_p, new_p = sys.argv[1], sys.argv[2]
    t0 = time.time()

    # ---- 收集需要的参照 name ----
    needed = set()
    for p in (old_p, new_p):
        for d in iter_rows(p):
            for L in RANKS:
                v = d.get(L)
                if v:
                    needed.add(v.lower())
    print("[%s] 需查参照 name 数: %d" % (time.strftime("%H:%M:%S"), len(needed)), flush=True)
    ref, dup = build_ref(needed)
    print("[%s] 参照命中 name 数: %d；同名多条(取首现)的 name 数: %d"
          % (time.strftime("%H:%M:%S"), len(ref), len(dup)), flush=True)

    # ---- 矛盾率（old / new）----
    res = {}
    for tag, p in (("old_v66f", old_p), ("new_v67", new_p)):
        n, ne, bc, nconf = conflict_stats(p, ref)
        res[tag] = (n, ne, bc, nconf)
        print("[%s] %s: 行=%d 非空contig=%d 矛盾行=%d 矛盾率=%.4f%%"
              % (time.strftime("%H:%M:%S"), tag, n, ne, nconf, 100.0 * nconf / n), flush=True)

    print("\n=== 七对相邻阶元矛盾分布 ===")
    hdr = ["口径", "行数"] + ["%s-%s" % (h, l) for h, l, _ in PAIRS] + ["矛盾行", "矛盾率"]
    print("\t".join(hdr))
    for tag in ("old_v66f", "new_v67"):
        n, ne, bc, nconf = res[tag]
        row = [tag, str(n)] + [str(bc["%s-%s" % (h, l)]) for h, l, _ in PAIRS] \
              + [str(nconf), "%.4f%%" % (100.0 * nconf / n)]
        print("\t".join(row))

    # ---- contig_id 集合 + 逐阶元变化 + Order/Family 得失 ----
    old_ids = set()
    old_vals = {}
    dup_old = 0
    for d in iter_rows(old_p):
        cid = d.get("contig_id") or ""
        if cid in old_vals:
            dup_old += 1
        old_ids.add(cid)
        old_vals[cid] = tuple(d.get(L) or "" for L in RANKS)
    new_ids = set()
    dup_new = 0
    chg = Counter()
    lost = Counter()
    missing = 0
    for d in iter_rows(new_p):
        cid = d.get("contig_id") or ""
        if cid in new_ids:
            dup_new += 1
        new_ids.add(cid)
        o = old_vals.get(cid)
        if o is None:
            missing += 1
            continue
        nv = tuple(d.get(L) or "" for L in RANKS)
        for i, L in enumerate(RANKS):
            if o[i] != nv[i]:
                chg[L] += 1
        for L in ("Order", "Family"):
            ov, nvv = o[IDX[L]], nv[IDX[L]]
            if ov and not nvv:
                lost[L] += 1
            elif not ov and nvv:
                lost[L + "_gain"] += 1

    print("\n=== contig_id 集合 ===")
    print("old 行数=%d 唯一id=%d 重复=%d | new 行数=%d 唯一id=%d 重复=%d"
          % (len(old_vals) + dup_old, len(old_ids), dup_old,
             len(new_ids) + dup_new, len(new_ids), dup_new))
    print("仅 old 有: %d | 仅 new 有: %d | new 中匹配不到 old 的行: %d"
          % (len(old_ids - new_ids), len(new_ids - old_ids), missing))
    if old_ids - new_ids:
        print("  仅old样例:", sorted(old_ids - new_ids)[:5])
    if new_ids - old_ids:
        print("  仅new样例:", sorted(new_ids - old_ids)[:5])

    print("\n=== 各阶元取值变化行数 (old -> new) ===")
    for L in RANKS:
        print("  %-9s %7d" % (L, chg[L]))
    print("=== Order/Family 信息损失 ===")
    print("  有值->空: Order %d / Family %d" % (lost["Order"], lost["Family"]))
    print("  空->有值: Order %d / Family %d" % (lost["Order_gain"], lost["Family_gain"]))

    print("\n[完成] 用时 %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
