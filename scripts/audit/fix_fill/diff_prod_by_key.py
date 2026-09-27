#!/usr/bin/env python3
"""按键(默认 contig_id)对齐对比两个成品分类表。

用法:
  python3 diff_prod_by_key.py <old.tsv> <new.tsv> [--key contig_id] [--top 12]

输出: 列集合差异、行数、仅旧/仅新键、共有键改动行数、各列改动数、方向(NA->值/值->NA/值->值)、
      值->NA 丢失明细、top 转移。按行 zip 的对比在有插入/删除行时会整体错位, 故此处一律按键对齐。
"""
import sys, csv, argparse

ap = argparse.ArgumentParser()
ap.add_argument("old")
ap.add_argument("new")
ap.add_argument("--key", default="contig_id")
ap.add_argument("--top", type=int, default=12)
ap.add_argument("--quiet", action="store_true")
args = ap.parse_args()

NAV = {"NA", "", "-", "nan", "None"}


def load(path):
    with open(path, newline="") as f:
        rd = csv.reader(f, delimiter="\t")
        cols = next(rd)
        ki = cols.index(args.key) if args.key in cols else 0
        rows, dups = {}, 0
        for r in rd:
            if not r:
                continue
            if len(r) < len(cols):
                r = r + [""] * (len(cols) - len(r))
            k = r[ki]
            if k in rows:
                dups += 1
            rows[k] = r
    return cols, rows, dups


cols_o, ro, dup_o = load(args.old)
cols_n, rn, dup_n = load(args.new)

print("列数: old=%d new=%d 一致=%s" % (len(cols_o), len(cols_n), cols_o == cols_n))
if cols_o != cols_n:
    print("  仅旧列: %s" % [c for c in cols_o if c not in cols_n])
    print("  仅新列: %s" % [c for c in cols_n if c not in cols_o])
print("行数: old=%d new=%d (重复键 old=%d new=%d)" % (len(ro), len(rn), dup_o, dup_n))

only_o = [k for k in ro if k not in rn]
only_n = [k for k in rn if k not in ro]
common = set(ro) & set(rn) if len(ro) < len(rn) else set(rn) & set(ro)
print("键: 仅旧=%d 仅新=%d 共有=%d" % (len(only_o), len(only_n), len(common)))
for k in only_o[:6]:
    print("   - 消失: %s" % k)
for k in only_n[:6]:
    print("   + 新增: %s" % k)

ci = {c: i for i, c in enumerate(cols_o)}
per, trans = {}, {}
kind = {"na2val": 0, "val2na": 0, "chg": 0}
loss, gain = {}, {}
changed = 0
for k in common:
    a, b = ro[k], rn[k]
    d = False
    for c in cols_o:
        if c == args.key:
            continue
        i = ci[c]
        va = (a[i] if i < len(a) else "").strip()
        j = cols_n.index(c) if c in cols_n else -1
        vb = (b[j] if j >= 0 and j < len(b) else "").strip()
        if va != vb:
            d = True
            per[c] = per.get(c, 0) + 1
            trans[(c, va, vb)] = trans.get((c, va, vb), 0) + 1
            if va in NAV and vb not in NAV:
                kind["na2val"] += 1
                gain[(c, vb)] = gain.get((c, vb), 0) + 1
            elif va not in NAV and vb in NAV:
                kind["val2na"] += 1
                loss[(c, va)] = loss.get((c, va), 0) + 1
            else:
                kind["chg"] += 1
    if d:
        changed += 1

print("共有键中有改动行 %d (%.2f%%)" % (changed, 100.0 * changed / max(1, len(common))))
print("单元格方向: NA->值(新增)=%d, 值->NA(丢失)=%d, 值->不同值=%d"
      % (kind["na2val"], kind["val2na"], kind["chg"]))
print("各列改动数: " + ", ".join("%s=%d" % (c, per[c]) for c in cols_o if per.get(c)))
if loss:
    print("!! 值->NA 明细 (可能回退):")
    for (c, va), n in sorted(loss.items(), key=lambda kv: -kv[1])[:10]:
        print("   %-12s %-34s -> NA x%d" % (c, va[:34], n))
if not args.quiet:
    print("top 转移:")
    for (c, va, vb), n in sorted(trans.items(), key=lambda kv: -kv[1])[:args.top]:
        print("   %-14s %-30s -> %-30s x%d" % (c, va[:30], vb[:30], n))
