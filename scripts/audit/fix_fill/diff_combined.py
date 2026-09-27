#!/usr/bin/env python3
# 对比两个 combined_taxonomy.tsv: 逐行逐阶元统计改动
# 用法: python3 /tmp/diff_combined.py old.tsv new.tsv [--top 20]
import argparse

RANKS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

ap = argparse.ArgumentParser()
ap.add_argument("a")
ap.add_argument("b")
ap.add_argument("--top", type=int, default=20)
ap.add_argument("--quiet", action="store_true")
args = ap.parse_args()


def load(p):
    with open(p) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in f if l.strip()]
    return hdr, rows


ha, ra = load(args.a)
hb, rb = load(args.b)
if ha != hb:
    print("!! header 不同: %s vs %s" % (ha, hb))
if len(ra) != len(rb):
    print("!! 行数不同: %d vs %d" % (len(ra), len(rb)))
ci = {r: ha.index(r) for r in RANKS if r in ha}
changed_rows = 0
per = {r: 0 for r in ci}
trans = {}
kind = {"na2val": 0, "val2na": 0, "chg": 0}
loss = {}
loss_keys = []
key_mismatch = 0
NAV = ("NA", "", "-")
for x, y in zip(ra, rb):
    if x[:2] != y[:2]:
        key_mismatch += 1
    d = False
    for r, i in ci.items():
        va = x[i].strip()
        vb = y[i].strip()
        if va != vb:
            d = True
            per[r] += 1
            trans[(r, va, vb)] = trans.get((r, va, vb), 0) + 1
            if va in NAV and vb not in NAV:
                kind["na2val"] += 1
            elif va not in NAV and vb in NAV:
                kind["val2na"] += 1
                loss[(r, va)] = loss.get((r, va), 0) + 1
                if len(loss_keys) < 12:
                    loss_keys.append((x[0], x[1], r, va))
            else:
                kind["chg"] += 1
    if d:
        changed_rows += 1

print("行数 %d (键不一致 %d)" % (len(ra), key_mismatch))
print("有改动行 %d (%.3f%%)" % (changed_rows, 100.0 * changed_rows / max(1, len(ra))))
print("单元格方向: NA->值(回填)=%d, 有值->NA(丢失)=%d, 值->不同值(覆盖)=%d" % (kind["na2val"], kind["val2na"], kind["chg"]))
print("各阶元改动数: " + ", ".join("%s=%d" % (r, per[r]) for r in RANKS if per[r]))
if loss:
    print("!! 有值->NA 明细 (可能回退):")
    for (r, va), c in sorted(loss.items(), key=lambda kv: -kv[1])[:10]:
        print("   %-8s %-32s -> NA x%d" % (r, va[:32], c))
    print("   定位 (前 %d 格):" % len(loss_keys))
    for k0, k1, r, va in loss_keys:
        print("     %-46s [%s] %s: %s -> NA" % (k0[:46], k1, r, va))
print("top 转移 (旧 -> 新):")
for (r, va, vb), c in sorted(trans.items(), key=lambda kv: -kv[1])[:args.top]:
    print("   %-8s %-32s -> %-32s x%d" % (r, va[:32], vb[:32], c))
if not args.quiet and changed_rows:
    print("前若干改动样本:")
    n = 0
    for x, y in zip(ra, rb):
        if x[:2] != y[:2] or any(x[ci[r]] != y[ci[r]] for r in ci):
            print("   %s [%s]" % (x[0][:60], x[1]))
            print("     old: %s" % " | ".join(x[ci[r]] for r in RANKS))
            print("     new: %s" % " | ".join(y[ci[r]] for r in RANKS))
            n += 1
            if n >= 6:
                break
