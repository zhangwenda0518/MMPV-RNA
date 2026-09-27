#!/usr/bin/env python3
"""深挖 RVDB HMM pilot: 换阈值/覆盖率后, 目标组与负对照能不能分开。"""
import csv
from collections import Counter, defaultdict

TSV = "/home/zhangwenda/rvdb_hmm_pilot/rvdb_hmm_per_contig.tsv"

rows = []
with open(TSV) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        rows.append(r)
print(f"总行 {len(rows)}")


def fn(x):
    try:
        return float(x)
    except Exception:
        return None


buckets = ["KEEP", "no_viral_domain_but_blast", "short_fragment", "DROP"]
by = defaultdict(list)
for r in rows:
    by[r["bucket"]].append(r)

print("\n=== 多阈值分离度 ===")
thrs = [("E<=1e-2", "e2", 1e-2, None), ("E<=1e-5", "e5", 1e-5, None),
        ("E<=1e-10", "e10", 1e-10, None),
        ("E<=1e-5 & cov>=0.5", "e5c", 1e-5, 0.5),
        ("E<=1e-10 & cov>=0.7", "e10c", 1e-10, 0.7),
        ("score>=50", "s50", None, None), ("score>=100", "s100", None, None)]
hdr = f"{'bucket':<30}{'n':<7}" + "".join(f"{t[0]:<20}" for t in thrs)
print(hdr)
for b in buckets:
    rs = by[b]
    line = f"{b:<30}{len(rs):<7}"
    for _, key, e, c in thrs:
        cnt = 0
        for r in rs:
            ev = fn(r["best_evalue"]) if r["best_evalue"] else None
            cv = fn(r["best_cov"]) if r["best_cov"] else None
            sc = fn(r["best_score"]) if r["best_score"] else None
            ok = False
            if key.startswith("e") and not key.startswith("s"):
                if ev is not None and ev <= e:
                    ok = True
                    if c is not None and (cv is None or cv < c):
                        ok = False
            elif key == "s50":
                ok = sc is not None and sc >= 50
            elif key == "s100":
                ok = sc is not None and sc >= 100
            cnt += ok
        line += f"{100*cnt/len(rs):<20.1f}"
    print(line)

print("\n=== 各 bucket best_hmm Top10 (看是否被个别万能HMM主导) ===")
for b in buckets:
    c = Counter(r["best_hmm"] for r in by[b] if r["best_hmm"])
    print(f"  {b} (有命中 {sum(c.values())}):")
    for k, v in c.most_common(10):
        print(f"      {v:<6}{k}")

print("\n=== 覆盖率分布 (best_cov, 仅命中者) ===")
for b in buckets:
    cvs = sorted(fn(r["best_cov"]) for r in by[b] if r["best_cov"] and fn(r["best_cov"]) is not None)
    if not cvs:
        continue
    n = len(cvs)
    print(f"  {b:<30} n={n:<6} 中位={cvs[n//2]:.2f} "
          f"cov>=0.5 占 {100*sum(1 for x in cvs if x>=0.5)/n:.1f}% "
          f"cov>=0.8 占 {100*sum(1 for x in cvs if x>=0.8)/n:.1f}%")

print("\n=== 命中数 n_hits 分布 ===")
for b in buckets:
    ns = [int(r["n_hits"]) for r in by[b] if r["n_hits"]]
    if not ns:
        continue
    print(f"  {b:<30} 中位={sorted(ns)[len(ns)//2]}  max={max(ns)}  >=3 占 {100*sum(1 for x in ns if x>=3)/len(ns):.1f}%")
