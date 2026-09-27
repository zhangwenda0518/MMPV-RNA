#!/usr/bin/env python3
"""对 值->NA 丢失格做「旧值是否有票」判定, 区分真丢失与 fill 层构造值。

流式扫描新旧 combined, 只收集丢失格对应的 (contig, rank) 键, 内存 O(丢失格数)。

用法:
  python3 probe_lost_backing.py --prod-old O --prod-new N --comb-old CO --comb-new CN [--label X] [--max 30]
判定:
  backed_old : 旧值在【旧】combined 该阶元里出现过(可能来自旧 fill 构造)
  backed_new : 旧值在【新】combined 该阶元里仍出现 = 真被工具投票, 属规则/门槛清空
  novote_new : 新 combined 该阶元一个有效取值都没有(含占位串) = 新值 NA 是票面必然
"""
import csv, argparse

RANKS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
NAV = {"", "na", "n/a", "-", "no rank", "undefined", "unknown", "null", "default", "unclassified"}

ap = argparse.ArgumentParser()
ap.add_argument("--prod-old", required=True)
ap.add_argument("--prod-new", required=True)
ap.add_argument("--comb-old", required=True)
ap.add_argument("--comb-new", required=True)
ap.add_argument("--label", default="")
ap.add_argument("--max", type=int, default=30)
a = ap.parse_args()


def norm(v):
    return (v or "").strip().strip('"')


def load(p):
    out = {}
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out[norm(r.get("contig_id"))] = r
    return out


def norm_rank(stripped):
    for rk in RANKS:
        if stripped == rk:
            return rk
    return None


old, new = load(a.prod_old), load(a.prod_new)
lost = {}          # (contig, rank) -> old value
for k, ro in old.items():
    rn = new.get(k)
    if not rn:
        continue
    for rk in RANKS:
        vo, vn = norm(ro.get(rk)), norm(rn.get(rk))
        if vo and vo.lower() not in NAV and (not vn or vn.lower() in NAV):
            lost[(k, rk)] = vo

co, cn = {}, {}
for path, store in ((a.comb_old, co), (a.comb_new, cn)):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        rd = csv.reader(f, delimiter="\t")
        cols = next(rd)
        ci = {c: i for i, c in enumerate(cols)}
        ki = ci[cols[0]]
        for r in rd:
            if not r:
                continue
            for rk in RANKS:
                key = (norm(r[ki]), rk)
                if key in lost:
                    store.setdefault(key, set()).add(norm(r[ci[rk]]))

backed_old = backed_new = novote_new = 0
detail = {"真丢失(新combined仍有票)": [], "旧值仅旧combined有(fill构造?)": [], "新票面为空": []}
for (k, rk), vo in lost.items():
    o = co.get((k, rk), set())
    n = cn.get((k, rk), set())
    n_valid = {v for v in n if v and v.lower() not in NAV}
    if any(v == vo for v in o):
        backed_old += 1
    if any(v == vo for v in n):
        backed_new += 1
    if not n_valid:
        novote_new += 1
    if any(v == vo for v in n):
        detail["真丢失(新combined仍有票)"].append((k, rk, vo, sorted(n_valid)[:3]))
    elif not n_valid:
        detail["新票面为空"].append((k, rk, vo, sorted(o)[:3]))
    else:
        detail["旧值仅旧combined有(fill构造?)"].append((k, rk, vo, sorted(n_valid)[:3]))

print("== %s 丢失格 %d" % (a.label, len(lost)))
print("   新旧 combined 都投过该值(= 真被工具投过票, 由规则/门槛清空): %d" % backed_new)
print("   仅旧 combined 有该值(新票面里没这个值): %d" % (backed_old - backed_new))
print("   新 combined 该阶元票面为空(NA 是必然): %d" % novote_new)
for tag, rows in detail.items():
    print("   -- %s: %d" % (tag, len(rows)))
    for r in rows[: a.max]:
        print("      %-58s %-8s %-30s | 新票面=%s" % (r[0][:58], r[1], r[2][:30], r[3]))
