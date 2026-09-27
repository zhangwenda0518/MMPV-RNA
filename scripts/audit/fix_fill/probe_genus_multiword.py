#!/usr/bin/env python3
"""统计成品 Genus 列里的多词值, 并查 rankedlineage.dmp 判断该名能否映射到真属。

用法:
  python3 probe_genus_multiword.py <prod.tsv> [<prod.tsv> ...] [--top 12] [--no-ref]

输出: 每个成品的多词 Genus 行数/不同值数、top 值、以及参照库能否给出其真实属/种。
"""
import csv, sys, argparse

RANKED = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"

ap = argparse.ArgumentParser()
ap.add_argument("prods", nargs="+")
ap.add_argument("--top", type=int, default=12)
ap.add_argument("--no-ref", action="store_true")
a = ap.parse_args()


def norm(v):
    if v is None:
        return ""
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    return v


tables = {}
need = set()
for p in a.prods:
    rows = list(csv.DictReader(open(p, newline="", encoding="utf-8", errors="replace"), delimiter="\t"))
    cnt = {}
    for r in rows:
        g = norm(r.get("Genus"))
        if g and g.upper() not in ("NA", "N/A") and " " in g:
            cnt[g] = cnt.get(g, 0) + 1
    tables[p] = (len(rows), cnt)
    need.update(k.lower() for k in cnt)

ref = {}
if not a.no_ref and need:
    with open(RANKED, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = [x.strip().strip("|").strip() for x in line.split("\t|\t")]
            if len(parts) < 4:
                continue
            nm = parts[1].lower()
            if nm in need and nm not in ref:
                ref[nm] = {"species": parts[2], "genus": parts[3] if len(parts) > 3 else ""}
    print("参照库(rankedlineage)命中 name: %d / %d" % (len(ref), len(need)))
    print()

for p in a.prods:
    nrows, cnt = tables[p]
    print("--- %s" % p)
    print("    行数 %d ; Genus 多词行 %d (%.2f%%) / 不同值 %d"
          % (nrows, sum(cnt.values()), 100.0 * sum(cnt.values()) / max(1, nrows), len(cnt)))
    resolv = 0
    for k in cnt:
        e = ref.get(k.lower())
        if e and e["genus"] and e["genus"].upper() not in ("", "NA"):
            resolv += 1
    print("    其中参照库能给出非空 genus 的: %d / %d" % (resolv, len(cnt)))
    for k, v in sorted(cnt.items(), key=lambda kv: -kv[1])[: a.top]:
        e = ref.get(k.lower())
        g = (e or {}).get("genus", "?")
        s = (e or {}).get("species", "?")
        print("      %-44s x%-5d ref: species=%s genus=%s" % (k[:44], v, s, g))
    print()
