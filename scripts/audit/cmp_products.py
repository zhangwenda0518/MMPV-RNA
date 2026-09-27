#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""只读对账：对任意数量的 final_integrated_classification.tsv 做
  (1) VMR 口径 科种 / 科属 不匹配率
  (2) 占位值残留清单（Family/Genus/Species）
  (3) 各阶元非空格数
并支持两两逐格 diff（changed / 变空 / 由空变有）。
用法: python3 cmp_products.py A.tsv [B.tsv ...]
"""
import csv, io, sys, os
from collections import defaultdict, Counter

VMR = "/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv"
TAX = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
QUAL = ["environmental", "uncultured", "unclassified", "unidentified",
        "unassigned", "not assigned", "no hit", "unknown", "undefined"]
# 产物里的空值写成字面量 NA（R 端 na="NA"）。不归一化的话「变空 / 由空变有」两列恒为 0，
# 非空格数也会虚高成全部行数。
NA_TOKENS = {"", "na", "n/a", "nan", "null", "none", "undefined"}


def norm(v):
    v = (v or "").strip()
    return "" if v.lower() in NA_TOKENS else v


def load_vmr():
    """(rank, lower(value)) -> set(Family)，只保留该行 Family 非空者"""
    out = defaultdict(set)
    with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
        rd = csv.reader(f, delimiter="\t", quotechar='"')
        hdr = next(rd)
        idx = {n: hdr.index(n) for n in TAX}
        for row in rd:
            if len(row) <= idx["Species"]:
                continue
            chain = [(row[idx[r]] or "").strip() for r in TAX]
            if not chain[5]:
                continue
            for i, rk in enumerate(TAX):
                if chain[i]:
                    out[(rk, chain[i].lower())].add(chain[5])
    return out


def load_product(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        raw = f.read()
    lines = [ln for ln in raw.splitlines() if not ln.startswith("#")]
    rd = csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t")
    rows, order = {}, []
    for r in rd:
        cid = r.get("contig_id")
        if not cid:
            continue
        order.append(cid)
        rows[cid] = {k: norm(r.get(k)) for k in TAX}
    return rows, order


def raw_empty_tokens(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        t = f.read()
    lines = [ln for ln in t.splitlines() if not ln.startswith("#")]
    rd = csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t")
    cnt = Counter()
    for r in rd:
        for k in TAX:
            v = (r.get(k) or "").strip()
            if v.lower() in NA_TOKENS:
                cnt[repr(v)] += 1
    return cnt.most_common(6)


def report(path, vmr):
    rows, order = load_product(path)
    print("\n" + "=" * 78)
    print("FILE %s" % path)
    print("  行数 %d" % len(order))
    cnt = {r: sum(1 for c in order if rows[c][r]) for r in TAX}
    print("  非空格数 " + " ".join("%s=%d" % (r, cnt[r]) for r in TAX))
    print("  空值字面量分布 %s" % raw_empty_tokens(path))

    n_ge = bad_ge = n_sp = bad_sp = 0
    for c in order:
        fam, ge, sp = rows[c]["Family"], rows[c]["Genus"], rows[c]["Species"]
        if fam and ge:
            f = vmr.get(("Genus", ge.lower()))
            if f:
                n_ge += 1
                if fam not in f:
                    bad_ge += 1
        if fam and sp:
            f = vmr.get(("Species", sp.lower()))
            if f:
                n_sp += 1
                if fam not in f:
                    bad_sp += 1

    def pct(a, b):
        return "%5d /%6d (%5.1f%%)" % (a, b, 100.0 * a / b if b else 0.0)
    print("  VMR口径 科属不匹配(Genus侧) %s" % pct(bad_ge, n_ge))
    print("  VMR口径 科种不匹配(Species侧) %s" % pct(bad_sp, n_sp))

    for col in ("Family", "Genus", "Species"):
        hits = defaultdict(int)
        for c in order:
            v = rows[c][col]
            if v and any(q in v.lower() for q in QUAL):
                hits[v] += 1
        print("  占位残留 %s: %d 格" % (col, sum(hits.values())))
        for v, n in sorted(hits.items(), key=lambda kv: -kv[1])[:8]:
            print("      %5d  %s" % (n, v))
    orphan = sum(1 for c in order if rows[c]["Genus"] and not rows[c]["Family"])
    print("  Genus有值但Family空: %d" % orphan)


def diff(pa, pb):
    ra, oa = load_product(pa)
    rb, ob = load_product(pb)
    ids = [c for c in oa if c in rb]
    print("\n" + "=" * 78)
    print("DIFF  A=%s\n      B=%s" % (pa, pb))
    print("  共同 contig %d（A 独有 %d, B 独有 %d）"
          % (len(ids), len(oa) - len(ids), len(ob) - len(ids)))
    print("  %-9s %8s %8s %8s" % ("阶元", "改变", "变空", "由空变有"))
    tot_ch = 0
    for r in TAX:
        ch = em = fl = 0
        for c in ids:
            a, b = ra[c][r], rb[c][r]
            if a != b:
                ch += 1
                if a and not b:
                    em += 1
                if b and not a:
                    fl += 1
        tot_ch += ch
        print("  %-9s %8d %8d %8d" % (r, ch, em, fl))
    print("  合计改变格数 %d" % tot_ch)


if __name__ == "__main__":
    vmr = load_vmr()
    print("VMR 谱系键: %d 个 (rank,value)" % len(vmr))
    paths = sys.argv[1:]
    for p in paths:
        if not os.path.isfile(p):
            print("缺失: %s" % p)
            continue
        report(p, vmr)
    for i in range(len(paths) - 1):
        if os.path.isfile(paths[i]) and os.path.isfile(paths[i + 1]):
            diff(paths[i], paths[i + 1])
