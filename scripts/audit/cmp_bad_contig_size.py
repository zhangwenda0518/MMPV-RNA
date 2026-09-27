#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""对比 barbarum 成品中「Phylum<->Class 参照矛盾行」与全体行的 contig 长度 / 覆盖度分布。

contig_id 形如 CRR1126134_clean_NODE_21_length_1900_cov_0.549535，可解析 length / cov。
用法: python3 cmp_bad_contig_size.py <product.tsv>
"""
import csv
import re
import statistics
import sys
from collections import Counter

RANKED = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
PAT = re.compile(r"_length_(\d+)_cov_([0-9.]+)")


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None
    return v


rows = []
with open(sys.argv[1], newline="", encoding="utf-8", errors="replace") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        rows.append({k.strip(): norm(v) for k, v in r.items()})

need = set()
for r in rows:
    for k in ("Phylum", "Class"):
        if r.get(k):
            need.add(r[k].lower())
ref = {}
with open(RANKED, encoding="utf-8", errors="replace") as fh:
    for line in fh:
        parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
        if len(parts) < 10:
            continue
        nm = parts[1].lower()
        if nm in need and nm not in ref:
            ref[nm] = parts[7]  # phylum


def bad_row(r):
    ph, cl = r.get("Phylum"), r.get("Class")
    if not (ph and cl):
        return False
    rp = ref.get(cl.lower())
    return bool(rp) and rp != ph


def stat(rs, name):
    lens, covs = [], []
    for r in rs:
        m = PAT.search(r.get("contig_id") or "")
        if m:
            lens.append(int(m.group(1)))
            covs.append(float(m.group(2)))
    if not lens:
        print("%s: 无法解析长度" % name)
        return
    lens.sort(); covs.sort()
    q = lambda a, p: a[int(p * (len(a) - 1))]
    print("%-14s n=%5d  length 中位=%6d 四分位=%6d/%6d  cov 中位=%.2f 四分位=%.2f/%.2f  cov<1 占比=%.1f%%" % (
        name, len(lens), statistics.median(lens), q(lens, 0.25), q(lens, 0.75),
        statistics.median(covs), q(covs, 0.25), q(covs, 0.75),
        100.0 * sum(1 for c in covs if c < 1.0) / len(covs)))


allr = [r for r in rows if PAT.search(r.get("contig_id") or "")]
badr = [r for r in allr if bad_row(r)]
print("成品: %s" % sys.argv[1])
stat(allr, "全体行")
stat(badr, "P<->C 矛盾行")
