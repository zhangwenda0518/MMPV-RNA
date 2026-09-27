#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""定性 v6.6f 成品里 Phylum<->Class 参照矛盾的成因。

对 v66f 成品中 Phylum != ref_parent(Class) 的行：
  1) 统计 agree 列形态（票数/工具名），区分「工具投出来的」与「参照回填的」
  2) 打印 top 组合的 5 个样例：v66f 行 / 线上行 / 该 contig 在 combined 表里的各工具原始报告
用法: python3 diag_phylum_class.py <combined.tsv> <online_product.tsv> <v66f_product.tsv>
"""
import csv
import sys
from collections import Counter, defaultdict

RANKED = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
NEED = {"class": 6, "phylum": 7}


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None
    return v


def read_tsv(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        return [{k.strip(): norm(v) for k, v in r.items()} for r in rd]


combined_p, online_p, v66f_p = sys.argv[1], sys.argv[2], sys.argv[3]
v66f = read_tsv(v66f_p)
online = {r["contig_id"]: r for r in read_tsv(online_p) if r.get("contig_id")}
combined = read_tsv(combined_p)

# 参照: name -> 父级
need = set()
for r in v66f:
    for k in ("Phylum", "Class"):
        v = r.get(k)
        if v:
            need.add(v.lower())
ref = {}
with open(RANKED, encoding="utf-8", errors="replace") as fh:
    for line in fh:
        parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
        if len(parts) < 10:
            continue
        nm = parts[1].lower()
        if nm in need and nm not in ref:
            ref[nm] = {"class": parts[6], "phylum": parts[7]}

bad = []
for r in v66f:
    ph, cl = r.get("Phylum"), r.get("Class")
    if not (ph and cl):
        continue
    e = ref.get(cl.lower())
    if e and e.get("phylum") and e["phylum"] != ph:
        bad.append((r, e["phylum"]))

print("v66f 行数: %d" % len(v66f))
print("Phylum<->Class 参照矛盾行: %d (%.2f%%)" % (len(bad), 100.0 * len(bad) / len(v66f)))
print()

# agree 形态
shape = Counter()
for r, _ in bad:
    pa = r.get("Phylum_agree") or ""
    ca = r.get("Class_agree") or ""
    shape[(bool(pa.split("/")[0] not in ("", "0")), bool(ca.split("/")[0] not in ("", "0")))] += 1
print("agree 形态 (Phylum有票, Class有票) -> 行数:")
for k, v in sorted(shape.items()):
    print("  P有票=%-5s C有票=%-5s %d" % (k[0], k[1], v))
print()

combo = Counter((r.get("Phylum"), r.get("Class"), online.get(r["contig_id"], {}).get("Phylum"),
                 online.get(r["contig_id"], {}).get("Class")) for r, _ in bad)
print("top 组合 (v66f P, v66f C | 线上 P, 线上 C) -> 行数:")
for k, v in combo.most_common(10):
    print("  %-22s %-22s | %-22s %-22s  %d" % (k[0], k[1], k[2], k[3], v))
print()

# 工具原始报告：按 contig 汇总
by_contig = defaultdict(list)
for r in combined:
    c = r.get("seq_name") or r.get("contig_id")
    if c:
        by_contig[c].append(r)

print("=== 样例（各 top 组合取 1 行）===")
seen = set()
for r, refph in bad:
    key = (r.get("Phylum"), r.get("Class"))
    if key in seen:
        continue
    seen.add(key)
    if len(seen) > 4:
        break
    cid = r["contig_id"]
    o = online.get(cid, {})
    print("-" * 100)
    print("contig=%s  primary_tool=%s" % (cid, r.get("primary_tool")))
    print("  v6.6f : P=%s(%s) C=%s(%s) F=%s(%s) G=%s  [参照 C 的父级=%s]" % (
        r.get("Phylum"), r.get("Phylum_agree"), r.get("Class"), r.get("Class_agree"),
        r.get("Family"), r.get("Family_agree"), r.get("Genus"), refph))
    print("  线上  : P=%s(%s) C=%s(%s) F=%s(%s) G=%s" % (
        o.get("Phylum"), o.get("Phylum_agree"), o.get("Class"), o.get("Class_agree"),
        o.get("Family"), o.get("Family_agree"), o.get("Genus")))
    for tr in by_contig.get(cid, []):
        tool = tr.get("tool") or tr.get("Tool") or tr.get("method") or "?"
        print("    tool=%-14s P=%-22s C=%-24s F=%-24s G=%s" % (
            tool, tr.get("Phylum"), tr.get("Class"), tr.get("Family"), tr.get("Genus")))
