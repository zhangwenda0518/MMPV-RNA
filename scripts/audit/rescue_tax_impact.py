#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""共识引擎改动对 rescue 的影响：只看 rescue 真正用得到的那一列（Genus）。

rescue_pipeline.py 的三处依赖：
  分支 C  按 Genus 分组（无属的进 no_genus_groups）
  分支 D  genus_len 拯救（无属 → status=no_taxonomy，直接不救）
  报告     skip reason 标注 no_taxonomy
三个都只读 05_Taxonomy/Votus.integrated/final_integrated_classification.tsv 的 Genus 列。

用法: python3 rescue_tax_impact.py <tag> <旧成品.tsv> <新成品.tsv> <rescue_centroids.fasta>
"""
import csv
import sys
from collections import Counter

EMPTY = {"", "na", "n/a", "nan", "-", "none", "null", "unclassified", "unknown"}


def nz(x):
    if x is None:
        return None
    x = x.strip().strip('"').strip()
    return None if x.lower() in EMPTY else x


def load(path, col="Genus"):
    out = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            c = nz(r.get("contig_id"))
            if c:
                out[c] = {k: nz(r.get(k)) for k in
                          ("Realm", "Phylum", "Class", "Order", "Family", "Genus", "Species")}
    return out


def fasta_ids(path):
    ids = set()
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line[1:].strip().split()[0])
    return ids


tag, old_p, new_p, resc = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
old, new = load(old_p), load(new_p)
rid = fasta_ids(resc)

# rescue 候选里有多少条在成品表里能对上
miss = len(rid - set(old))
inold = rid & set(old)

g_old = {c for c in inold if old[c]["Genus"]}
g_new = {c for c in inold if new.get(c, {}).get("Genus")}
lost = g_old - g_new
gain = g_new - g_old

# 丢属的那些还有没有科（有科=只是属被挖掉）
lost_fam = sum(1 for c in lost if old[c]["Family"])

print("### %s" % tag)
print("| 项 | 值 |")
print("|---|---|")
print("| 成品行数（旧 / 新） | %d / %d |" % (len(old), len(new)))
print("| 全表有属（旧 / 新） | %d / %d |" % (
    sum(1 for v in old.values() if v["Genus"]), sum(1 for v in new.values() if v["Genus"])))
print("| rescue 候选条数 | %d |" % len(rid))
print("| 其中能在成品表对上 | %d（对不上 %d，多为 known/过滤阶段外） |" % (len(inold), miss))
print("| rescue 候选里有属（旧 / 新） | %d / %d |" % (len(g_old), len(g_new)))
print("| **丢属** | **%d**（其中仍有科 %d） |" % (len(lost), lost_fam))
print("| 得属 | %d |" % len(gain))

if lost:
    # 丢属的进 rescue 分支时会变成 no_taxonomy（分支 D 不救、分支 C 无属分组）
    print()
    print("丢属样例（前 5 条，含旧科属）:")
    for c in sorted(lost)[:5]:
        print("  %s  Family=%s Genus=%s" % (c, old[c]["Family"] or "NA", old[c]["Genus"]))
print()
