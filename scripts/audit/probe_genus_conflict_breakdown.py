#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
probe_genus_conflict_breakdown.py —— 把「Species 首词 != 最终 Genus」的行按严重度分层

分组口径（统一小写比较；首词须以 virus 结尾且非 viridae/virinae 才视为独立属主张）：
  A 同科不同属  : 首词属在参照表里有定型科，且该科 == 行内 Family -> 科一致、属不一致
  B 跨科        : 首词属的定型科与行内 Family 都不同 -> 连科都对不上（两道闸门都没复核）
  C 行 Family 空: 行内 Family 本身为空，没有可比对的科
  D 首词不在参照表: 参照表查不到该属（多为新属/拼写变体）
  E 带限定词    : Species 里含 sp./cf./strain 等占位限定词（名称本身就不完整，另计）

用法：python3 probe_genus_conflict_breakdown.py <final_integrated_classification.tsv> [...]
"""
import csv
import os
import re
import sys

REF = os.path.expanduser("~/database/taxonomy/genus_family_ref.tsv")
BAD_SUFFIX = re.compile(r"viridae$|virinae$", re.I)
# 注意：限定词里的 "sp." 以句点结尾，句点是非单词字符，后面接 \b 永远匹配不到，
# 必须用 (?=\s|$) 收尾（首轮写 \b 导致占位行统计恒 0）。
QUAL = re.compile(r"(?:^|\s)(sp\.|cf\.|aff\.|strain|isolate|uncultured|unclassified|environmental)(?=\s|$)", re.I)


def load_ref():
    d = {}
    with open(REF, "r", encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = (r.get("Genus") or "").strip().strip('"').lower()
            if not g:
                continue
            fams = set()
            for col in ("NCBI_Family", "VMR_Family"):
                v = (r.get(col) or "").strip().strip('"').lower()
                for piece in v.split(";"):
                    piece = piece.strip()
                    if piece:
                        fams.add(piece)
            d[g] = fams
    return d


def run(path, ref):
    if not os.path.exists(path):
        print("  [跳过] 不存在: %s" % path)
        return
    with open(path, "r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    buckets = {"A 同科不同属": 0, "B 跨科": 0, "C 行 Family 空": 0, "D 首词不在参照表": 0}
    n_qual = 0
    n_agree11 = 0
    ex = {k: [] for k in buckets}
    n_binom = 0
    for r in rows:
        sp = (r.get("Species") or "").strip().strip('"')
        ge = (r.get("Genus") or "").strip().strip('"')
        fa = (r.get("Family") or "").strip().strip('"')
        if " " not in sp:
            continue
        first = sp.split(" ")[0]
        if not first.lower().endswith("virus") or BAD_SUFFIX.search(first):
            continue
        if not ge or ge.upper() in ("NA", "N/A", "-"):
            continue
        n_binom += 1
        if first.lower() == ge.lower():
            continue
        if QUAL.search(sp):
            n_qual += 1
        sa = (r.get("Species_agree") or "").strip().strip('"')
        if sa.startswith("1/1"):
            n_agree11 += 1
        fam_l = fa.lower()
        if fam_l in ("", "na", "n/a", "-"):
            key = "C 行 Family 空"
        elif first.lower() not in ref:
            key = "D 首词不在参照表"
        elif fam_l in ref[first.lower()]:
            key = "A 同科不同属"
        else:
            key = "B 跨科"
        buckets[key] += 1
        if len(ex[key]) < 4:
            ex[key].append("%s | F=%s | G=%s | S=%s | S_agree=%s" %
                           (r.get("contig_id", "")[:46], fa or "NA", ge, sp, sa))
    n_conf = sum(buckets.values())
    print("  %s" % path)
    print("    可比对行(双名法且 Genus 非空) %d ; 冲突 %d 行" % (n_binom, n_conf))
    for k in ("A 同科不同属", "B 跨科", "C 行 Family 空", "D 首词不在参照表"):
        print("      %-14s %5d 行 (%.1f%%)" % (k, buckets[k],
              100.0 * buckets[k] / n_conf if n_conf else 0.0))
        for e in ex[k]:
            print("          %s" % e)
    print("      其中带占位限定词(sp./cf./strain 等) %d 行 ; Species_agree 仍写 1/1 的 %d 行"
          % (n_qual, n_agree11))
    print()


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    ref = load_ref()
    print("参照表属数: %d" % len(ref))
    for p in sys.argv[1:]:
        run(p, ref)
