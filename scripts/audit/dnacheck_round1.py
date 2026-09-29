#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""分析 PAF：每 query 取最佳比对（按对齐碱基×一致性），给三档判定档位。"""
import csv
from collections import defaultdict

OUT = "/home/zhangwenda/goji_dnacheck_20260929"
GENOMES = ["zhonghua", "ningxia", "heiguo"]

best = {}  # qname -> (genome, pident, qcov, tname, ts, te)
for g in GENOMES:
    with open(f"{OUT}/mm_{g}.paf") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            q, qlen = p[0], int(p[1])
            qs, qe = int(p[2]), int(p[3])
            t, ts, te = p[5], int(p[7]), int(p[8])
            nmatch, alen = int(p[9]), int(p[10])
            de = dict(kv.split(":", 2)[0:2] for kv in [] )  # placeholder
            de_val = None
            for kv in p[12:]:
                parts = kv.split(":")
                if parts[0] == "de":
                    de_val = float(parts[2])
            pid = 1.0 - de_val if de_val is not None else nmatch / alen
            qcov = (qe - qs) / qlen
            score = (qe - qs) * pid
            cur = best.get(q)
            if cur is None or score > cur[0]:
                best[q] = (score, g, pid, qcov, t, ts, te, qlen)

rows = []
with open(f"{OUT}/dna_manifest.tsv") as f:
    rd = csv.DictReader(f, delimiter="\t")
    for r in rd:
        q = f"{r['Cohort']}__{r['contig_id']}"
        b = best.get(q)
        if r["length"] == "MISSING" or b is None:
            r.update(genome_hit="none", pid="", qcov="", hit_locus="")
        else:
            _, g, pid, qcov, t, ts, te, qlen = b
            r["genome_hit"] = g
            r["pid"] = f"{pid*100:.1f}"
            r["qcov"] = f"{qcov*100:.1f}"
            r["hit_locus"] = f"{t}:{ts}-{te}"
        rows.append(r)


def tier(r):
    if r["genome_hit"] == "none":
        return "T3_no_hit"
    pid, qcov = float(r["pid"]), float(r["qcov"])
    if pid >= 90 and qcov >= 80:
        return "T1_recent_EVE"      # 近期整合/同源拷贝，直接定案
    if pid >= 90 and qcov >= 30:
        return "T2_partial_match"
    if pid >= 55 and qcov >= 50:
        return "T2_partial_match"
    return "T3_no_hit" if (pid < 55 or qcov < 30) else "T2_partial_match"


for r in rows:
    r["tier"] = tier(r)

with open(f"{OUT}/genome_match_round1.tsv", "w") as fo:
    fo.write("Cohort\tcontig_id\tGenus\tSpecies\tFamily_true\tCategory\tn_samples\tprevalence_pct\tcheckv\tlength\tgenome_hit\tpid\tqcov\thit_locus\ttier\n")
    for r in rows:
        fo.write("\t".join(str(r.get(c, "")) for c in
                 ["Cohort", "contig_id", "Genus", "Species", "Family_true", "Category", "n_samples", "prevalence_pct", "checkv", "length", "genome_hit", "pid", "qcov", "hit_locus", "tier"]) + "\n")

import collections
cnt = collections.Counter(r["tier"] for r in rows)
print("round1 档位:", dict(cnt))
print()
print("T3_no_hit（待敏感补扫）按队列:")
for k, v in collections.Counter(r["Cohort"] for r in rows if r["tier"] == "T3_no_hit").most_common():
    print(" ", k, v)
print()
print("T3_no_hit 按科:")
for k, v in collections.Counter(r["Family_true"] for r in rows if r["tier"] == "T3_no_hit").most_common():
    print(" ", k, v)
