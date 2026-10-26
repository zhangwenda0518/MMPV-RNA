#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""score_mock_benchmark.py — mock 社区三臂真值评分 (在 246 运行)"""
import collections
import os
import subprocess

# benchmark 工作根目录。可用环境变量 MMPV_GP_BENCH_ROOT 覆盖; 默认值与服务器现网路径一致。
ROOT = os.environ.get("MMPV_GP_BENCH_ROOT", "/home/zhangwenda/gp_benchmark")
NAMES = {"NC_001367": "TMV", "NC_011532": "LYMov", "NC_001841": "SPFMV", "NC_001495": "AMV"}
TRUTH_LEN = {"TMV": 6395, "LYMov": 12926, "SPFMV": 10820, "AMV": 3644}
VORDER = ("TMV", "LYMov", "SPFMV", "AMV")
ASM = ROOT + "/asm/mock1/mock1_megahit.contig.fasta"
TSV = ROOT + "/graph/mock1/mock1.graph_verdict.tsv"

verdicts, votes = {}, {}
for line in open(TSV):
    if line.startswith("contig\t"):
        continue
    f = line.rstrip("\n").split("\t")
    verdicts[f[0]] = f[11]
    votes[f[0]] = int(f[5])

cand = {c for c, v in votes.items() if v >= 1}
demoted_all = {c for c, v in verdicts.items() if v == "graph_nonviral"}
derived = cand - demoted_all
gv = {c for c, v in verdicts.items() if v == "graph_viral"}
seeds = {c for c, v in votes.items() if v >= 2}
arms = {"baseline(候选池)": cand, "graph(派生候选)": derived, "graph_full(含0票图救)": cand | gv}

# 相关 contig 全集 vs truth db
allids = cand | gv | demoted_all
with open("/tmp/ids.txt", "w") as fh:
    fh.write("\n".join(sorted(allids)) + "\n")
subprocess.run("seqkit grep -f /tmp/ids.txt %s > /tmp/rel.fa" % ASM,
               shell=True, check=True)
out = subprocess.run(
    "blastn -query /tmp/rel.fa -db %s/score/truthdb -outfmt '6 qseqid sseqid' "
    "-evalue 1e-10 -max_target_seqs 1 -num_threads 16" % ROOT,
    shell=True, capture_output=True, text=True).stdout
truth_of = {}
for line in out.splitlines():
    q, s = line.split("\t")[:2]
    truth_of[q] = next((v for k, v in NAMES.items() if s.startswith(k)), None)

hdr = ("臂", "条数", "真病毒", "背景FP", "精确率") + VORDER
print("%-16s%6s%6s%6s%8s | %9s%9s%9s%9s" % (hdr[0], *hdr[1:]))
for arm, ids in arms.items():
    tp = sum(1 for c in ids if truth_of.get(c))
    fp = len(ids) - tp
    cov = {v: 0 for v in VORDER}
    if ids:
        with open("/tmp/arm_ids.txt", "w") as fh:
            fh.write("\n".join(sorted(ids)) + "\n")
        subprocess.run("seqkit grep -f /tmp/arm_ids.txt %s > /tmp/arm.fa" % ASM,
                       shell=True, check=True)
        o2 = subprocess.run(
            "blastn -query /tmp/arm.fa -db %s/score/truthdb -outfmt "
            "'6 qseqid sseqid qstart qend' -evalue 1e-10 -num_threads 16" % ROOT,
            shell=True, capture_output=True, text=True).stdout
        covered = collections.defaultdict(set)
        for line in o2.splitlines():
            q, s, qs, qe = line.split("\t")
            v = next((v for k, v in NAMES.items() if s.startswith(k)), None)
            if v:
                a, b = min(int(qs), int(qe)), max(int(qs), int(qe))
                covered[v].update(range(max(1, a - 25), b + 25))
        cov = {v: len(covered[v]) for v in VORDER}
    prec = tp / len(ids) if ids else 0.0
    covs = "".join(("%.1f%%" % (100 * cov[v] / TRUTH_LEN[v])).rjust(9) for v in VORDER)
    print("%-16s%6d%6d%6d%8s | %s" % (
        arm, len(ids), tp, fp, "%.1f%%" % (100 * prec), covs))

dem_hit = sum(1 for c in demoted_all if truth_of.get(c))
gv_hit = sum(1 for c in gv if truth_of.get(c))
print()
print("被降级 %d 条(全为0票)中真病毒: %d → TP 误伤 %d" % (len(demoted_all), dem_hit, dem_hit))
print("graph_viral %d 条中真病毒: %d (%.1f%%)" % (len(gv), gv_hit, 100 * gv_hit / max(1, len(gv))))
print("候选池票数: ≥2票 %d (种子), =1票 %d (unknown 交图裁决)" % (len(seeds), len(cand) - len(seeds)))
