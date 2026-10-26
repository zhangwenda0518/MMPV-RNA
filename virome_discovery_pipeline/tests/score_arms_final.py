#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""score_arms_final.py — 恢复主链契约后的四臂真值评分 (在 246 运行)"""
import collections
import os
import subprocess

# benchmark 工作根目录。可用环境变量 MMPV_GP_BENCH_ROOT 覆盖; 默认值与服务器现网路径一致。
ROOT = os.environ.get("MMPV_GP_BENCH_ROOT", "/home/zhangwenda/gp_benchmark")
NAMES = {"NC_001367": "TMV", "NC_011532": "LYMov", "NC_001841": "SPFMV", "NC_001495": "AMV"}
TRUTH_LEN = {"TMV": 6395, "LYMov": 12926, "SPFMV": 10820, "AMV": 3644}
VORDER = ("TMV", "LYMov", "SPFMV", "AMV")

ARMS = [
    ("0.旧主链(原始候选,9-30以来)", ROOT + "/ident/mock1/mock1_virus.all.candidate.fasta"),
    ("1.filter基线(无graph)", ROOT + "/filter_base/mock1/candidates/mock1/mock1_virus.all.candidate.fasta"),
    ("2.filter+graph(新主链镜像)", ROOT + "/filter_graph/mock1/candidates/mock1/mock1_virus.all.candidate.fasta"),
    ("3.graph直通(旧路径,对照)", ROOT + "/graph/candidates/mock1/mock1_virus.all.candidate.fasta"),
]

all_ids = set()
arm_ids = {}
for name, fa in ARMS:
    ids = set()
    cur = None
    for line in open(fa):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            ids.add(cur)
    arm_ids[name] = ids
    all_ids |= ids

with open("/tmp/all_ids.txt", "w") as fh:
    fh.write("\n".join(sorted(all_ids)) + "\n")
subprocess.run("seqkit grep -f /tmp/all_ids.txt %s/asm/mock1/mock1_megahit.contig.fasta > /tmp/all_rel.fa"
               % ROOT, shell=True, check=True)
out = subprocess.run(
    "blastn -query /tmp/all_rel.fa -db %s/score/truthdb -outfmt '6 qseqid sseqid' "
    "-evalue 1e-10 -max_target_seqs 1 -num_threads 16" % ROOT,
    shell=True, capture_output=True, text=True).stdout
truth_of = {}
for line in out.splitlines():
    q, s = line.split("\t")[:2]
    truth_of[q] = next((v for k, v in NAMES.items() if s.startswith(k)), None)

hdr = ("臂", "条数", "真病毒", "背景FP", "精确率") + VORDER
print("%-26s%6s%6s%6s%8s | %9s%9s%9s%9s" % (hdr[0], *hdr[1:]))
summary = {}
for name, _ in ARMS:
    ids = arm_ids[name]
    tp = sum(1 for c in ids if truth_of.get(c))
    fp = len(ids) - tp
    cov = {v: set() for v in VORDER}
    if ids:
        with open("/tmp/arm_ids.txt", "w") as fh:
            fh.write("\n".join(sorted(ids)) + "\n")
        subprocess.run("seqkit grep -f /tmp/arm_ids.txt %s/asm/mock1/mock1_megahit.contig.fasta > /tmp/arm.fa"
                       % ROOT, shell=True, check=True)
        o2 = subprocess.run(
            "blastn -query /tmp/arm.fa -db %s/score/truthdb -outfmt '6 qseqid sseqid qstart qend' "
            "-evalue 1e-10 -num_threads 16" % ROOT,
            shell=True, capture_output=True, text=True).stdout
        for line in o2.splitlines():
            q, s, qs, qe = line.split("\t")
            v = next((v for k, v in NAMES.items() if s.startswith(k)), None)
            if v:
                a, b = min(int(qs), int(qe)), max(int(qs), int(qe))
                cov[v].update(range(max(1, a - 25), b + 25))
    prec = tp / len(ids) if ids else 0.0
    covs = "".join(("%.1f%%" % (100 * len(cov[v]) / TRUTH_LEN[v])).rjust(9) for v in VORDER)
    print("%-26s%6d%6d%6d%8s | %s" % (name, len(ids), tp, fp, "%.1f%%" % (100 * prec), covs))
    summary[name] = (len(ids), tp, fp, prec)

print()
b = summary[ARMS[0][0]]
f1 = summary[ARMS[1][0]]
f2 = summary[ARMS[2][0]]
print("净贡献分解 (真值):")
print("  filter 自身: 候选 %d→%d, FP %d→%d (削 %d), 真病毒 %d→%d (误伤 %d)"
      % (b[0], f1[0], b[2], f1[2], b[2] - f1[2], b[1], f1[1], b[1] - f1[1]))
print("  graph 增量 (filter→filter+graph): 候选 %d→%d, FP %d→%d (削 %d), 真病毒 %d→%d (误伤 %d)"
      % (f1[0], f2[0], f1[2], f2[2], f1[2] - f2[2], f2[1], f2[1], f1[1] - f2[1]))
