#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""summarize_graph_run.py — 246 实跑结果跨样本汇总 + RVDB 验证"""
import collections
import os
import random
import subprocess
import sys

# 项目根 (246 实跑结果 + 影子组装目录都在其下)。
# 可用环境变量 MMPV_GP_PROJECT_ROOT 覆盖; 默认值与服务器现网路径一致。
_PROJECT_ROOT = os.environ.get("MMPV_GP_PROJECT_ROOT", "/home/zhangwenda/N20220441-06-01")
ROOT = _PROJECT_ROOT + "/02a_GraphProp_test"
SHADOW = _PROJECT_ROOT + "/graph_asm_shadow"
# 病毒数据库根。可用环境变量 MMPV_DB_ROOT 覆盖; 默认值与服务器现网路径一致。
_DB_ROOT = os.environ.get("MMPV_DB_ROOT", "/home/zhangwenda/database/virus-db")
RVDB = _DB_ROOT + "/RVDB-v31/RVDB_nucl"
WORK = "/tmp/gp_summary"
cols = ["contig", "source_tool", "length", "gc", "depth", "votes", "seed_label",
        "community_id", "community_size", "community_viral", "community_nonviral",
        "verdict", "reason"]

tot = collections.Counter()
per = {}
for s in sorted(os.listdir(ROOT)):
    tsv = os.path.join(ROOT, s, s + ".graph_verdict.tsv")
    if not os.path.isfile(tsv):
        continue
    c = collections.Counter()
    gv_ids, nonviral_ids = [], []
    for line in open(tsv):
        if line.startswith("contig\t"):
            continue
        f = line.rstrip("\n").split("\t")
        d = dict(zip(cols, f))
        v = int(d["votes"])
        c["contigs"] += 1
        if v >= 1:
            c["cand"] += 1
            c["cand2" if v >= 2 else "cand1"] += 1
        c[d["verdict"]] += 1
        if d["verdict"] == "graph_nonviral" and v >= 1:
            c["demoted_in_cand"] += 1
        if d["verdict"] == "graph_viral" and v >= 1:
            c["gv_in_cand"] += 1
        if d["verdict"] == "graph_viral" and v == 0:
            c["gv_out"] += 1
            gv_ids.append(d["contig"])
        if d["verdict"] in ("nonviral_seed", "graph_nonviral"):
            nonviral_ids.append(d["contig"])
    per[s] = (c, gv_ids, nonviral_ids)
    for k in c:
        tot[k] += c[k]

hdr = ("样本", "contigs", "候选池", "≥2票", "=1票", "viral_seed", "graph_viral",
       "nonviral_seed", "graph_nonviral", "conflicting", "unknown",
       "候选被降", "候选被救", "0票图救")
keys = ("contigs", "cand", "cand2", "cand1", "viral_seed", "graph_viral",
        "nonviral_seed", "graph_nonviral", "conflicting", "unknown",
        "demoted_in_cand", "gv_in_cand", "gv_out")
w = (12, 8, 7, 6, 6, 11, 12, 14, 15, 11, 8, 9, 9, 9)
print("".join(h.rjust(x) for h, x in zip(hdr, w)))
for s, (c, _, _) in per.items():
    print(s[:-6].rjust(w[0]) + "".join(str(c[k]).rjust(x) for k, x in zip(keys, w[1:])))
print("-" * sum(w))
print("合计".rjust(w[0]) + "".join(str(tot[k]).rjust(x) for k, x in zip(keys, w[1:])))
print()
print("候选池 FP 削减 = {}/{} = {:.2%}".format(
    tot["demoted_in_cand"], tot["cand"],
    tot["demoted_in_cand"] / max(1, tot["cand"])))
print("0票图救片段(记录未注入) =", tot["gv_out"])

# ── RVDB blastn 验证: graph_viral vs 非病毒对照 ──
os.makedirs(WORK, exist_ok=True)
random.seed(42)


def extract(fa, ids, out):
    want = set(ids)
    n = 0
    with open(fa) as fin, open(out, "w") as fout:
        keep = False
        for line in fin:
            if line.startswith(">"):
                cid = line[1:].split()[0]
                keep = cid in want
                if keep:
                    n += 1
            if keep:
                fout.write(line)
    return n


gv_fa = os.path.join(WORK, "graph_viral.fa")
ctl_fa = os.path.join(WORK, "control_nonviral.fa")
n_gv = n_ctl = 0
for s, (c, gv_ids, nonviral_ids) in per.items():
    if not gv_ids and not nonviral_ids:
        continue
    fa = os.path.join(SHADOW, s, s + "_refineC.contig.fasta")
    if gv_ids:
        n_gv += extract(fa, gv_ids, gv_fa + "." + s)
    if nonviral_ids:
        k = min(len(nonviral_ids), max(1, len(gv_ids)))
        ctl = random.sample(nonviral_ids, k)
        n_ctl += extract(fa, ctl, ctl_fa + "." + s)

# 合并分片
for prefix, parts in ((gv_fa, [gv_fa + "." + s for s in per]),
                      (ctl_fa, [ctl_fa + "." + s for s in per])):
    with open(prefix, "w") as fout:
        for p in parts:
            if os.path.isfile(p):
                fout.write(open(p).read())
                os.remove(p)

def blast_hits(fa):
    out = subprocess.run(
        ["blastn", "-query", fa, "-db", RVDB, "-evalue", "1e-5",
         "-outfmt", "6 qseqid", "-max_target_seqs", "1", "-num_threads", "16"],
        capture_output=True, text=True).stdout
    return len({q for q in out.splitlines() if q})


h_gv, h_ctl = blast_hits(gv_fa), blast_hits(ctl_fa)
print()
print("RVDB blastn 验证 (e<=1e-5):")
print("  graph_viral (图拓扑救回, 0票): {}/{} = {:.1%} 有病毒参考命中".format(
    h_gv, n_gv, h_gv / max(1, n_gv)))
print("  非病毒对照 (nonviral_seed+graph_nonviral 抽样): {}/{} = {:.1%}".format(
    h_ctl, n_ctl, h_ctl / max(1, n_ctl)))
