#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_eval_chain.py — 评测集按新主链顺序跑 filter→graph 并按类评分 (在 246 运行)"""
import collections
import os
import subprocess
import sys

# 仓库根 (filter_virus.py / graph_propagation.py 所在)。
# 可用环境变量 MMPV_REPO_ROOT 覆盖; 默认值与服务器现网路径一致。
REPO_ROOT = os.environ.get("MMPV_REPO_ROOT", "/home/zhangwenda/MMPV-RNA")
# 评测工作根目录。可用环境变量 MMPV_GP_EVAL_ROOT 覆盖; 默认值与服务器现网路径一致。
ROOT = os.environ.get("MMPV_GP_EVAL_ROOT", "/home/zhangwenda/gp_eval")
IDENT = ROOT + "/ident/eval"
ASM = ROOT + "/asm/eval/eval_megahit.contig.fasta"
NAMES = {"positive": "positive", "negative_A": "negA(随机)",
         "negative_B": "negB(保守域陷阱)", "negative_C": "negC(EVE)"}

# ── 1) 重建 ensemble 并集候选 ──
import glob
tool_files = [f for f in glob.glob(IDENT + "/eval_virus.*.result.id")
              if not f.endswith("all.result.id")]
votes = collections.Counter()
for f in tool_files:
    for line in open(f):
        cid = line.strip()
        if cid:
            votes[cid] += 1
union = sorted(votes)
with open(IDENT + "/eval_virus.all.result.id", "w") as fh:
    fh.write("\n".join(union) + "\n")
want = set(union)
n = 0
with open(ASM) as fin, open(IDENT + "/eval_virus.all.candidate.fasta", "w") as fout:
    keep = False
    for line in fin:
        if line.startswith(">"):
            cid = line[1:].split()[0]
            keep = cid in want
            if keep:
                n += 1
        if keep:
            fout.write(line)
hist = collections.Counter(votes[c] for c in union)
print("工具数:", len(tool_files), "候选池:", n, "票数分布:", dict(sorted(hist.items())))

# ── 2) filter (纯过滤; graph 未跑, 直通镜像) ──
os.makedirs(ROOT + "/filter/eval", exist_ok=True)
r = subprocess.run(
    ["python", REPO_ROOT + "/virome_discovery_pipeline/filter_virus.py",
     "-i", IDENT + "/eval_virus.all.candidate.fasta",
     "-o", ROOT + "/filter/eval", "-m", "filter", "-t", "16"],
    capture_output=True, text=True)
print("filter rc:", r.returncode)
print(r.stdout[-600:])

# ── 3) graph (裁决 filter 无证据保留结果, 更新 COBRA 镜像) ──
r = subprocess.run(
    ["python", REPO_ROOT + "/virome_discovery_pipeline/graph_propagation.py",
     "--asm-dir", ROOT + "/asm", "--ident-dir", ROOT + "/ident",
     "--reads-dir", ROOT + "/reads", "--filter-dir", ROOT + "/filter",
     "--output-dir", ROOT + "/graph", "-t", "16", "--force"],
    capture_output=True, text=True)
print("graph rc:", r.returncode)
print(r.stdout[-800:])

# ── 4) 按类评分 ──
idmap = {}
for line in open(ROOT + "/id_map.tsv"):
    f = line.rstrip("\n").split("\t")
    if f[0] != "eval_id":
        idmap[f[0]] = (f[2], int(f[4]) if len(f) > 4 else 0)


def ids_of(path):
    ids = []
    if os.path.isfile(path):
        for line in open(path):
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


ARMS = [
    ("0.旧主链(原始候选)", IDENT + "/eval_virus.all.candidate.fasta"),
    ("1.filter(纯过滤)", ROOT + "/filter/eval/eval_virus.all.candidate_filtered.fasta"),
    ("2.filter+graph(新主链)", ROOT + "/graph/eval/eval_virus.candidate.final.fasta"),
]
cls_of = {c: i for c, i in idmap.items()}

print()
hdr = ("臂",) + tuple(NAMES[c] for c in ("positive", "negative_A",
                                         "negative_B", "negative_C"))
print("%-26s %8s %8s %10s %10s" % hdr)
for name, fa in ARMS:
    ids = ids_of(fa)
    cnt = collections.Counter(cls_of.get(c, "?") for c in ids)
    print("%-26s %8d %8d %10d %10d" % (
        name, len(ids), cnt["positive"], cnt["negative_A"],
        cnt["negative_B"], cnt["negative_C"]))

# 每类的保留率 / 移除数
base_ids = set(ids_of(ARMS[0][1]))
fin_ids = set(ids_of(ARMS[2][1]))
print()
per = collections.Counter()
for c in base_ids:
    per[(cls_of.get(c, "?"), c in fin_ids)] += 1
print("新主链相对原始候选的变动 (按类):")
for cls in ("positive", "negative_A", "negative_B", "negative_C"):
    kept = per[(cls, True)]
    removed = per[(cls, False)]
    tot = kept + removed
    print("  %s: %d 条在候选池 → 保留 %d / 移除 %d" % (NAMES[cls], tot, kept, removed))

# graph 裁决审计
audit = ROOT + "/graph/eval/eval.graph_adjust_calls.tsv"
if os.path.isfile(audit):
    act = collections.Counter()
    for line in open(audit):
        f = line.rstrip("\n").split("\t")
        if f[0] != "contig":
            act[f[4]] += 1
    print("graph 裁决审计:", dict(act))
