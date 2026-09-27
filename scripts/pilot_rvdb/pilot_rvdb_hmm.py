#!/usr/bin/env python3
"""RVDB HMM pilot — 步骤1: 抽路径②的 1291 条 + Prodigal ORF + 计时探针。"""
import csv
import os
import subprocess
import sys
import time

OUT = "/home/zhangwenda/rvdb_hmm_pilot"
os.makedirs(OUT, exist_ok=True)
S = "/home/zhangwenda/data-test/out/09b_Analysis_Verify/virus_validation/rescue_evidence_scored.tsv"
HQ = "/home/zhangwenda/data-test/out/09_Virome_Analysis/HQ_plant_viruses.fasta"
HMM = "/home/zhangwenda/database/virus-db/RVDB-v31/U-RVDBv31.0-prot.hmm"

target = set()
keep = set()
drop = set()
with open(S) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        v = r["verdict"]
        fr = r.get("fp_reason") or ""
        cid = r["contig_id"]
        if v == "REVIEW" and fr.startswith("no_viral_domain_but_blast"):
            target.add(cid)
        elif v == "KEEP":
            keep.add(cid)
        elif v == "DROP":
            drop.add(cid)
print(f"target(路径②)={len(target)}  keep={len(keep)}  drop={len(drop)}", flush=True)


def extract(ids, path):
    n = 0
    with open(path, "w") as fo:
        with open(HQ) as fi:
            w = False
            for l in fi:
                if l.startswith(">"):
                    cid = l[1:].split()[0]
                    w = cid in ids
                    if w:
                        n += 1
                if w:
                    fo.write(l)
    return n


n = extract(target, f"{OUT}/pilot_1291.fasta")
print(f"抽取 target 序列: {n} 条 → pilot_1291.fasta", flush=True)

# Prodigal ORF
t0 = time.time()
r = subprocess.run(["prodigal", "-p", "meta", "-c", "-i", f"{OUT}/pilot_1291.fasta",
                    "-a", f"{OUT}/proteins.faa", "-o", "/dev/null"],
                   capture_output=True, text=True)
print(f"prodigal rc={r.returncode} ({time.time()-t0:.1f}s)", flush=True)
if r.returncode != 0:
    print(r.stderr[-500:], flush=True)
    sys.exit(1)

# 清理非标准字符 (prodigal -a 可能带 '*')
raw = open(f"{OUT}/proteins.faa").read().split("\n")
clean = []
for l in raw:
    clean.append(l.replace("*", "").replace(" ", "") if not l.startswith(">") else l)
open(f"{OUT}/proteins.clean.faa", "w").write("\n".join(clean))
n_orf = sum(1 for l in clean if l.startswith(">"))
print(f"ORF 数: {n_orf}", flush=True)

# 计时探针: 前 500 个 HMM 搜全部 ORF
import pyhmmer
from pyhmmer.easel import SequenceFile

with SequenceFile(f"{OUT}/proteins.clean.faa", digital=True) as sf:
    seqs = sf.read_block()
print(f"载入 ORF block: {len(seqs)} 条 ({sum(len(s) for s in seqs)} aa)", flush=True)

t1 = time.time()
N = 500
hmms = []
with pyhmmer.plan7.HMMFile(HMM) as hf:
    for i, h in enumerate(hf):
        if i >= N:
            break
        hmms.append(h)
print(f"载入 {len(hmms)} 个 HMM: {time.time()-t1:.1f}s", flush=True)

t2 = time.time()
cnt = 0
nhit = 0
with pyhmmer.plan7.HMMFile(HMM) as hf:
    pass
for hits in pyhmmer.hmmer.hmmsearch(hmms, seqs, cpus=32):
    cnt += 1
    nhit += sum(1 for _ in hits)
dt = time.time() - t2
print(f"搜索 {cnt} 个 HMM 用时 {dt:.1f}s, 命中(hit) {nhit}", flush=True)
if dt > 0:
    print(f"外推全库 13758 个 HMM ≈ {dt*13758/cnt/60:.1f} 分钟", flush=True)
