#!/usr/bin/env python3
"""零假设: 打乱 ORF 氨基酸顺序(保留组成), 重跑 RVDB HMM, 对比命中率。"""
import csv
import random
import time
from collections import defaultdict

import pyhmmer
from pyhmmer.easel import SequenceFile

OUT = "/home/zhangwenda/rvdb_hmm_pilot"
V3 = f"{OUT}/rvdb_hmm_per_contig_v3.tsv"
HMM = "/home/zhangwenda/database/virus-db/RVDB-v31/U-RVDBv31.0-prot.hmm"
clean_path = f"{OUT}/all_proteins.clean.faa"
shuf_path = f"{OUT}/all_proteins.shuffled.faa"

random.seed(42)
n = 0
with open(shuf_path, "w") as fo, open(clean_path) as fi:
    for l in fi:
        if l.startswith(">"):
            fo.write(l)
            n += 1
        else:
            s = list(l.strip())
            random.shuffle(s)
            fo.write("".join(s) + "\n")
print(f"打乱 {n} 条 ORF", flush=True)

meta = {}
with open(V3) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        meta[r["contig_id"]] = r["bucket"]

with SequenceFile(shuf_path, digital=True) as sf:
    seqs = sf.read_block()

hit = defaultdict(int)
tot = defaultdict(int)
t1 = time.time()
with pyhmmer.plan7.HMMFile(HMM) as hf:
    for hits in pyhmmer.hmmer.hmmsearch(hf, seqs, cpus=32):
        for h in hits:
            if h.evalue > 1e-5:
                continue
            cid = h.name.rsplit("_", 1)[0]
            hit[cid] += 1
print(f"打乱后搜索完成 {time.time()-t1:.1f}s", flush=True)

for r in csv.DictReader(open(V3), delimiter="\t"):
    tot[r["bucket"]] += 1

print("\n=== 真实 vs 打乱 (E<=1e-5) ===", flush=True)
print(f"{'bucket':<30}{'n':<7}{'真实命中%':<12}{'打乱命中%':<12}", flush=True)
real = defaultdict(int)
for r in csv.DictReader(open(V3), delimiter="\t"):
    if r["best_evalue"]:
        try:
            if float(r["best_evalue"]) <= 1e-5:
                real[r["bucket"]] += 1
        except Exception:
            pass
for b in ["KEEP", "no_viral_domain_but_blast", "short_fragment", "DROP"]:
    if not tot[b]:
        continue
    print(f"{b:<30}{tot[b]:<7}{100*real[b]/tot[b]:<12.1f}{100*hit[b]/tot[b]:<12.1f}", flush=True)
