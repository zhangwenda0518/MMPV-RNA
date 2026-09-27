#!/usr/bin/env python3
"""RVDB HMM pilot v3: 覆盖率改用 aln.hmm_length / aln.hmm_name (不用 hits.hmm)。"""
import csv
import time
from collections import defaultdict

import pyhmmer
from pyhmmer.easel import SequenceFile

OUT = "/home/zhangwenda/rvdb_hmm_pilot"
S = "/home/zhangwenda/data-test/out/09b_Analysis_Verify/virus_validation/rescue_evidence_scored.tsv"
HMM = "/home/zhangwenda/database/virus-db/RVDB-v31/U-RVDBv31.0-prot.hmm"
clean_path = f"{OUT}/all_proteins.clean.faa"

meta = {}
with open(S) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        fr = r.get("fp_reason") or ""
        if r["verdict"] == "REVIEW":
            bucket = ("no_viral_domain_but_blast" if fr.startswith("no_viral_domain_but_blast")
                      else "short_fragment" if fr.startswith("short_fragment") else "other_REVIEW")
        else:
            bucket = r["verdict"]
        meta[r["contig_id"]] = (r["verdict"], r.get("cdd_evidence", ""), bucket)

with SequenceFile(clean_path, digital=True) as sf:
    seqs = sf.read_block()

per = defaultdict(lambda: {"n": 0, "n_cov50": 0, "best_ev": 1.0, "best_sc": 0.0,
                           "best_cov": 0.0, "best_hmm": "", "best_orf": 0})
t1 = time.time()
with pyhmmer.plan7.HMMFile(HMM) as hf:
    for hits in pyhmmer.hmmer.hmmsearch(hf, seqs, cpus=32):
        for hit in hits:
            ev = hit.evalue
            if ev > 1e-2:
                continue
            cid = hit.name.rsplit("_", 1)[0]
            cov, ol, hn = 0.0, 0, ""
            try:
                doms = list(hit.domains)
                if doms:
                    bd = max(doms, key=lambda x: getattr(x, "score", 0))
                    aln = bd.alignment
                    hl = aln.hmm_length
                    cov = (aln.hmm_to - aln.hmm_from + 1) / hl if hl else 0.0
                    ol = aln.target_length or 0
                    hn = aln.hmm_name
                    if isinstance(hn, bytes):
                        hn = hn.decode()
            except Exception:
                pass
            d = per[cid]
            d["n"] += 1
            if cov >= 0.5:
                d["n_cov50"] += 1
            if ev < d["best_ev"]:
                d["best_ev"] = ev
                d["best_hmm"], d["best_cov"], d["best_orf"] = hn, cov, ol
            if hit.score > d["best_sc"]:
                d["best_sc"] = hit.score
print(f"搜索完成 {time.time()-t1:.1f}s, 涉及 contig {len(per)}", flush=True)

tsv = f"{OUT}/rvdb_hmm_per_contig_v3.tsv"
with open(tsv, "w", newline="") as fo:
    w = csv.writer(fo, delimiter="\t")
    w.writerow(["contig_id", "verdict", "cdd_evidence", "bucket", "n_hits", "n_hits_cov50",
                "best_evalue", "best_score", "best_cov", "best_orf_len", "best_hmm"])
    for cid in sorted(meta):
        v, ce, b = meta[cid]
        d = per.get(cid, {"n": 0, "n_cov50": 0, "best_ev": 1.0, "best_sc": 0.0,
                          "best_cov": 0.0, "best_hmm": "", "best_orf": 0})
        w.writerow([cid, v, ce, b, d["n"], d["n_cov50"],
                    f"{d['best_ev']:.2e}" if d["n"] else "",
                    f"{d['best_sc']:.1f}" if d["n"] else "",
                    f"{d['best_cov']:.3f}" if d["n"] else "",
                    d["best_orf"], d["best_hmm"]])

print("\n=== 分离度 (覆盖率已修) ===", flush=True)
g = defaultdict(lambda: defaultdict(int))
for cid, (v, ce, b) in meta.items():
    d = per.get(cid)
    gg = g[b]
    gg["n"] += 1
    if not d:
        continue
    if d["best_ev"] <= 1e-5:
        gg["e5"] += 1
    if d["best_ev"] <= 1e-5 and d["best_cov"] >= 0.5:
        gg["e5c50"] += 1
    if d["best_ev"] <= 1e-10 and d["best_cov"] >= 0.5:
        gg["e10c50"] += 1
    if d["best_ev"] <= 1e-5 and d["best_cov"] >= 0.8:
        gg["e5c80"] += 1
    if d["n_cov50"] >= 1 and d["best_ev"] <= 1e-5:
        gg["anycov50"] += 1
print(f"{'bucket':<30}{'n':<7}{'E5':<8}{'E5&cov50':<10}{'E5&cov80':<10}{'E10&cov50':<11}{'有cov50命中':<12}", flush=True)
for b in ["KEEP", "no_viral_domain_but_blast", "short_fragment", "DROP"]:
    gg = g.get(b)
    if not gg:
        continue
    n = gg["n"]
    print(f"{b:<30}{n:<7}{100*gg['e5']/n:<8.1f}{100*gg['e5c50']/n:<10.1f}"
          f"{100*gg['e5c80']/n:<10.1f}{100*gg['e10c50']/n:<11.1f}{100*gg['anycov50']/n:<12.1f}", flush=True)
print(f"\n输出 {tsv}", flush=True)
