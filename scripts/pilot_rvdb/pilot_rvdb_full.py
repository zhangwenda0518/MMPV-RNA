#!/usr/bin/env python3
"""RVDB HMM pilot 正式跑: 全 8915 条 (路径②为目标, KEEP 正对照, DROP 负对照)。
产出: per-contig 命中 TSV + 分组列联表。复用已有的 ORF 中间文件。
"""
import csv
import os
import subprocess
import time
from collections import defaultdict

import pyhmmer
from pyhmmer.easel import SequenceFile

OUT = "/home/zhangwenda/rvdb_hmm_pilot"
S = "/home/zhangwenda/data-test/out/09b_Analysis_Verify/virus_validation/rescue_evidence_scored.tsv"
HQ = "/home/zhangwenda/data-test/out/09_Virome_Analysis/HQ_plant_viruses.fasta"
HMM = "/home/zhangwenda/database/virus-db/RVDB-v31/U-RVDBv31.0-prot.hmm"
clean_path = f"{OUT}/all_proteins.clean.faa"
allfa = f"{OUT}/all_8915.fasta"

# ── 1. 元信息 ──
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
print(f"总条数 {len(meta)}", flush=True)

# ── 2. 中间文件复用 ──
if not (os.path.exists(allfa) and os.path.getsize(allfa) > 1000):
    ids = set(meta)
    n = 0
    with open(allfa, "w") as fo, open(HQ) as fi:
        w = False
        for l in fi:
            if l.startswith(">"):
                w = l[1:].split()[0] in ids
                if w:
                    n += 1
            if w:
                fo.write(l)
    print(f"抽取 {n} 条", flush=True)

if not (os.path.exists(clean_path) and os.path.getsize(clean_path) > 1000):
    t0 = time.time()
    r = subprocess.run(["prodigal", "-p", "meta", "-c", "-i", allfa,
                        "-a", f"{OUT}/all_proteins.faa", "-o", "/dev/null"],
                       capture_output=True, text=True)
    print(f"prodigal rc={r.returncode} ({time.time()-t0:.1f}s)", flush=True)
    with open(clean_path, "w") as fo, open(f"{OUT}/all_proteins.faa") as fi:
        for l in fi:
            fo.write(l.rstrip("\n") if l.startswith(">") else l.replace("*", "").replace(" ", ""))
n_orf = sum(1 for l in open(clean_path) if l.startswith(">"))
print(f"ORF 数 {n_orf}", flush=True)

# ── 3. 全库搜索 ──
with SequenceFile(clean_path, digital=True) as sf:
    seqs = sf.read_block()

per = defaultdict(lambda: {"n": 0, "best_ev": 1.0, "best_sc": 0.0, "best_cov": 0.0, "best_hmm": ""})
t1 = time.time()
ntot = 0
with pyhmmer.plan7.HMMFile(HMM) as hf:
    for hits in pyhmmer.hmmer.hmmsearch(hf, seqs, cpus=32):
        qname, qlen = "", 0
        try:
            h = hits.hmm
            if h is not None:
                nm = h.name
                qname = nm.decode() if isinstance(nm, bytes) else str(nm)
                qlen = h.M
        except Exception:
            pass
        if not qname:
            try:
                q = hits.query
                qname = q.decode() if isinstance(q, bytes) else str(q)
            except Exception:
                qname = "?"
        for hit in hits:
            ev = hit.evalue
            if ev > 1e-2:
                continue
            cid = hit.name.rsplit("_", 1)[0]
            ntot += 1
            d = per[cid]
            d["n"] += 1
            if ev < d["best_ev"]:
                d["best_ev"] = ev
                d["best_hmm"] = qname
                try:
                    aln = hit.alignment
                    d["best_cov"] = (aln.hmm_to - aln.hmm_from + 1) / qlen if qlen else 0.0
                except Exception:
                    d["best_cov"] = 0.0
            if hit.score > d["best_sc"]:
                d["best_sc"] = hit.score
print(f"搜索完成 {time.time()-t1:.1f}s, 命中记录 {ntot}, 涉及 contig {len(per)}", flush=True)

# ── 4. per-contig 输出 ──
tsv = f"{OUT}/rvdb_hmm_per_contig.tsv"
with open(tsv, "w", newline="") as fo:
    w = csv.writer(fo, delimiter="\t")
    w.writerow(["contig_id", "verdict", "cdd_evidence", "bucket", "n_hits",
                "best_evalue", "best_score", "best_cov", "best_hmm"])
    for cid in sorted(meta):
        v, ce, b = meta[cid]
        d = per.get(cid, {"n": 0, "best_ev": 1.0, "best_sc": 0.0, "best_cov": 0.0, "best_hmm": ""})
        w.writerow([cid, v, ce, b, d["n"],
                    f"{d['best_ev']:.2e}" if d["n"] else "",
                    f"{d['best_sc']:.1f}" if d["n"] else "",
                    f"{d['best_cov']:.2f}" if d["n"] else "",
                    d["best_hmm"]])

# ── 5. 列联表 ──
print("\n=== 分组命中率 (E<=1e-2 记为命中) ===", flush=True)
print(f"{'bucket':<30}{'n':<7}{'命中n':<8}{'命中率%':<10}{'E<=1e-5%':<10}{'E<=1e-10%':<10}", flush=True)
groups = defaultdict(lambda: [0, 0, 0, 0])
for cid, (v, ce, b) in meta.items():
    d = per.get(cid)
    g = groups[b]
    g[0] += 1
    if d:
        g[1] += 1
        if d["best_ev"] <= 1e-5:
            g[2] += 1
        if d["best_ev"] <= 1e-10:
            g[3] += 1
for b in ["KEEP", "no_viral_domain_but_blast", "short_fragment", "DROP", "other_REVIEW"]:
    if b not in groups:
        continue
    g = groups[b]
    if g[0] == 0:
        continue
    print(f"{b:<30}{g[0]:<7}{g[1]:<8}{100*g[1]/g[0]:<10.1f}{100*g[2]/g[0]:<10.1f}{100*g[3]/g[0]:<10.1f}", flush=True)
print(f"\n输出: {tsv}", flush=True)
