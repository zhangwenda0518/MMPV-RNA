#!/usr/bin/env python3
"""
rescue_pipeline.py -- 病毒基因组三支路级联拯救 v3.0
==================================================

独立拯救脚本: 接收已聚类的 centroids + clusters + 拆分文件, 直接执行三支路拯救。
不重新聚类, 不重新过滤。聚类由 cluster_pipeline.py 完成, rescue 只做拯救。

流程:
  分支 A: CheckV 并行评估 → completeness>90% 输出; 失败 → 分支 B
  分支 B: Virseqimprover reads 延伸 (cluster 内多样本聚合) → CheckV; 失败 → 分支 C
  分支 C: BLASTN 参考比对

依赖: checkv, blastn, Virseqimprover.py
"""

import argparse, subprocess, sys, os, threading, shutil
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm


# ══════════════════════════════════════════════════════════════
# 工具函数
# ══════════════════════════════════════════════════════════════

def run(cmd, desc="", timeout=None, check=True):
    label = desc or " ".join(str(c) for c in cmd)
    print(f"  [RUN] {label}", flush=True)
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=check)
    except subprocess.CalledProcessError as e:
        print(f"  [FAIL] 退出码 {e.returncode}")
        if e.stderr:
            print(f"        {e.stderr.strip()[:500]}")
        if check:
            sys.exit(1)
        return None


def load_fasta_info(fasta_path):
    """加载 FASTA, 返回 {id: {"header": str, "seq": str, "length": int}}"""
    info = {}
    cur_id, cur_seq, cur_header = None, [], ""
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id:
                    s = "".join(cur_seq)
                    info[cur_id] = {"header": cur_header, "seq": s, "length": len(s)}
                cur_header, cur_id, cur_seq = line, line[1:].split()[0], []
            else:
                cur_seq.append(line)
        if cur_id:
            s = "".join(cur_seq)
            info[cur_id] = {"header": cur_header, "seq": s, "length": len(s)}
    return info


# ══════════════════════════════════════════════════════════════
# CheckV
# ══════════════════════════════════════════════════════════════

def run_checkv(fasta, out_dir, db, threads=4):
    d = Path(out_dir); d.mkdir(parents=True, exist_ok=True)
    run(["checkv", "completeness", str(fasta), str(d), "-d", str(db), "-t", str(threads)],
        f"CheckV: {Path(fasta).name}", check=False)
    return d / "completeness.tsv"


def parse_checkv(qs_path, threshold=90.0):
    """解析 completeness.tsv (aai_completeness 列); 返回 (pass_ids, fail_ids, skip_ids)
    skip_ids: aai_completeness=NA 的 contig (无法评估, 非病毒/无参考, 不进分支B)"""
    if not Path(qs_path).is_file():
        return set(), set(), set()
    df = pd.read_csv(qs_path, sep='\t')
    comp_col = 'aai_completeness' if 'aai_completeness' in df.columns else 'completeness'
    pass_ids = set()
    fail_ids = set()
    skip_ids = set()
    for _, row in df.iterrows():
        cid = str(row.get('contig_id', ''))
        val = row.get(comp_col, 0)
        if pd.isna(val) or str(val).strip() in ('NA', 'Not-determined', ''):
            skip_ids.add(cid)
        elif float(val) >= threshold:
            pass_ids.add(cid)
        else:
            fail_ids.add(cid)
    return pass_ids, fail_ids


# ══════════════════════════════════════════════════════════════
# BLASTN
# ══════════════════════════════════════════════════════════════

def run_blastn(query, db, out, threads=4):
    run(["blastn", "-task", "megablast", "-query", str(query), "-db", str(db),
         "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
         "-num_threads", str(threads), "-max_target_seqs", "10", "-out", str(out)],
        f"BLASTN: {Path(query).name}", check=False)


# ══════════════════════════════════════════════════════════════
# Virseqimprover 调用
# ══════════════════════════════════════════════════════════════

def run_vsi(ref_fa, r1, r2, out_dir, threads, vsi_path, salmon_bin, checkv_db=None):
    d = Path(out_dir); d.mkdir(parents=True, exist_ok=True)
    vsi = vsi_path if vsi_path and Path(vsi_path).is_file() else None
    if not vsi:
        print("  [WARN] Virseqimprover.py 未找到")
        return None, False
    cmd = ["python", str(vsi), "-1", str(r1), "-2", str(r2),
           "-scaffold", str(ref_fa), "-o", str(d),
           "-salmon", str(salmon_bin), "-t", str(threads)]
    if checkv_db:
        cmd += ["-checkv_db", str(checkv_db)]
    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=14400, check=False)
    except Exception:
        pass
    for name in ["scaffold.fasta", "pilon_out.fasta"]:
        f = d / name
        if f.is_file() and f.stat().st_size > 0:
            return str(f), True
    return None, False


# ══════════════════════════════════════════════════════════════
# Reads 匹配与聚合
# ══════════════════════════════════════════════════════════════

def _find_reads(fastq_dir, sample):
    """自动匹配 reads: _1/_2 或 .R1./.R2. 格式 + 多种扩展名 + 单端回退"""
    # PE: _1 / _2
    for suffix in ['', '_clean']:
        for ext in ['.fastq.gz', '.fq.gz', '.fa.gz', '.fasta.gz', '.fastq', '.fq', '.fa', '.fasta']:
            r1 = os.path.join(fastq_dir, f"{sample}{suffix}_1{ext}")
            r2 = os.path.join(fastq_dir, f"{sample}{suffix}_2{ext}")
            if os.path.isfile(r1) and os.path.isfile(r2):
                return r1, r2
    # PE: .R1. / .R2.
    for suffix in ['', '_clean']:
        for ext in ['.fastq.gz', '.fq.gz', '.fa.gz', '.fasta.gz']:
            r1 = os.path.join(fastq_dir, f"{sample}{suffix}.R1{ext}")
            r2 = os.path.join(fastq_dir, f"{sample}{suffix}.R2{ext}")
            if os.path.isfile(r1) and os.path.isfile(r2):
                return r1, r2
    # SE fallback
    for suffix in ['', '_clean']:
        for ext in ['.fastq.gz', '.fq.gz', '.fa.gz', '.fasta.gz', '.fastq', '.fq', '.fa', '.fasta']:
            r1 = os.path.join(fastq_dir, f"{sample}{suffix}{ext}")
            if os.path.isfile(r1):
                return r1, None  # SE: R2=None
    return None, None


def _gather_cluster_reads(fastq_dir, samples, work_dir, prefix="multi"):
    """收集 cluster 内所有样本的 reads, 合并为 R1/R2。返回 (r1_merged, r2_merged, sample_count)。
    支持 PE/SE 混合: SE 样本仅合并 R1, PE 样本合并 R1+R2。"""
    r1_files, r2_files = [], []
    paired_count, single_count, not_found = 0, 0, []
    for sample in sorted(samples):
        r1, r2 = _find_reads(fastq_dir, sample)
        if r1 and r2:
            r1_files.append(r1)
            r2_files.append(r2)
            paired_count += 1
        elif r1 and not r2:
            r1_files.append(r1)
            single_count += 1
        else:
            not_found.append(sample)

    if not_found:
        print(f"    [WARN] reads 未找到: {', '.join(not_found[:8])}{'...' if len(not_found)>8 else ''} "
              f"(在 {fastq_dir} 中, {len(not_found)}/{len(samples)} 样本)")

    if not r1_files:
        return None, None, 0

    n_total = paired_count + single_count
    if n_total == 1:
        # 单样本: 不合并, 直接返回
        return r1_files[0], r2_files[0] if r2_files else None, n_total

    work_dir = Path(work_dir); work_dir.mkdir(parents=True, exist_ok=True)
    r1_cat = str(work_dir / f"{prefix}_merged_R1.fastq.gz")
    r2_cat = str(work_dir / f"{prefix}_merged_R2.fastq.gz") if r2_files else None
    with open(r1_cat, "wb") as out:
        for f in r1_files:
            with open(f, "rb") as inf:
                shutil.copyfileobj(inf, out)
    if r2_files and r2_cat:
        with open(r2_cat, "wb") as out:
            for f in r2_files:
                with open(f, "rb") as inf:
                    shutil.copyfileobj(inf, out)
    se_info = f" + {single_count} SE" if single_count else ""
    print(f"    合并 {paired_count} PE{se_info} 样本 reads → {r1_cat}")
    return r1_cat, r2_cat, n_total


def _contig_cluster_samples(contig_id, clusters, fasta_info=None, max_samples=10):
    """返回 contig 所在 cluster 的样本 ID 集合 (去重, 按 contig 长度取 top-N)。
    若 contig 不在任何 cluster 中, 返回 None (参考序列, 无需 rescue)。"""
    for cname, info in clusters.items():
        if contig_id in info["members"]:
            if max_samples > 0 and fasta_info and len(info["members"]) > max_samples:
                ranked = sorted(info["members"],
                                key=lambda m: fasta_info.get(m, {}).get("length", 0),
                                reverse=True)
            else:
                ranked = list(info["members"])
            seen = set(); samples = []
            for m in ranked:
                s = m.split('_')[0]
                if s not in seen:
                    seen.add(s); samples.append(s)
                if max_samples > 0 and len(samples) >= max_samples: break
            return sorted(samples)
    # 不在任何 cluster → 可能是参考序列, 返回 None
    parts = contig_id.split('_')
    if len(parts) < 2 or parts[0].startswith('NC_') or parts[0].startswith('ref|'):
        return None
    return [parts[0]]


# ══════════════════════════════════════════════════════════════
# 分支 A: CheckV 并行评估
# ══════════════════════════════════════════════════════════════

def branch_a(ref_fasta, work_dir, checkv_db, threads, jobs):
    """CheckV 并行评估 centroids。返回 (pass_fa, fail_fa, pass_count, fail_count)"""
    d = Path(work_dir) / "branch_a"; d.mkdir(parents=True, exist_ok=True)
    records = list(SeqIO.parse(ref_fasta, "fasta"))
    total = len(records)
    if not records:
        print("  分支 A: 无 centroids, 跳过"); return None, None, 0, 0

    # 分块并行
    chunk_size = max(1, total // min(jobs, total))
    chunks = [records[i:i+chunk_size] for i in range(0, total, chunk_size)]
    print(f"  分支 A: {total} 条, {len(chunks)} 块 (jobs={jobs})")

    t_per = max(1, threads // min(jobs, len(chunks)))
    lock = threading.Lock()
    pass_records = []
    fail_records = []
    skip_records = []  # aai_completeness=NA, 不进B/C

    def _do(chunk):
        lp, lf, ls = [], [], []
        tmp_dir = d / f"chunk_{hash(str(chunk[0].id))}"
        tmp_dir.mkdir(exist_ok=True)
        chunk_fa = tmp_dir / "chunk.fasta"
        SeqIO.write(chunk, chunk_fa, "fasta")
        qs = run_checkv(chunk_fa, tmp_dir / "checkv_out", checkv_db, t_per)
        pids, fids, sids = parse_checkv(qs)
        for r in chunk:
            if r.id in pids:
                lp.append(r)
            elif r.id in sids:
                ls.append(r)
            else:
                lf.append(r)
        with lock:
            pass_records.extend(lp)
            fail_records.extend(lf)
            skip_records.extend(ls)

    if jobs > 1 and total > 1:
        with ThreadPoolExecutor(max_workers=min(jobs, len(chunks))) as ex:
            list(ex.map(_do, chunks))
    else:
        for chunk in chunks:
            _do(chunk)

    pass_fa = d / "branchA_pass.fasta"
    fail_fa = d / "branchA_fail.fasta"
    SeqIO.write(pass_records, pass_fa, "fasta")
    SeqIO.write(fail_records, fail_fa, "fasta")  # fail ≠ skip
    n_pass, n_fail, n_skip = len(pass_records), len(fail_records), len(skip_records)
    print(f"  分支 A: pass={n_pass}  fail={n_fail}  skip(NA)={n_skip}")
    return str(pass_fa) if n_pass > 0 else None, \
           str(fail_fa) if n_fail > 0 else None, \
           n_pass, n_fail


# ══════════════════════════════════════════════════════════════
# 分支 B: Virseqimprover reads 延伸
# ══════════════════════════════════════════════════════════════

def branch_b(fail_fa, fastq_dir, work_dir, checkv_db, threads, jobs, vsi_path, salmon_bin, clusters=None, fasta_info=None, max_vsi_samples=10, min_vsi_len=2000):
    """对 A 失败序列并行 Virseqimprover (cluster 内多样本 reads 聚合)。返回 (pass_fa, fail_fa, pass_count, fail_count)"""
    d = Path(work_dir) / "branch_b"; d.mkdir(parents=True, exist_ok=True)
    log_dir = d / "logs"; log_dir.mkdir(exist_ok=True)
    if not fail_fa or not os.path.isfile(fail_fa):
        print("  分支 B: 无失败序列, 跳过"); return None, None, 0, 0
    fail_records = list(SeqIO.parse(fail_fa, "fasta"))
    total = len(fail_records)
    if not fail_records:
        print("  分支 B: 无失败序列, 跳过"); return None, None, 0, 0

    # 长度过滤: < min_vsi_len 直升 branch C
    vsi_records = [r for r in fail_records if len(r.seq) >= min_vsi_len]
    short_records = [r for r in fail_records if len(r.seq) < min_vsi_len]
    if short_records:
        print(f"  分支 B: {len(short_records)} 条 <{min_vsi_len}bp → 跳过VSI直升分支C")
    if not vsi_records:
        print(f"  分支 B: 全部 <{min_vsi_len}bp, 全部直升分支C")
        return None, fail_fa, 0, total

    fail_records = vsi_records  # 仅 VSI 候选
    total = len(fail_records)

    scaffold_dir = d / "scaffolds"; scaffold_dir.mkdir(exist_ok=True)
    merged_reads_dir = d / "merged_reads"; merged_reads_dir.mkdir(exist_ok=True)
    threads_per_job = threads
    lock = threading.Lock()
    stats = {"ok": 0, "fail": 0, "skip": 0}

    def _do(rec):
        sid = rec.id
        sample = sid.split('_')[0]

        out_dir = d / f"out_{sid}"
        vsi_fa = out_dir / "scaffold.fasta"
        if not vsi_fa.is_file():
            vsi_fa = out_dir / "pilon_out.fasta"

        if vsi_fa.is_file() and vsi_fa.stat().st_size > 0:
            with lock: stats["ok"] += 1
            return sid, str(vsi_fa), "VSI 完成 (resume)"

        sample_ids = _contig_cluster_samples(sid, clusters, fasta_info, max_vsi_samples) if clusters else [sample]
        if sample_ids is None:
            with lock: stats["skip"] += 1
            return sid, None, "参考序列(不在cluster中)"
        r1, r2, nsamp = _gather_cluster_reads(fastq_dir, sample_ids, merged_reads_dir, prefix=sid[:80])
        if not r1:
            with lock: stats["skip"] += 1
            return sid, None, f"reads 未找到 (samples: {sample_ids[:5]}{'...' if len(sample_ids)>5 else ''})"

        sf = scaffold_dir / f"{sid}.fasta"
        SeqIO.write([rec], sf, "fasta")

        task_log = log_dir / f"{sid}.log"
        with open(task_log, "w") as lf:
            lf.write(f"# VSI multi-sample task: {sid}\n"
                     f"# Cluster samples ({nsamp}): {', '.join(sample_ids)}\n"
                     f"# Merged Reads: {r1}, {r2}\n# Scaffold: {sf}\n")

        _, ok = run_vsi(sf, r1, r2, out_dir, threads_per_job, vsi_path, salmon_bin, checkv_db)

        vsi_fa = out_dir / "scaffold.fasta"
        if not vsi_fa.is_file():
            vsi_fa = out_dir / "pilon_out.fasta"
        if ok and vsi_fa.is_file():
            with lock: stats["ok"] += 1
            print(f"    [C-OK] {sid[:60]}  (reads from {nsamp} samples)  ✓", flush=True)
            return sid, str(vsi_fa), f"VSI 完成 ({nsamp} samples)"
        with lock: stats["fail"] += 1
        return sid, None, "VSI 失败" if ok else "VSI 无输出"

    print(f"  并行启动 {jobs} 个 Virseqimprover (各 {threads_per_job} 线程)...")
    results = {}
    if jobs > 1 and total > 1:
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futures = {ex.submit(_do, rec): rec.id for rec in fail_records}
            for f in tqdm(as_completed(futures), total=total, desc="  分支 B", unit="task"):
                sid, fa, msg = f.result()
                if fa:
                    results[sid] = fa
    else:
        for rec in tqdm(fail_records, desc="  分支 B", unit="task"):
            sid, fa, msg = _do(rec)
            if fa:
                results[sid] = fa

    print(f"  分支 B 完成: ok={stats['ok']} fail={stats['fail']} skip={stats['skip']}")

    vsi_records = []
    for sid, fa in results.items():
        for rec in SeqIO.parse(fa, "fasta"):
            rec.id = sid; rec.description = ""
            vsi_records.append(rec)

    if not vsi_records:
        return None, fail_fa, 0, total

    extended_fa = d / "branchB_extended.fasta"
    SeqIO.write(vsi_records, extended_fa, "fasta")

    qs = run_checkv(extended_fa, d / "checkv_out", checkv_db, threads)
    pass_ids, fail_ids, _ = parse_checkv(qs)

    pass_fa = d / "branchB_pass.fasta"
    fail_fa_out = d / "branchB_fail.fasta"
    SeqIO.write([r for r in vsi_records if r.id in pass_ids], pass_fa, "fasta")
    fail_to_write = []
    for rec in fail_records:
        if rec.id in fail_ids:
            vsi_match = [r for r in vsi_records if r.id == rec.id]
            fail_to_write.append(vsi_match[0] if vsi_match else rec)
    # short_records (< min_vsi_len) 也加入 fail → 直升分支 C
    if short_records:
        fail_to_write.extend(short_records)
    SeqIO.write(fail_to_write, fail_fa_out, "fasta")

    n_fail_total = len(fail_ids) + len(short_records)
    print(f"  分支 B: pass={len(pass_ids):,}  fail={n_fail_total} (含{len(short_records)}短序列直升C)")
    return str(pass_fa) if Path(pass_fa).is_file() else None, \
           str(fail_fa_out) if Path(fail_fa_out).is_file() else None, \
           len(pass_ids), n_fail_total


# ══════════════════════════════════════════════════════════════
# 分支 C: BLASTN + CheckV
# ══════════════════════════════════════════════════════════════

def branch_c(fail_fa, fastq_dir, work_dir,
             checkv_db, blast_db, threads, jobs, vsi_path=None, salmon_bin=None, clusters=None, fasta_info=None, max_vsi_samples=10):
    """B 失败 → BLASTN + CheckV (纯比对评级, 不跑 VSI)"""
    d = Path(work_dir) / "branch_c"; d.mkdir(parents=True, exist_ok=True)

    if not fail_fa or not os.path.isfile(fail_fa):
        print("  分支 C: 无待处理序列"); return None, 0
    fail_records = list(SeqIO.parse(fail_fa, "fasta"))
    if not fail_records:
        print("  分支 C: 无待处理序列"); return None, 0

    if not blast_db:
        print("  分支 C: 无 BLAST 数据库, 跳过"); return None, 0

    lock = threading.Lock()
    complete = {}

    def _do(rec):
        ref = rec.id
        sub_dir = d / f"tmp_{ref[:30]}"; sub_dir.mkdir(exist_ok=True)
        tmp = sub_dir / f"{ref[:30]}.fa"
        SeqIO.write([rec], tmp, "fasta")

        run_blastn(tmp, blast_db, sub_dir / "blastn.tsv", threads)

        qs = run_checkv(tmp, sub_dir / "cv", checkv_db, threads)
        pids, _, _ = parse_checkv(qs)
        if ref in pids:
            with lock:
                complete[ref] = rec
                print(f"  [BLASTN-HQ] {ref[:50]}", flush=True)

    print(f"  分支 C: BLASTN {len(fail_records)} 条 (jobs={jobs})")
    if jobs > 1 and len(fail_records) > 1:
        with ThreadPoolExecutor(max_workers=min(jobs, len(fail_records))) as ex:
            list(tqdm(ex.map(_do, fail_records), total=len(fail_records), desc="  分支 C", unit="task"))
    else:
        for rec in tqdm(fail_records, desc="  分支 C", unit="task"):
            _do(rec)

    if not complete:
        return None, 0

    pass_fa = d / "branchC_pass.fasta"
    SeqIO.write(list(complete.values()), pass_fa, "fasta")
    print(f"  分支 C: pass={len(complete)}")
    return str(pass_fa), len(complete)



# ══════════════════════════════════════════════════════════════
# vclust (仅用于最终去重)
# ══════════════════════════════════════════════════════════════

def run_vclust_prefilter_align_cluster(fasta_in, out_dir, ani=0.95, qcov=0.85, threads=4, algorithm="leiden"):
    """vclust 三件套: prefilter → align → cluster。返回 clusters.tsv 路径"""
    d = Path(out_dir); d.mkdir(parents=True, exist_ok=True)
    pre, atsv, ctsv, ids = d / "vclust_prefilter.txt", d / "vclust_ani.tsv", d / "vclust_clusters.tsv", d / "vclust_ani.ids.tsv"
    run(["vclust", "prefilter", "-i", str(fasta_in), "-o", str(pre), "--min-ident", str(ani), "--threads", str(threads)], "vclust prefilter")
    run(["vclust", "align", "--filter", str(pre), "-i", str(fasta_in), "--out-ani", str(ani), "--out-qcov", str(qcov), "-o", str(atsv), "--out-aln", str(d / "vclust_ani.aln.tsv"), "--threads", str(threads)], "vclust align")
    run(["vclust", "cluster", "-i", str(atsv), "-o", str(ctsv), "--ids", str(ids), "--algorithm", algorithm, "--metric", "ani", "--ani", str(ani), "--qcov", str(qcov), "--out-repr"], "vclust cluster")
    return str(ctsv)


def parse_vclust_clusters(tsv):
    raw = defaultdict(list)
    with open(tsv) as f:
        f.readline()
        for line in f:
            p = line.strip().split()
            if len(p) >= 2:
                raw[p[1].strip()].append(p[0].strip())
    return {f"cluster_{cid}": {"ref": m[0], "members": list(set(m))} for cid, m in raw.items()}


# ══════════════════════════════════════════════════════════════
# 主流程
# ══════════════════════════════════════════════════════════════

def main():
    p = argparse.ArgumentParser(description="病毒基因组三支路级联拯救 v3.0")
    p.add_argument("--centroids", "-c", required=True, help="输入 centroids FASTA")
    p.add_argument("--clusters-tsv", required=True, help="vclust 聚类结果 TSV")
    p.add_argument("--split-dir", required=True, help="per-cluster 拆分文件目录 (split_fastas/)")
    p.add_argument("--output-dir", "-o", required=True)
    p.add_argument("--fastq-dir", "-fq", required=True, help="原始 reads 目录")
    p.add_argument("--checkv-db", "-cv", required=True)
    p.add_argument("--blast-db", "-db", default=None)
    p.add_argument("--virseqimprover-path", default=None)
    p.add_argument("--salmon-bin", default=os.path.expanduser("~/mambaforge/envs/Virseqimprover/bin/salmon"))
    p.add_argument("--max-vsi-samples", type=int, default=10, help="VSI 最大合并样本数 (0=不限制, 默认10)")
    p.add_argument("--min-vsi-len", type=int, default=2000, help="VSI 最小 contig 长度 bp (短于此值直升分支C, 默认2000)")
    p.add_argument("--threads", "-t", type=int, default=64)
    p.add_argument("--jobs", "-j", type=int, default=4, help="Virseqimprover 并行数")
    p.add_argument("--ani", type=float, default=0.95, help="最终 vclust ANI")
    p.add_argument("--qcov", type=float, default=0.85, help="最终 vclust QCOV")
    p.add_argument("--resume", action="store_true", help="断点续传")
    args = p.parse_args()

    start = datetime.now()
    print(f"开始: {start.strftime('%Y-%m-%d %H:%M:%S')}")

    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)

    # 1. 加载 centroids + clusters + fasta_info
    print("\n── 加载输入 ──")
    centroids_fa = Path(args.centroids)
    ctsv = Path(args.clusters_tsv)
    split_dir = Path(args.split_dir)

    clusters = parse_vclust_clusters(ctsv)
    print(f"  clusters: {len(clusters)} 个")

    # 加载 centroids FASTA
    centroids_records = list(SeqIO.parse(centroids_fa, "fasta"))
    print(f"  centroids: {len(centroids_records)} 条")

    # 加载所有 cluster 成员的序列 (从 split_fastas/)
    fasta_info = {}
    for fa in split_dir.glob("cluster_*.all.fasta"):
        fasta_info.update(load_fasta_info(str(fa)))
    print(f"  fasta_info: {len(fasta_info)} 条")

    # 统计长短序列分布 (短序列后续跳过B直升C)
    n_long = sum(1 for r in centroids_records if len(r.seq) >= args.min_vsi_len)
    n_short = len(centroids_records) - n_long
    if n_short > 0:
        print(f"  长度分布: ≥{args.min_vsi_len}bp={n_long} 条, <{args.min_vsi_len}bp={n_short} 条 (A后跳过B直升C)")

    # 2. 分支 A: CheckV 并行评估 (全部序列)
    print("\n── Step 1: 分支 A (CheckV) ──")
    fa_a_pass = out / "branch_a" / "branchA_pass.fasta"
    fa_a_fail = out / "branch_a" / "branchA_fail.fasta"
    if args.resume and fa_a_pass.is_file() and fa_a_pass.stat().st_size > 0:
        print("  [RESUME] 分支 A 已有结果, 跳过")
        cnt_a = sum(1 for _ in SeqIO.parse(fa_a_pass, "fasta"))
        cnt_a_fail = sum(1 for _ in SeqIO.parse(fa_a_fail, "fasta")) if fa_a_fail.is_file() else 0
        fa_a_pass, fa_a_fail = str(fa_a_pass), str(fa_a_fail)
    else:
        # 写 centroids 到临时文件供 branch_a
        tmp_centroids = out / "input_centroids.fasta"
        SeqIO.write(centroids_records, tmp_centroids, "fasta")
        fa_a_pass, fa_a_fail, cnt_a, cnt_a_fail = branch_a(str(tmp_centroids), out, args.checkv_db, args.threads, args.jobs)

    # 3. 分支 B: Virseqimprover
    print("\n── Step 2: 分支 B (VSI) ──")
    fa_b_pass = out / "branch_b" / "branchB_pass.fasta"
    fa_b_fail = out / "branch_b" / "branchB_fail.fasta"
    if args.resume and fa_b_pass.is_file() and fa_b_pass.stat().st_size > 0:
        print("  [RESUME] 分支 B 已有结果, 跳过")
        cnt_b = sum(1 for _ in SeqIO.parse(fa_b_pass, "fasta"))
        cnt_b_fail = sum(1 for _ in SeqIO.parse(fa_b_fail, "fasta")) if fa_b_fail.is_file() else 0
        fa_b_pass, fa_b_fail = str(fa_b_pass), str(fa_b_fail)
    else:
        fa_b_pass, fa_b_fail, cnt_b, cnt_b_fail = branch_b(fa_a_fail, args.fastq_dir, out, args.checkv_db, args.threads, args.jobs, args.virseqimprover_path, args.salmon_bin, clusters, fasta_info, args.max_vsi_samples, args.min_vsi_len)

    # 4. 分支 C: BLASTN (纯比对)
    print("\n── Step 3: 分支 C (BLASTN) ──")
    fa_c_pass = out / "branch_c" / "branchC_pass.fasta"
    if args.resume and fa_c_pass.is_file() and fa_c_pass.stat().st_size > 0:
        print("  [RESUME] 分支 C 已有结果, 跳过")
        cnt_c = sum(1 for _ in SeqIO.parse(fa_c_pass, "fasta"))
        fa_c_pass = str(fa_c_pass)
    elif args.blast_db:
        fail_for_c = fa_b_fail if (fa_b_fail and Path(fa_b_fail).is_file() and Path(fa_b_fail).stat().st_size > 0) else fa_a_fail
        fa_c_pass, cnt_c = branch_c(fail_for_c, args.fastq_dir, out, args.checkv_db, args.blast_db, args.threads, args.jobs)
    else:
        fa_c_pass, cnt_c = None, 0
        print("  [SKIP] 分支 C — 无 BLAST 数据库")

    # 5. 合并
    print("\n── Step 4: 合并 ──")
    d4 = out / "merged"; d4.mkdir(exist_ok=True)
    merged = d4 / "all_HQ.fasta"
    centroids_final = out / "centroids" / "final_centroids.fasta"

    if args.resume and centroids_final.is_file() and centroids_final.stat().st_size > 0:
        print("  [RESUME] Step 4+5 已有最终结果, 跳过")
        total_m = sum(1 for _ in SeqIO.parse(merged, "fasta")) if merged.is_file() else cnt_a + cnt_b + cnt_c
        n_final = sum(1 for _ in SeqIO.parse(centroids_final, "fasta"))
    else:
        with open(merged, "w") as mf:
            for fp in [fa_a_pass, fa_b_pass, fa_c_pass]:
                if fp and os.path.isfile(fp):
                    with open(fp) as inf:
                        mf.write(inf.read())
        total_m = sum(1 for _ in SeqIO.parse(merged, "fasta"))
        print(f"  A:{cnt_a}  B:{cnt_b}  C:{cnt_c}  →  {total_m} 条  →  {merged}")

        if total_m < 2:
            print("  [SKIP] <2 条, 跳过最终 vclust")
            centroids_final.parent.mkdir(parents=True, exist_ok=True)
            with open(centroids_final, "w") as cf, open(merged) as mf2:
                cf.write(mf2.read())
            n_final = total_m
        else:
            # 6. 最终 vclust 去重
            print("\n── Step 5: vclust 最终去重 ──")
            centroids_final.parent.mkdir(parents=True, exist_ok=True)
            fctsv = run_vclust_prefilter_align_cluster(str(merged), out / "centroids", args.ani, args.qcov, args.threads)
            fclusters = parse_vclust_clusters(fctsv)
            all_seqs = SeqIO.to_dict(SeqIO.parse(merged, "fasta"))
            cen_records = [all_seqs[ci["ref"]] for ci in fclusters.values() if ci["ref"] in all_seqs]
            SeqIO.write(cen_records, centroids_final, "fasta")
            n_final = len(cen_records)

    # 清理 merged_reads 释放磁盘
    for bd in out.glob("branch_b/merged_reads"):
        if bd.is_dir():
            shutil.rmtree(bd, ignore_errors=True)
            print(f"\n  清理: {bd}")

    elapsed = (datetime.now() - start).total_seconds()
    print(f"\n{'=' * 60}")
    print(f"  分支 A (CheckV):      {cnt_a:,}")
    print(f"  分支 B (VSI):         {cnt_b:,}")
    print(f"  分支 C (BLASTN):  {cnt_c:,}")
    print(f"  最终无冗余:           {n_final:,}")
    print(f"  最终输出:             {centroids_final}")
    print(f"  耗时:                 {elapsed / 60:.1f} min")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
