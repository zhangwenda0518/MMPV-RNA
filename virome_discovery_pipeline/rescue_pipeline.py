#!/usr/bin/env python3
"""
rescue_pipeline.py -- 病毒基因组三支路级联拯救 v3.0
==================================================

独立拯救脚本: 接收已聚类的 centroids + clusters + 拆分文件, 直接执行三支路拯救。
不重新聚类, 不重新过滤。聚类由 cluster_pipeline.py 完成, rescue 只做拯救。

流程:
  分支 A: CheckV 并行评估 → completeness>90% 输出; 失败 → 分支 B
  分支 B: Virseqimprover reads 延伸 (cluster 内多样本聚合) → CheckV; 失败 → 分支 C
  分支 C: BLASTN 参考比对 + ragtag 参考引导延伸 → CheckV; 失败 → 分支 D
  分支 D: genus_len 属水平长度拯救 (同属物种长度 ±15%)

依赖: checkv, blastn, Virseqimprover.py
"""

import argparse, subprocess, sys, os, threading, shutil, json
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
    return pass_ids, fail_ids, skip_ids


# ══════════════════════════════════════════════════════════════
# BLASTN
# ══════════════════════════════════════════════════════════════

def run_blastn(query, db, out, threads=4):
    """BLASTN 搜索: dc-megablast → outfmt 6 std qlen slen (供 WVDB 坐标合并)"""
    run(["blastn", "-task", "dc-megablast", "-query", str(query), "-db", str(db),
         "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen",
         "-num_threads", str(threads), "-max_target_seqs", "5",
         "-evalue", "1e-5", "-word_size", "11",
         "-out", str(out)],
        f"BLASTN: {Path(query).name}", check=False)


# ══════════════════════════════════════════════════════════════
# Virseqimprover 调用
# ══════════════════════════════════════════════════════════════

def run_vsi(ref_fa, r1, r2, out_dir, threads, vsi_path, salmon_bin, checkv_db=None, genus_avg_len=0):
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
    if genus_avg_len > 0:
        cmd += ["-genus_avg_len", str(genus_avg_len)]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=14400, check=False)
    # 保存 VSI 输出
    vsi_log = d / "vsi_stdout.log"
    with open(vsi_log, 'w') as lf:
        lf.write(result.stdout or "")
        if result.stderr:
            lf.write("\n=== STDERR ===\n")
            lf.write(result.stderr)
    # 仅当 VSI 正常退出 (returncode=0) 且输出文件存在才认可
    if result.returncode == 0:
        for name in ["scaffold.fasta", "pilon_out.fasta"]:
            f = d / name
            if f.is_file() and f.stat().st_size > 0:
                return str(f), True
    return None, False


# ══════════════════════════════════════════════════════════════
# Reads 匹配与聚合
# ══════════════════════════════════════════════════════════════

def _find_reads(fastq_dir, sample):
    """自动匹配 reads: 标准命名 + glob 模糊匹配 (兼容 S2_L001_R1 等 Illumina 格式)"""
    import glob as _glob
    # 1. 标准 PE: _1 / _2
    for suffix in ['', '_clean']:
        for ext in ['.fastq.gz', '.fq.gz', '.fa.gz', '.fasta.gz', '.fastq', '.fq', '.fa', '.fasta']:
            r1 = os.path.join(fastq_dir, f"{sample}{suffix}_1{ext}")
            r2 = os.path.join(fastq_dir, f"{sample}{suffix}_2{ext}")
            if os.path.isfile(r1) and os.path.isfile(r2):
                return r1, r2
    # 2. PE: .R1. / .R2.
    for suffix in ['', '_clean']:
        for ext in ['.fastq.gz', '.fq.gz', '.fa.gz', '.fasta.gz']:
            r1 = os.path.join(fastq_dir, f"{sample}{suffix}.R1{ext}")
            r2 = os.path.join(fastq_dir, f"{sample}{suffix}.R2{ext}")
            if os.path.isfile(r1) and os.path.isfile(r2):
                return r1, r2
    # 3. Glob 模糊匹配: Illumina 格式 {sample}_S*_L*_R1*.fastq.gz 等
    for ext in ['.fastq.gz', '.fq.gz', '.fa.gz', '.fasta.gz', '.fastq.gz', '.gz']:
        r1_files = sorted(_glob.glob(os.path.join(fastq_dir, f"{sample}*_R1*{ext}")))
        r2_files = sorted(_glob.glob(os.path.join(fastq_dir, f"{sample}*_R2*{ext}")))
        if r1_files and r2_files:
            return r1_files[0], r2_files[0]
        r1_files = sorted(_glob.glob(os.path.join(fastq_dir, f"{sample}*_1*{ext}")))
        r2_files = sorted(_glob.glob(os.path.join(fastq_dir, f"{sample}*_2*{ext}")))
        if r1_files and r2_files:
            return r1_files[0], r2_files[0]
    # 4. SE fallback
    for suffix in ['', '_clean']:
        for ext in ['.fastq.gz', '.fq.gz', '.fa.gz', '.fasta.gz', '.fastq', '.fq', '.fa', '.fasta']:
            r1 = os.path.join(fastq_dir, f"{sample}{suffix}{ext}")
            if os.path.isfile(r1):
                return r1, None
    # 5. Co-assembly fallback: 找不到样本 reads → 用 ALL_merged 文件
    merged_r1 = os.path.join(fastq_dir, "ALL_merged_R1.fq.gz")
    merged_r2 = os.path.join(fastq_dir, "ALL_merged_R2.fq.gz")
    if os.path.isfile(merged_r1) and os.path.isfile(merged_r2):
        return merged_r1, merged_r2
    if os.path.isfile(merged_r1):
        return merged_r1, None
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

    # 检测 reads 格式 (FASTA vs FASTQ)
    is_fasta = r1_files[0].endswith(('.fa.gz', '.fasta.gz', '.fa', '.fasta'))
    ext = '.fa.gz' if is_fasta else '.fastq.gz'

    work_dir = Path(work_dir); work_dir.mkdir(parents=True, exist_ok=True)
    r1_cat = str(work_dir / f"{prefix}_merged_R1{ext}")
    r2_cat = str(work_dir / f"{prefix}_merged_R2{ext}") if r2_files else None
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


def _contig_cluster_samples(contig_id, clusters, fasta_info=None, max_samples=10, flye_sample_map=None):
    """返回 contig 所在 cluster 的样本 ID 集合 (去重, 按 contig 长度取 top-N)。
    Flye contig (contig_XXX) 优先使用 flye_sample_map 精确映射;
    原始 contig 从 vclust clusters 解析样本前缀。
    NC_/ref| 前缀返回 None (参考序列, 无需 rescue)。"""
    # Flye contig → 精确映射优先
    if flye_sample_map and contig_id in flye_sample_map:
        source_samples = flye_sample_map[contig_id]
        if max_samples > 0 and len(source_samples) > max_samples:
            source_samples = source_samples[:max_samples]
        return source_samples
    # 原始 contig → 查 clusters
    if clusters:
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
    # 从 contig_id 解析样本前缀
    parts = contig_id.split('_')
    if len(parts) < 2 or parts[0].startswith('NC_') or parts[0].startswith('ref|'):
        return None
    return [parts[0]]


# ══════════════════════════════════════════════════════════════
# 分支 A: CheckV 并行评估
# ══════════════════════════════════════════════════════════════

def branch_a(ref_fasta, work_dir, checkv_db, threads, jobs, threshold=90.0):
    """CheckV 并行评估 centroids。返回 (pass_fa, fail_fa, pass_count, fail_count)。支持断点续传。"""
    d = Path(work_dir) / "branch_a"; d.mkdir(parents=True, exist_ok=True)

    # 断点续传: 最终输出文件已存在 → 直接恢复
    pass_fa = d / "branchA_pass.fasta"
    fail_fa = d / "branchA_fail.fasta"
    if pass_fa.is_file() and fail_fa.is_file():
        n_pass = sum(1 for _ in SeqIO.parse(str(pass_fa), "fasta"))
        n_fail = sum(1 for _ in SeqIO.parse(str(fail_fa), "fasta"))
        if n_pass + n_fail > 0:
            print(f"  分支 A: [resume] pass={n_pass} fail={n_fail} → 跳过 CheckV")
            return str(pass_fa) if n_pass > 0 else None, \
                   str(fail_fa) if n_fail > 0 else None, \
                   n_pass, n_fail

    records = list(SeqIO.parse(ref_fasta, "fasta"))
    total = len(records)
    if not records:
        print("  分支 A: 无 centroids, 跳过"); return None, None, 0, 0

    # 分块并行 (用 chunk 索引确保目录名稳定 → 断点续传)
    chunk_size = max(1, total // min(jobs, total))
    chunks = [(i, records[i:i+chunk_size]) for i in range(0, total, chunk_size)]
    print(f"  分支 A: {total} 条, {len(chunks)} 块 (jobs={jobs})")

    t_per = max(1, threads // min(jobs, len(chunks)))
    lock = threading.Lock()
    pass_records = []
    fail_records = []
    skip_records = []  # aai_completeness=NA, 不进B/C

    def _do(args):
        chunk_idx, chunk = args
        lp, lf, ls = [], [], []
        tmp_dir = d / f"chunk_{chunk_idx}"
        tmp_dir.mkdir(exist_ok=True)
        chunk_fa = tmp_dir / "chunk.fasta"
        cv_tsv = tmp_dir / "checkv_out" / "completeness.tsv"

        # 断点续传: completeness.tsv 已存在 → 直接解析
        if cv_tsv.is_file():
            pids, fids, sids = parse_checkv(cv_tsv, threshold)
        else:
            SeqIO.write(chunk, chunk_fa, "fasta")
            qs = run_checkv(chunk_fa, tmp_dir / "checkv_out", checkv_db, t_per)
            pids, fids, sids = parse_checkv(qs, threshold)
        for r in chunk:
            if r.id in pids:
                lp.append(r)
            elif r.id in sids:
                # NA 但 ≥2000bp → 还有机会, 当 fail 处理进分支 B/C
                # NA 且 <2000bp → 短片段无蛋白基因, 跳过
                if len(r.seq) >= 2000:
                    lf.append(r)
                else:
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

    SeqIO.write(pass_records, pass_fa, "fasta")
    SeqIO.write(fail_records, fail_fa, "fasta")  # fail = completeness<阈值 + NA且≥2000bp
    n_pass, n_fail, n_skip = len(pass_records), len(fail_records), len(skip_records)
    print(f"  分支 A: pass={n_pass}  fail={n_fail}  skip(NA+short)={n_skip}")
    return str(pass_fa) if n_pass > 0 else None, \
           str(fail_fa) if n_fail > 0 else None, \
           n_pass, n_fail


# ══════════════════════════════════════════════════════════════
# 分支 B: Virseqimprover reads 延伸
# ══════════════════════════════════════════════════════════════

def branch_b(fail_fa, fastq_dir, work_dir, checkv_db, threads, jobs, vsi_path, salmon_bin, clusters=None, fasta_info=None, max_vsi_samples=10, min_vsi_len=2000, threshold=90.0, genus_map=None, genus_lens=None, flye_sample_map=None):
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
        # 优先取 scaffold-truncated (完整延伸版), 回退 scaffold.fasta → pilon_out.fasta
        vsi_fa = out_dir / "scaffold-truncated" / "scaffold.fasta"
        if not vsi_fa.is_file():
            vsi_fa = out_dir / "scaffold.fasta"
        if not vsi_fa.is_file():
            vsi_fa = out_dir / "pilon_out.fasta"

        if vsi_fa.is_file() and vsi_fa.stat().st_size > 0:
            # 验证 VSI 确实正常完成 (stdout 有 Finished 且无 Traceback)
            vsi_log = out_dir / "vsi_stdout.log"
            vsi_completed = False
            if vsi_log.is_file():
                with open(vsi_log) as lf:
                    content = lf.read()
                    if "Finished growing scaffold" in content and "Traceback" not in content:
                        vsi_completed = True
            if vsi_completed:
                with lock: stats["ok"] += 1
                return sid, str(vsi_fa), "VSI 完成 (resume)"
            else:
                # VSI 之前崩溃了, 删除残留重新跑
                import shutil as _shutil
                _shutil.rmtree(str(out_dir), ignore_errors=True)
                print(f"    [RE-RUN] {sid[:50]} | VSI 之前崩溃, 重新运行", flush=True)

        sample_ids = _contig_cluster_samples(sid, clusters, fasta_info, max_vsi_samples, flye_sample_map=flye_sample_map) if clusters else [sample]
        if sample_ids is None:
            with lock: stats["skip"] += 1
            return sid, None, "参考序列(不在cluster中)"
        r1, r2, nsamp = _gather_cluster_reads(fastq_dir, sample_ids, merged_reads_dir, prefix=sid[:80])
        if not r1:
            with lock: stats["skip"] += 1
            return sid, None, f"reads 未找到 (samples: {sample_ids[:5]}{'...' if len(sample_ids)>5 else ''})"

        sf = scaffold_dir / f"{sid}.fasta"
        SeqIO.write([rec], sf, "fasta")

        # 查找属平均长度 (CheckV NA 时的备选截止条件)
        gal = 0
        has_genus = False
        if genus_map and genus_lens:
            g = genus_map.get(sid) or genus_map.get(sample, "")
            g_clean = g.replace("g__", "").replace("G__", "")
            gal = genus_lens.get(g_clean, genus_lens.get(g, 0))
            has_genus = (gal > 0)

        # CheckV NA 且无 genus → VSI 无法获知截止条件, 跳过直接进分支 C
        if not has_genus:
            with lock:
                stats["skip"] += 1
            print(f"    [SKIP-VSI] {sid[:60]} | CheckV NA + no genus => branch C", flush=True)
            return sid, None, "无属分类(B分支C)"

        task_log = log_dir / f"{sid}.log"
        with open(task_log, "w") as lf:
            lf.write(f"# VSI multi-sample task: {sid}\n"
                     f"# Cluster samples ({nsamp}): {', '.join(sample_ids)}\n"
                     f"# Merged Reads: {r1}, {r2}\n# Scaffold: {sf}\n"
                     f"# Genus: {g if genus_map else 'N/A'} → avg_len: {gal}bp\n")

        _, ok = run_vsi(sf, r1, r2, out_dir, threads_per_job, vsi_path, salmon_bin, checkv_db, genus_avg_len=gal)

        vsi_fa = out_dir / "scaffold-truncated" / "scaffold.fasta"
        if not vsi_fa.is_file():
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
    pass_ids, fail_ids, skip_ids = parse_checkv(qs, threshold)

    # ── CheckV=NA 或 <threshold 的 genus 回退拯救 ──
    # VSI 自身已停止延伸, CheckV 数值不可靠时, 用属平均长度兜底
    genus_rescued = set()
    genus_tol = 0.85  # 属长度容忍度 (同分支D)
    rescue_candidates = set()
    if skip_ids: rescue_candidates.update(skip_ids)
    if fail_ids: rescue_candidates.update(fail_ids)
    if rescue_candidates and genus_map and genus_lens:
        for rec in vsi_records:
            if rec.id in rescue_candidates:
                gal = 0
                for key in (rec.id, rec.id.split('_')[0]):
                    g = genus_map.get(key, '')
                    g_clean = g.replace('g__', '').replace('G__', '')
                    gal = genus_lens.get(g_clean, genus_lens.get(g, 0))
                    if gal > 0:
                        break
                if gal > 0 and len(rec.seq) >= gal * genus_tol:
                    pass_ids.add(rec.id)
                    genus_rescued.add(rec.id)
        if genus_rescued:
            source = 'NA' if (genus_rescued & skip_ids) else ''
            if genus_rescued & fail_ids:
                source += ('+' if source else '') + 'fail'
            print(f"  [genus rescue] CheckV={source} 但长度达标: {len(genus_rescued)} 条 (≥{genus_tol*100:.0f}% 属平均长度)")

    pass_fa = d / "branchB_pass.fasta"
    fail_fa_out = d / "branchB_fail.fasta"
    SeqIO.write([r for r in vsi_records if r.id in pass_ids], pass_fa, "fasta")
    fail_to_write = []
    for rec in fail_records:
        in_fail = rec.id in fail_ids and rec.id not in genus_rescued
        in_skip = rec.id in skip_ids and rec.id not in genus_rescued
        if in_fail or in_skip:
            vsi_match = [r for r in vsi_records if r.id == rec.id]
            fail_to_write.append(vsi_match[0] if vsi_match else rec)
    # short_records (< min_vsi_len) 也加入 fail → 直升分支 C
    if short_records:
        fail_to_write.extend(short_records)
    SeqIO.write(fail_to_write, fail_fa_out, "fasta")

    n_fail_total = len(fail_ids - genus_rescued) + len(skip_ids - genus_rescued) + len(short_records)
    na_fail = len(skip_ids - genus_rescued)
    num_fail = len(fail_ids - genus_rescued)
    print(f"  分支 B: pass={len(pass_ids):,}  fail={n_fail_total} (含{len(short_records)}短序列直升C, {na_fail}条NA直升C, {num_fail}条<阈值直升C)")
    return str(pass_fa) if Path(pass_fa).is_file() else None, \
           str(fail_fa_out) if Path(fail_fa_out).is_file() else None, \
           len(pass_ids), n_fail_total


# ══════════════════════════════════════════════════════════════
# 分支 C: BLASTN + ragtag 多contig参考引导拼接 (N填充)
# ══════════════════════════════════════════════════════════════

def _load_clusters(clusters_tsv):
    """加载 vclust 聚类结果 → {contig_id: cluster_rep_id}"""
    cmap = {}
    if not clusters_tsv or not Path(clusters_tsv).is_file():
        return cmap
    with open(clusters_tsv) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                cmap[parts[0].strip()] = parts[1].strip()
    return cmap


def branch_c(fail_fa, fastq_dir, work_dir,
             checkv_db, blast_db, threads, jobs, vsi_path=None, salmon_bin=None,
             clusters=None, clustermap_tsv=None, centroids_fa=None, split_dir=None, fasta_info=None,
             max_vsi_samples=10, threshold=90.0, taxonomy_tsv=None, cdhit_fa=None, flye_sample_map=None):
    """B 失败 → 按属分组多contig + BLASTN自建库参考 + ragtag拼接 + CheckV质量标注"""
    d = Path(work_dir) / "branch_c"; d.mkdir(parents=True, exist_ok=True)

    if not fail_fa or not os.path.isfile(fail_fa):
        print("  分支 C: 无待处理序列"); return None, 0
    fail_records = list(SeqIO.parse(fail_fa, "fasta"))
    if not fail_records:
        print("  分支 C: 无待处理序列"); return None, 0

    # 加载 cluster 信息 + centroids
    cmap = {}
    if clustermap_tsv:
        cmap = _load_clusters(clustermap_tsv)
    centroids_all = {}
    cluster_split = {}  # {cluster_id: [records]} from split_fastas
    if centroids_fa and Path(centroids_fa).is_file():
        for rec in SeqIO.parse(centroids_fa, "fasta"):
            centroids_all[rec.id] = rec
    if split_dir and Path(split_dir).is_dir():
        for sp in Path(split_dir).glob("cluster_*.all.fasta"):
            cid = sp.stem.replace("cluster_", "").replace(".all", "")
            cluster_split[cid] = list(SeqIO.parse(sp, "fasta"))
    if cmap and centroids_all:
        print(f"  分支 C: clusters={len(set(cmap.values()))}, centroids={len(centroids_all)}")

    # 查找工具 (建库前)
    ragtag_bin = shutil.which("ragtag.py") or os.path.expanduser("~/mambaforge/bin/ragtag.py")
    blastdbcmd_bin = shutil.which("blastdbcmd") or "blastdbcmd"
    bwa_bin = shutil.which("bwa-mem2") or shutil.which("bwa") or "bwa-mem2"
    vc_bin = shutil.which("viral_consensus")
    has_ragtag = os.path.isfile(ragtag_bin) if ragtag_bin else False
    has_vc = vc_bin is not None
    print(f"  分支 C: ragtag={'✓' if has_ragtag else '✗ (仅BLASTN)'}  viral_consensus={'✓' if has_vc else '✗ (仅ragtag/BLASTN)'}  bwa={bwa_bin}")

    # ── BLAST DB: 优先用户指定, 否则自动构建 (cdhit_combined + 分支B成功补充) ──
    if blast_db:
        print(f"  分支 C: 使用指定 BLAST DB = {blast_db}")
    else:
        all_cluster_seqs = []
        seen_ids = set()
        # 主库: cdhit_combined.fasta (用户指定 → 自动搜索)
        if not cdhit_fa:
            for guess in [
                Path(work_dir).parent.parent / "04_CLUSTER" / "2_cdhit" / "cdhit_combined.fasta",
                Path(work_dir).parent / "04_CLUSTER" / "2_cdhit" / "cdhit_combined.fasta",
            ]:
                if guess.is_file():
                    cdhit_fa = str(guess)
                    break
            if not cdhit_fa:
                cdhit_glob = list(Path(work_dir).parent.parent.glob("*/2_cdhit/cdhit_combined.fasta"))
                if cdhit_glob:
                    cdhit_fa = str(cdhit_glob[0])
        cdhit_added = 0
        if cdhit_fa and Path(cdhit_fa).is_file():
            for rec in SeqIO.parse(str(cdhit_fa), "fasta"):
                if rec.id not in seen_ids:
                    all_cluster_seqs.append(rec)
                    seen_ids.add(rec.id)
                    cdhit_added += 1
            if cdhit_added:
                print(f"  分支 C: 加载 cdhit_combined → +{cdhit_added:,} 条 (主参考库)")
        # 补充: 分支 B 成功的基因组 (VSI reads 延伸后通过 CheckV/genus 的 contig)
        branch_b_pass = Path(work_dir) / "branch_b" / "branchB_pass.fasta"
        bb_added = 0
        if branch_b_pass.is_file() and branch_b_pass.stat().st_size > 0:
            for rec in SeqIO.parse(str(branch_b_pass), "fasta"):
                if rec.id not in seen_ids:
                    all_cluster_seqs.append(rec)
                    seen_ids.add(rec.id)
                    bb_added += 1
            if bb_added:
                print(f"  分支 C: 加载分支B成功基因组 → +{bb_added} 条")
        if not all_cluster_seqs:
            print("  分支 C: 无可用的序列构建 BLAST DB, 跳过"); return None, 0

        auto_db_dir = d / "auto_blast_db"
        auto_db_dir.mkdir(exist_ok=True)
        auto_db_fa = auto_db_dir / "cluster_centroids.fasta"
        SeqIO.write(all_cluster_seqs, auto_db_fa, "fasta")

        subprocess.run(["makeblastdb", "-in", str(auto_db_fa), "-dbtype", "nucl",
                        "-out", str(auto_db_dir / "auto_db"), "-title", "Cluster_Centroids"],
                       capture_output=True, check=False)
        blast_db = str(auto_db_dir / "auto_db")
        db_final_count = sum(1 for _ in SeqIO.parse(auto_db_fa, "fasta"))
        print(f"  分支 C: 自动构建 BLAST DB → {db_final_count} 条")

    # ── VirASCA 风格: 先 BLAST → 三档分类 → 属分组 → 分流 ──

    # Phase 1: BLASTN 所有 fail contigs (并行, 断点续传)
    blast_done_dir = d / "blast_results"; blast_done_dir.mkdir(exist_ok=True)
    contig_blast = {}  # {contig_id: {"best_acc": str, "best_pident": float, "best_qcov": float, "best_bitscore": float, "status": str}}
    blast_lock = threading.Lock()

    def _blast_one(rec):
        bid = rec.id
        out_tsv = blast_done_dir / f"{bid[:60]}.tsv"
        if not out_tsv.is_file():
            tmp_fa = blast_done_dir / f"{bid[:60]}.fa"
            SeqIO.write([rec], tmp_fa, "fasta")
            run_blastn(tmp_fa, blast_db, out_tsv, threads)
        if not out_tsv.is_file() or out_tsv.stat().st_size == 0:
            with blast_lock: contig_blast[bid] = {"status": "no_hit"}
            return
        hits = [l.strip() for l in open(out_tsv) if l.strip() and not l.startswith('#')]
        if not hits:
            with blast_lock: contig_blast[bid] = {"status": "no_hit"}
            return
        # 对每个 hit 计算 qcov, 排除自匹配, 选最佳
        best = None
        qlen = len(rec.seq)
        for line in hits:
            cols = line.split('\t')
            if len(cols) < 12: continue
            try:
                pid = float(cols[2]); alen = int(cols[3]); qlen_hit = int(cols[10])
                acc = cols[1]; bitscore = float(cols[9])
            except: continue
            if acc == bid: continue  # 自匹配
            qcov = alen / qlen if qlen > 0 else 0
            if best is None or bitscore > best["best_bitscore"]:
                best = {"best_acc": acc, "best_pident": pid, "best_qcov": qcov, "best_bitscore": bitscore}
        if best is None:
            with blast_lock: contig_blast[bid] = {"status": "no_hit"}
            return
        # 三档分类 (VirASCA: Completo/Intermediario/Incompleto)
        if best["best_qcov"] >= 0.98:
            best["status"] = "completo"
        elif best["best_pident"] >= 90:
            best["status"] = "intermediario"
        elif best["best_pident"] >= 60:
            best["status"] = "incompleto"
        else:
            best["status"] = "desconhecido"
        with blast_lock: contig_blast[bid] = best

    print(f"  分支 C: BLASTN {len(fail_records)} 个 contigs (并行 {jobs} jobs)...")
    if jobs > 1 and len(fail_records) > 1:
        with ThreadPoolExecutor(max_workers=min(jobs, len(fail_records))) as ex:
            list(tqdm(ex.map(_blast_one, fail_records), total=len(fail_records), desc="  BLASTN", unit="contig"))
    else:
        for rec in tqdm(fail_records, desc="  BLASTN", unit="contig"):
            _blast_one(rec)

    # 统计
    n_completo = sum(1 for v in contig_blast.values() if v.get("status") == "completo")
    n_inter = sum(1 for v in contig_blast.values() if v.get("status") == "intermediario")
    n_incomp = sum(1 for v in contig_blast.values() if v.get("status") == "incompleto")
    n_des = sum(1 for v in contig_blast.values() if v.get("status") == "desconhecido")
    n_nohit = sum(1 for v in contig_blast.values() if v.get("status") == "no_hit")

    # Phase 2: 按 Genus 分组 (从 05_Taxonomy), completo 直达 pass
    genus_map = {}
    if taxonomy_tsv and Path(taxonomy_tsv).is_file():
        genus_map = _load_taxonomy(taxonomy_tsv)
    genus_groups = defaultdict(set)
    no_genus_groups = set()  # 无属分类的组名集合
    completed_pass = []  # Completo contigs 直接通过
    for rec in fail_records:
        bi = contig_blast.get(rec.id, {})
        if bi.get("status") in ("no_hit", "desconhecido"):
            continue  # 无 BLASTN 命中或 identity<60%, 直接进分支D
        if bi.get("status") == "completo":
            # Completo: 直接放入 completed_pass, 不进入 genus 组
            completed_pass.append(rec)
            continue
        g = genus_map.get(rec.id, "")
        g = g.replace('g__', '').replace('G__', '')
        if not g:
            g = f"NOGENUS_{cmap.get(rec.id, rec.id)[:30]}"
            no_genus_groups.add(g)
        genus_groups[g].add(rec.id)
    if completed_pass:
        print(f"  分支 C: {len(completed_pass)} 条 Completo → 直达 pass (无需组内组装)")
    n_nogenus = sum(1 for g in genus_groups if g in no_genus_groups)
    n_valid = len(genus_groups) - n_nogenus
    print(f"  分支 C: {n_valid} 个属组 + {n_nogenus} 个无属组 (将直接进fail)")
    print(f"  分支 C 统计: Completo={n_completo}  Intermed={n_inter}  Incompleto={n_incomp}  Desconh={n_des}  NoHit={n_nohit}")

    # Phase 3: 断点续传
    pass_fa = d / "branchC_pass.fasta"
    ckpt_file = d / "branchC_checkpoint.txt"
    done_groups = set()
    if pass_fa.is_file() and pass_fa.stat().st_size > 0:
        for rec in SeqIO.parse(str(pass_fa), "fasta"):
            done_groups.add(rec.id)
    if ckpt_file.is_file():
        with open(ckpt_file) as cf:
            for line in cf: done_groups.add(line.strip())
    if done_groups:
        print(f"  分支 C: [resume] 已完成 {len(done_groups)} 个组")

    lock = threading.Lock()
    pass_records = list(completed_pass)  # Completo contigs 直达
    fail_group_ids = set()

    # NOGENUS 组标记为 fail (无属分类无法 blastdbcmd 提取参考)
    for g in list(genus_groups.keys()):
        if g.startswith("NOGENUS_"):
            fail_group_ids.add(g)
    pending_groups = {g: ms for g, ms in genus_groups.items()
                      if g not in done_groups and not g.startswith("NOGENUS_")}
    if fail_group_ids:
        n_nogenus_fail = sum(1 for g in fail_group_ids if g.startswith("NOGENUS_"))
        print(f"  分支 C: {n_nogenus_fail} 个无属组 → 直接进 fail (无参考可提取)")
    print(f"  分支 C: {len(pending_groups)} 个属组待处理")

    # ── 参考引导组装 ──
    def _ref_guided_assembly(ref_fa, member_ids, out_dir, rga_threads=8, log_fh=None):
        """member_ids: contig ID 列表, 从 cluster 反查样本 → 合并 reads → bwa-mem2+viral_consensus"""
        d = Path(out_dir); d.mkdir(parents=True, exist_ok=True)
        consensus_fa = d / "consensus.fa"
        log = log_fh or open(d / "rga.log", "a")
        # 断点续传
        if consensus_fa.is_file() and consensus_fa.stat().st_size > 100:
            return str(consensus_fa)
        sample_ids = set()
        for mid in member_ids:
            samples = _contig_cluster_samples(mid, clusters, fasta_info, max_vsi_samples, flye_sample_map=flye_sample_map)
            if samples:
                sample_ids.update(samples)
        if not sample_ids:
            log.write(f"[RGA] 无样本信息 (contigs: {list(member_ids)[:5]})\n"); log.flush()
            return None
        r1, r2, _ = _gather_cluster_reads(fastq_dir, list(sample_ids), d, prefix="rga")
        if not r1:
            log.write(f"[RGA] reads 未找到 (samples: {list(sample_ids)[:5]})\n"); log.flush()
            return None
        log.write(f"[RGA] bwa index {ref_fa}\n"); log.flush()
        subprocess.run([bwa_bin, "index", str(ref_fa)], stdout=log, stderr=subprocess.STDOUT, check=False)
        mapped_bam = d / "mapped.bam"
        bwa_args = [bwa_bin, "mem", "-t", str(rga_threads), str(ref_fa), str(r1)]
        if r2 and Path(r2).is_file() and Path(r2).stat().st_size > 0:
            bwa_args.append(str(r2))
        log.write(f"[RGA] bwa mem + samtools sort\n"); log.flush()
        p1 = subprocess.Popen(bwa_args, stdout=subprocess.PIPE, stderr=log)
        p2 = subprocess.Popen(["samtools", "sort", "-o", str(mapped_bam), "-"],
                              stdin=p1.stdout, stdout=log, stderr=subprocess.STDOUT)
        if p1.stdout: p1.stdout.close()
        p2.communicate()
        if not mapped_bam.is_file() or mapped_bam.stat().st_size < 100:
            log.write("[RGA] mapped.bam empty or missing\n"); log.flush()
            return None
        log.write("[RGA] samtools index\n"); log.flush()
        subprocess.run(["samtools", "index", str(mapped_bam)], stdout=log, stderr=subprocess.STDOUT, check=False)
        log.write("[RGA] viral_consensus\n"); log.flush()
        subprocess.run([vc_bin if has_vc else "viral_consensus",
                 "-i", str(mapped_bam),
                 "-r", str(ref_fa),
                 "-o", str(consensus_fa),
                 "-d", "1",
                 "-f", "0.5",
                 "-q", "20",
                 "--ambig", "N"],
                stdout=log, stderr=subprocess.STDOUT, check=False)
        if consensus_fa.is_file() and consensus_fa.stat().st_size > 100:
            log.write(f"[RGA] done → {consensus_fa}\n"); log.flush()
            return str(consensus_fa)
        log.write("[RGA] consensus.fa empty or missing\n"); log.flush()
        return None

    def _do_group(group_id):
        """输入: genus 名. 从 contig_blast 取所有成员的 BLASTN 结果, 选参考, 分流."""
        member_ids = genus_groups[group_id]
        sub_dir = d / f"group_{group_id[:50].replace('/','_').replace(' ','_')}"
        sub_dir.mkdir(exist_ok=True)
        group_log = open(sub_dir / "group.log", "w")
        group_log.write(f"# Group: {group_id}\n# Members: {len(member_ids)}\n\n")
        group_log.flush()

        # 收集同属 contigs 的 BLASTN 最佳 hit
        group_hits = {}  # {acc: {"total_pident": 0, "total_bitscore": 0, "count": 0, "contigs": set()}}
        group_contigs = []
        for mid in member_ids:
            bi = contig_blast.get(mid, {})
            if bi.get("status") in ("no_hit", "desconhecido", None):
                continue
            acc = bi["best_acc"]
            if acc not in group_hits:
                group_hits[acc] = {"total_pident": 0, "total_bitscore": 0, "count": 0, "contigs": set()}
            group_hits[acc]["total_pident"] += bi["best_pident"]
            group_hits[acc]["total_bitscore"] += bi["best_bitscore"]
            group_hits[acc]["count"] += 1
            group_hits[acc]["contigs"].add(mid)
            # 收集 contig 序列
            if mid in centroids_all:
                group_contigs.append(centroids_all[mid])
            else:
                match = [r for r in fail_records if r.id == mid]
                if match: group_contigs.append(match[0])

        if not group_hits:
            group_log.write("[SKIP] no valid BLASTN hits\n"); group_log.close()
            fail_group_ids.add(group_id)
            return

        # 选最佳参考 (平均 pident 最高, 次选 bitscore)
        best_acc = max(group_hits, key=lambda a: (group_hits[a]["total_pident"] / group_hits[a]["count"], group_hits[a]["total_bitscore"]))
        best_info = group_hits[best_acc]
        avg_pident = best_info["total_pident"] / best_info["count"]
        n_contigs = best_info["count"]
        print(f"  [{group_id[:40]}] | {n_contigs} contigs → ref={best_acc} | avg_pident={avg_pident:.1f}%", flush=True)

        # 提取参考序列
        ref_fa = sub_dir / "reference.fa"
        if not ref_fa.is_file():
            top_acc_clean = best_acc.replace('ref|', '').replace('|', ' ')
            group_log.write(f"[blastdbcmd] extracting {best_acc}\n"); group_log.flush()
            result = subprocess.run([blastdbcmd_bin, "-db", str(blast_db), "-entry", best_acc,
                                     "-out", str(ref_fa)], capture_output=True, text=True)
            if result.returncode != 0:
                result = subprocess.run([blastdbcmd_bin, "-db", str(blast_db), "-entry",
                                         top_acc_clean, "-out", str(ref_fa)], capture_output=True, text=True)
                if result.returncode != 0 or not ref_fa.is_file() or ref_fa.stat().st_size <= 100:
                    group_log.write(f"[FAIL] blastdbcmd: {best_acc} → {result.stderr}\n"); group_log.flush()
                    print(f"  [FAIL] blastdbcmd: {best_acc}", flush=True)
                    fail_group_ids.add(group_id)
                    group_log.close()
                    return

        # 路径分流 (VirASCA: Intermediario/Incompleto, Completo已在Phase2直达pass)
        final_fa = None; method = None
        query_fa = sub_dir / "query_contigs.fa"
        SeqIO.write(group_contigs, query_fa, "fasta")
        scaffold_fa = sub_dir / "ragtag_out" / "ragtag.scaffold.fasta"
        cv_out = sub_dir / "cv" / "completeness.tsv"

        # Intermediario (avg_pident >= 90%): 参考引导组装
        if avg_pident >= 90.0:
            group_log.write(f"[method] Intermediario → RGA\n"); group_log.flush()
            consensus = _ref_guided_assembly(ref_fa, member_ids, sub_dir / "rga", log_fh=group_log)
            if consensus:
                final_fa = Path(consensus)
                method = "RGA"
            else:
                method = "RGA-failed"

        # Incompleto (avg_pident 60-90%): ragtag (>=2 contigs) 或 RGA (<2 contigs)
        elif n_contigs >= 2:
            group_log.write(f"[method] Incompleto ({n_contigs} contigs) → RagTag\n"); group_log.flush()
            if not scaffold_fa.is_file():
                ragtag_out = sub_dir / "ragtag_out"
                ragtag_out.mkdir(exist_ok=True)
                subprocess.run([ragtag_bin, "scaffold", str(ref_fa), str(query_fa),
                              "-o", str(ragtag_out), "-t", str(threads), "-w", "-u"],
                             stdout=group_log, stderr=subprocess.STDOUT, check=False)
            if scaffold_fa.is_file() and scaffold_fa.stat().st_size > 0:
                final_fa = scaffold_fa
                method = "ragtag"
        else:
            # Incompleto 单 contig: 回退 RGA
            group_log.write(f"[method] Incompleto (1 contig) → RGA-solo\n"); group_log.flush()
            consensus = _ref_guided_assembly(ref_fa, member_ids, sub_dir / "rga", log_fh=group_log)
            if consensus:
                final_fa = Path(consensus)
                method = "RGA-solo"

        if final_fa is None:
            group_log.write(f"[FAIL] {method}\n"); group_log.flush()
            print(f"  [FAIL] {group_id[:40]} | {method}", flush=True)
            fail_group_ids.add(group_id)
            group_log.close()
            return

        # 输出 + CheckV 标注
        final_size = final_fa.stat().st_size if final_fa.is_file() else 0
        comp_val = "NA"
        if checkv_db and final_fa.is_file() and final_size > 0:
            if not cv_out.is_file() or cv_out.stat().st_mtime < final_fa.stat().st_mtime:
                run_checkv(final_fa, sub_dir / "cv", checkv_db, threads)
            if cv_out.is_file():
                try:
                    df = pd.read_csv(cv_out, sep='\t')
                    if 'aai_completeness' in df.columns and len(df) > 0:
                        cv = df['aai_completeness'].iloc[0]
                        comp_val = f"{cv:.1f}%" if not pd.isna(cv) else "NA"
                except: pass
        group_log.write(f"[DONE] method={method} size={final_size}bp CheckV={comp_val}\n"); group_log.flush()
        print(f"  [DONE] {group_id[:40]} | method={method} | size={final_size}bp | CheckV={comp_val}", flush=True)
        group_log.close()
        with lock:
            for srec in SeqIO.parse(str(final_fa), "fasta"):
                srec.id = group_id
                srec.description = f"method={method} CheckV={comp_val}"
                pass_records.append(srec)
            with open(ckpt_file, "a") as cf: cf.write(group_id + "\n")

    # Phase 4: 执行
    if jobs > 1 and len(pending_groups) > 1:
        with ThreadPoolExecutor(max_workers=min(jobs, len(pending_groups))) as ex:
            list(tqdm(ex.map(_do_group, pending_groups), total=len(pending_groups),
                     desc="  分支 C", unit="group"))
    else:
        for gid in tqdm(list(pending_groups), desc="  分支 C", unit="group"):
            _do_group(gid)

    if pass_records:
        SeqIO.write(pass_records, pass_fa, "fasta")

    # 输出失败的 contigs → 供分支 D 使用
    fail_fa_out = d / "branchC_fail.fasta"
    fail_contig_ids = set()
    for gid in fail_group_ids:
        if gid in genus_groups:
            fail_contig_ids |= genus_groups[gid]
    # 加入 no_hit + desconhecido contigs
    for rec in fail_records:
        bi = contig_blast.get(rec.id, {})
        if bi.get("status") in ("no_hit", "desconhecido"):
            fail_contig_ids.add(rec.id)
    fail_contigs = [r for r in fail_records if r.id in fail_contig_ids]
    if fail_contigs:
        SeqIO.write(fail_contigs, fail_fa_out, "fasta")

    n_pass = len(pass_records)
    n_fail = len(fail_contigs)
    print(f"  分支 C: pass={n_pass}  fail={n_fail}")

    # 写出 genus→contig 映射 (供最终报告追溯)
    map_tsv = d / "branchC_genus_map.tsv"
    with open(map_tsv, "w") as mf:
        mf.write("genus_group\tcontig_ids\n")
        for gid, mids in genus_groups.items():
            mf.write(f"{gid}\t{','.join(mids)}\n")

    if n_pass == 0:
        return None, str(fail_fa_out) if fail_contigs else None, 0, n_fail
    return str(pass_fa), str(fail_fa_out) if fail_contigs else None, n_pass, n_fail


# ══════════════════════════════════════════════════════════════
# 分支 D: genus_len 属水平长度拯救
# ══════════════════════════════════════════════════════════════

def _load_taxonomy(tax_tsv):
    """加载 taxonomy TSV, 返回 {contig_id: genus}"""
    mapping = {}
    if not tax_tsv or not Path(tax_tsv).is_file():
        return mapping
    df = pd.read_csv(tax_tsv, sep='\t')
    for _, row in df.iterrows():
        cid = str(row.get('contig_id', ''))
        genus = str(row.get('Genus', '')).strip()
        if cid and genus and genus not in ('NA', '', 'nan'):
            # genus 值可能带前缀 g__ 或 G_, 统一去掉
            genus = genus.replace('g__', '').replace('G__', '').replace('g_', '').replace('G_', '')
            mapping[cid] = genus
    print(f"  taxonomy 加载: {len(mapping)} 个 contig→genus 映射")
    return mapping


def _load_genus_len(genus_len_path, format_prefix='g__'):
    """加载 genus_len 文件, 返回 {genus: avg_length}"""
    mapping = {}
    if not genus_len_path or not Path(genus_len_path).is_file():
        return mapping
    with open(genus_len_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('genus'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            genus = parts[0].strip()
            try:
                avg_len = float(parts[1])
            except ValueError:
                continue
            # 去掉前缀 (g__), 也保留原始形式
            genus_clean = genus.replace('g__', '').replace('G__', '')
            mapping[genus_clean] = avg_len
            if genus_clean != genus:
                mapping[genus] = avg_len
    print(f"  genus_len 加载: {len(mapping)} 个 genus")
    return mapping


def branch_d(fail_fa, tax_tsv, genus_len_path, output_dir,
             genus_tol=0.85, min_len=500):
    """分支 D: genus_len 属水平长度拯救

    对 ABC 都失败的 contig, 如果其长度 ≥ 同属平均长度的 genus_tol (默认85%),
    则认为其具有属水平完整长度, 拯救为 HQ。

    Args:
        fail_fa: 最终失败的 FASTA 路径 (来自 branch C 或 B)
        tax_tsv: taxonomy 分类 TSV (final_integrated_classification.tsv)
        genus_len_path: genus→avg_length 参考文件
        output_dir: 输出根目录 (通常为 08_Rescue/Plant/)
        genus_tol: 容忍度 (默认 0.85, 即 contig >= 85% 属平均长度则拯救)
        min_len: 最小长度阈值 bp (低于此值不参与 genus 拯救)
    Returns:
        (pass_fa, pass_count): 拯救的 FASTA 路径 & 数量
    """
    d = Path(output_dir) / "branch_d"; d.mkdir(parents=True, exist_ok=True)

    if not fail_fa or not Path(fail_fa).is_file():
        print("  分支 D: 无失败序列"); return None, 0

    if not tax_tsv or not Path(tax_tsv).is_file():
        print("  分支 D: 无 taxonomy 文件, 跳过"); return None, 0

    if not genus_len_path or not Path(genus_len_path).is_file():
        print("  分支 D: 无 genus_len 文件, 跳过"); return None, 0

    fail_records = list(SeqIO.parse(str(fail_fa), "fasta"))
    total = len(fail_records)
    if not fail_records:
        print("  分支 D: 无失败序列"); return None, 0

    # 加载 taxonomy 和 genus_len
    tax_map = _load_taxonomy(tax_tsv)
    genus_len_map = _load_genus_len(genus_len_path)

    rescued = []
    no_tax = 0
    no_genus_ref = 0
    too_short = 0
    below_threshold = 0

    for rec in fail_records:
        cid = rec.id
        contig_len = len(rec.seq)

        # 长度过滤
        if contig_len < min_len:
            too_short += 1
            continue

        # 查 taxonomy → genus
        genus = tax_map.get(cid)
        if not genus:
            no_tax += 1
            continue

        # 查 genus → avg_length
        avg_len = genus_len_map.get(genus)
        if not avg_len or avg_len <= 0:
            no_genus_ref += 1
            continue

        # 计算 fraction
        fraction = contig_len / avg_len
        if fraction >= genus_tol:
            rec.description = f"genus_rescued genus={genus} len={contig_len}bp avg={avg_len:.0f}bp fraction={fraction:.2f}"
            rescued.append(rec)
        else:
            below_threshold += 1

    pass_fa = d / "branchD_pass.fasta"
    report_tsv = d / "branchD_report.tsv"
    fail_fa_out = d / "branchD_fail.fasta"
    fail_records_out = []

    if rescued:
        SeqIO.write(rescued, str(pass_fa), "fasta")

    # 写出详细报告（无论是否有成功拯救）
    with open(report_tsv, "w") as rf:
        rf.write("contig_id\tlength\tgenus\tref_avg_len\tfraction\tstatus\n")
        for rec in fail_records:
            cid = rec.id; clen = len(rec.seq)
            g = tax_map.get(cid, "")
            al = genus_len_map.get(g, 0) if g else 0
            frac = clen / al if al > 0 else 0
            if clen < min_len:
                status = "too_short"
            elif not g:
                status = "no_taxonomy"
            elif not al:
                status = "no_genus_ref"
            elif frac >= genus_tol:
                status = "rescued"
            else:
                status = "below_threshold"
                fail_records_out.append(rec)
            rf.write(f"{cid}\t{clen}\t{g}\t{al:.0f}\t{frac:.4f}\t{status}\n")

    # 写出失败序列 FASTA（供调试）
    if fail_records_out:
        SeqIO.write(fail_records_out, str(fail_fa_out), "fasta")

    print(f"  分支 D: rescue={len(rescued)}  no_tax={no_tax}  "
          f"no_ref={no_genus_ref}  short={too_short}  below_thr={below_threshold}")
    return str(pass_fa) if rescued else None, len(rescued)


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
# 最终报告
# ══════════════════════════════════════════════════════════════

def _write_rescue_report(out_dir, centroids_records, clusters,
                         fa_a_pass, fa_b_pass, fa_c_pass, fa_d_pass,
                         cnt_a, cnt_b, cnt_c, cnt_d, n_final,
                         taxonomy_tsv=None, genus_lens_path=None, min_vsi_len=2000):
    """生成 rescue 最终追踪报告 TSV + Markdown"""
    d = Path(out_dir); d.mkdir(parents=True, exist_ok=True)

    def _load_pass_ids(fa_path):
        ids = {}
        if fa_path and os.path.isfile(fa_path):
            for rec in SeqIO.parse(fa_path, "fasta"):
                ids[rec.id] = {"length": len(rec.seq), "desc": rec.description}
        return ids

    pass_a = _load_pass_ids(fa_a_pass)
    pass_b = _load_pass_ids(fa_b_pass)
    pass_c = _load_pass_ids(fa_c_pass)
    pass_d = _load_pass_ids(fa_d_pass)

    # 加载 taxonomy + genus_lens (补充属平均长度信息)
    tax_map = {}
    genus_len_map = {}
    if taxonomy_tsv and Path(taxonomy_tsv).is_file():
        tax_map = _load_taxonomy(taxonomy_tsv)
    if genus_lens_path and Path(genus_lens_path).is_file():
        genus_len_map = _load_genus_len(genus_lens_path)

    # 分支 C 的 genus 组 → 原始 contig 映射
    genus_to_contigs = {}
    map_tsv = d / "branch_c" / "branchC_genus_map.tsv"
    if map_tsv.is_file():
        with open(map_tsv) as mf:
            mf.readline()
            for line in mf:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    genus_to_contigs[parts[0]] = parts[1].split(',')

    contig_to_genus = {}
    for gid, cids in genus_to_contigs.items():
        for cid in cids:
            contig_to_genus[cid] = gid

    # 写 TSV 报告
    tsv_path = d / "rescue_report.tsv"
    # 失败原因统计
    fail_reasons = {}

    def _lookup_genus(cid):
        g = tax_map.get(cid, "")
        return g.replace('g__', '').replace('G__', '')

    with open(tsv_path, "w") as tf:
        tf.write("contig_id\tinput_len\tfinal_len\tbranch\tmethod\tgenus\tgenus_avg_len\tpct_of_genus\tcheckv\tnote\n")
        for rec in centroids_records:
            cid = rec.id; clen = len(rec.seq)
            genus = _lookup_genus(cid)
            gal = genus_len_map.get(genus, 0)
            pct = f"{clen/gal*100:.1f}%" if gal > 0 else "-"

            if cid in pass_a:
                info = pass_a[cid]
                tf.write(f"{cid}\t{clen}\t{info['length']}\tA\tcheckv\t{genus}\t{gal:.0f}\t{pct}\t≥90%\tprotein_complete\n")

            elif cid in pass_b:
                info = pass_b[cid]
                flen = info['length']
                pct2 = f"{flen/gal*100:.1f}%" if gal > 0 else "-"
                # 从 desc 中提取 CheckV 值
                cv = "-"
                desc = info.get('desc', '')
                for part in desc.split():
                    if part.startswith('checkv='):
                        cv = part.replace('checkv=', '')
                        break
                rescue_type = "genus_rescued" if "genus_rescued" in desc else "vsi_extended"
                tf.write(f"{cid}\t{clen}\t{flen}\tB\t{rescue_type}\t{genus}\t{gal:.0f}\t{pct2}\t{cv}\t{desc}\n")

            elif cid in pass_d:
                info = pass_d[cid]
                flen = info['length']
                pct2 = f"{flen/gal*100:.1f}%" if gal > 0 else "-"
                desc = info.get('desc', '')
                tf.write(f"{cid}\t{clen}\t{flen}\tD\tgenus_len\t{genus}\t{gal:.0f}\t{pct2}\t-\t{desc}\n")

            elif cid in contig_to_genus:
                gid = contig_to_genus[cid]
                if gid in pass_c:
                    info = pass_c[gid]
                    flen = info['length']
                    desc = info.get('desc', '')
                    method = '?'; cv = '-'
                    for part in desc.split():
                        if part.startswith('method='):
                            method = part.replace('method=', '')
                        elif part.startswith('CheckV='):
                            cv = part.replace('CheckV=', '')
                    gc = _lookup_genus(gid) if genus == _lookup_genus(gid) else genus
                    pct2 = f"{flen/gal*100:.1f}%" if gal > 0 else "-"
                    tf.write(f"{cid}\t{clen}\t{flen}\tC\t{method}\t{gid}\t{gal:.0f}\t{pct2}\t{cv}\t组组装成功\n")
                else:
                    reason = "genus_group_failed"
                    fail_reasons[reason] = fail_reasons.get(reason, 0) + 1
                    tf.write(f"{cid}\t{clen}\t-\tfail\t-\t{gid}\t-\t-\t-\tgenus组处理失败\n")

            else:
                # 分类失败原因
                reason = "unknown"
                if clen < min_vsi_len:
                    reason = "too_short"
                elif not genus:
                    reason = "no_taxonomy"
                elif gal > 0 and clen / gal >= 0.85:
                    reason = "should_rescue_but_skipped"
                else:
                    reason = "below_all_thresholds"
                fail_reasons[reason] = fail_reasons.get(reason, 0) + 1
                note = {
                    "too_short": f"<{min_vsi_len}bp",
                    "no_taxonomy": "无属分类",
                    "should_rescue_but_skipped": "长度达标但未被拯救",
                    "below_all_thresholds": "未通过任何分支",
                    "unknown": "未拯救",
                }.get(reason, reason)
                tf.write(f"{cid}\t{clen}\t-\tfail\t-\t{genus}\t{gal:.0f}\t{pct}\t-\t{note}\n")

    rescued = cnt_a + cnt_b + cnt_c + cnt_d
    total = len(centroids_records)
    failed = total - rescued

    # 失败原因描述
    reason_labels = {
        "too_short": f"过短 (<{min_vsi_len}bp)",
        "no_taxonomy": "无属分类信息",
        "should_rescue_but_skipped": "长度≥85%属平均但未拯救",
        "below_all_thresholds": "未通过任何分支判定",
        "genus_group_failed": "属组处理失败 (blastdbcmd/RGA/ragtag)",
        "unknown": "其他原因",
    }

    md_path = d / "rescue_summary.md"
    with open(md_path, "w") as mf:
        mf.write("# 病毒基因组拯救报告\n\n")
        mf.write(f"## 分支统计\n\n")
        mf.write(f"| 分支 | 策略 | 拯救数 | 占比 |\n")
        mf.write(f"|------|------|--------|------|\n")
        mf.write(f"| A | CheckV ≥90% (蛋白完整) | {cnt_a} | {cnt_a/total*100:.1f}% |\n")
        mf.write(f"| B | VSI reads延伸 + genus回退 | {cnt_b} | {cnt_b/total*100:.1f}% |\n")
        mf.write(f"| C | BLASTN + RGA/ragtag | {cnt_c} | {cnt_c/total*100:.1f}% |\n")
        mf.write(f"| D | genus_len 属水平兜底 | {cnt_d} | {cnt_d/total*100:.1f}% |\n")
        mf.write(f"| **合计** | | **{rescued}** | **{rescued/total*100:.1f}%** |\n")
        mf.write(f"| 未拯救 | | {failed} | {failed/total*100:.1f}% |\n\n")

        if fail_reasons:
            mf.write(f"## 未拯救原因分布\n\n")
            mf.write(f"| 原因 | 数量 | 占比 |\n")
            mf.write(f"|------|------|------|\n")
            for reason, count in sorted(fail_reasons.items(), key=lambda x: -x[1]):
                label = reason_labels.get(reason, reason)
                mf.write(f"| {label} | {count} | {count/failed*100:.1f}% |\n")
            mf.write("\n")

        mf.write(f"最终无冗余 vOTU: **{n_final}**\n\n")
        mf.write(f"详细: [{tsv_path.name}]({tsv_path.name})\n")

    print(f"  报告: {tsv_path}")
    print(f"  摘要: {md_path}")


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
    p.add_argument("--checkv-threshold", type=float, default=90.0, help="CheckV completeness 通过阈值 (默认90, Plant建议80)")
    p.add_argument("--taxonomy-tsv", default=None, help="05_Taxonomy 分类结果 final_integrated_classification.tsv (分支D genus_len 拯救)")
    p.add_argument("--cdhit-fa", default=None, help="cdhit_combined.fasta 路径 (供分支C自动建库, 默认自动搜索 04_CLUSTER/)")
    p.add_argument("--genus-len", default=None, help="genus_len 参考文件 (属平均长度, 分支D用)")
    p.add_argument("--genus-tolerance", type=float, default=0.85, help="分支D genus 拯救容忍度 (默认0.85, 即contig≥85%%属平均长度则拯救)")
    p.add_argument("--threads", "-t", type=int, default=64)
    p.add_argument("--jobs", "-j", type=int, default=4, help="Virseqimprover 并行数")
    p.add_argument("--ani", type=float, default=0.95, help="最终 vclust ANI")
    p.add_argument("--qcov", type=float, default=0.85, help="最终 vclust QCOV")
    p.add_argument("--resume", action="store_true", help="断点续传")
    p.add_argument("--flye-sample-map", default=None,
                   help="Flye 共组装 contig→样本 映射 JSON (merge 阶段产出, 供 VSI/RGA reads 回退)")
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

    # 加载 Flye contig→样本 映射 (merge 阶段产出)
    flye_sample_map = None
    if args.flye_sample_map and os.path.isfile(args.flye_sample_map):
        with open(args.flye_sample_map) as fm:
            flye_sample_map = json.load(fm)
        print(f"  flye_sample_map: {len(flye_sample_map)} 条 Flye contig 映射")

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
        fa_a_pass, fa_a_fail, cnt_a, cnt_a_fail = branch_a(str(tmp_centroids), out, args.checkv_db, args.threads, args.jobs, args.checkv_threshold)

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
        # 加载属分类信息 (供 VSI 的 genus_avg_len 备选截止)
        genus_map = {}
        genus_lens = {}
        if getattr(args, 'taxonomy_tsv', None):
            genus_map = _load_taxonomy(args.taxonomy_tsv)
        if getattr(args, 'genus_len', None):
            genus_lens = _load_genus_len(args.genus_len)
        fa_b_pass, fa_b_fail, cnt_b, cnt_b_fail = branch_b(fa_a_fail, args.fastq_dir, out, args.checkv_db, args.threads, args.jobs, args.virseqimprover_path, args.salmon_bin, clusters, fasta_info, args.max_vsi_samples, args.min_vsi_len, args.checkv_threshold, genus_map=genus_map, genus_lens=genus_lens, flye_sample_map=flye_sample_map)

    # 4. 分支 C: BLASTN (纯比对)
    print("\n── Step 3: 分支 C (BLASTN) ──")
    fa_c_pass = out / "branch_c" / "branchC_pass.fasta"
    fa_c_fail = out / "branch_c" / "branchC_fail.fasta"
    cnt_c = 0; cnt_c_fail = 0
    if args.resume and fa_c_pass.is_file() and fa_c_pass.stat().st_size > 0:
        print("  [RESUME] 分支 C 已有结果, 跳过")
        cnt_c = sum(1 for _ in SeqIO.parse(fa_c_pass, "fasta"))
        cnt_c_fail = sum(1 for _ in SeqIO.parse(fa_c_fail, "fasta")) if fa_c_fail.is_file() else 0
        fa_c_pass, fa_c_fail = str(fa_c_pass), str(fa_c_fail) if fa_c_fail.is_file() else None
    else:
        fail_for_c = fa_b_fail if (fa_b_fail and Path(fa_b_fail).is_file() and Path(fa_b_fail).stat().st_size > 0) else fa_a_fail
        fa_c_pass, fa_c_fail, cnt_c, cnt_c_fail = branch_c(fail_for_c, args.fastq_dir, out, args.checkv_db,
            getattr(args, 'blast_db', None),
            args.threads, args.jobs,
            clusters=clusters,
            clustermap_tsv=args.clusters_tsv,
            centroids_fa=args.centroids,
            split_dir=args.split_dir,
            fasta_info=fasta_info,
            taxonomy_tsv=getattr(args, 'taxonomy_tsv', None),
            cdhit_fa=getattr(args, 'cdhit_fa', None),
            flye_sample_map=flye_sample_map)

    # 5. 分支 D: genus_len 属水平长度拯救 (ABC 全失败后的最终兜底)
    print("\n── Step 4: 分支 D (genus_len) ──")
    fa_d_pass = out / "branch_d" / "branchD_pass.fasta"
    cnt_d = 0
    if args.resume and fa_d_pass.is_file() and fa_d_pass.stat().st_size > 0:
        print("  [RESUME] 分支 D 已有结果, 跳过")
        cnt_d = sum(1 for _ in SeqIO.parse(fa_d_pass, "fasta"))
        fa_d_pass = str(fa_d_pass)
    else:
        # 取最后一个可用的 fail_fa: C > B > A
        fail_for_d = fa_c_fail
        if not fail_for_d or not Path(fail_for_d).is_file():
            fail_for_d = fa_b_fail if (fa_b_fail and Path(fa_b_fail).is_file()) else None
        if not fail_for_d:
            fail_for_d = fa_a_fail if (fa_a_fail and Path(fa_a_fail).is_file()) else None
        if fail_for_d and Path(fail_for_d).is_file():
            fa_d_pass, cnt_d = branch_d(fail_for_d, args.taxonomy_tsv, args.genus_len,
                                         out, args.genus_tolerance)
        else:
            print("  分支 D: 无失败序列可处理")

    # 6. 合并
    print("\n── Step 5: 合并 ──")
    d4 = out / "merged"; d4.mkdir(exist_ok=True)
    merged = d4 / "all_HQ.fasta"
    centroids_final = out / "centroids" / "final_centroids.fasta"

    if args.resume and centroids_final.is_file() and centroids_final.stat().st_size > 0:
        print("  [RESUME] Step 5+6 已有最终结果, 跳过")
        total_m = sum(1 for _ in SeqIO.parse(merged, "fasta")) if merged.is_file() else cnt_a + cnt_b + cnt_c + cnt_d
        n_final = sum(1 for _ in SeqIO.parse(centroids_final, "fasta"))
    else:
        with open(merged, "w") as mf:
            for fp in [fa_a_pass, fa_b_pass, fa_c_pass, fa_d_pass]:
                if fp and os.path.isfile(fp):
                    with open(fp) as inf:
                        mf.write(inf.read())
        total_m = sum(1 for _ in SeqIO.parse(merged, "fasta"))
        print(f"  A:{cnt_a}  B:{cnt_b}  C:{cnt_c}  D:{cnt_d}  →  {total_m} 条  →  {merged}")

        if total_m < 2:
            print("  [SKIP] <2 条, 跳过最终 vclust")
            centroids_final.parent.mkdir(parents=True, exist_ok=True)
            with open(centroids_final, "w") as cf, open(merged) as mf2:
                cf.write(mf2.read())
            n_final = total_m
        else:
            # 7. 最终 vclust 去重
            print("\n── Step 6: vclust 最终去重 ──")
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

    # ── 生成最终追踪报告 ──
    print("\n── Step 7: 生成报告 ──")
    _write_rescue_report(out, centroids_records, clusters,
                          fa_a_pass, fa_b_pass, fa_c_pass, fa_d_pass,
                          cnt_a, cnt_b, cnt_c, cnt_d, n_final,
                          taxonomy_tsv=getattr(args, 'taxonomy_tsv', None),
                          genus_lens_path=getattr(args, 'genus_len', None),
                          min_vsi_len=args.min_vsi_len)

    elapsed = (datetime.now() - start).total_seconds()
    print(f"\n{'=' * 60}")
    print(f"  分支 A (CheckV):      {cnt_a:,}")
    print(f"  分支 B (VSI):         {cnt_b:,}")
    print(f"  分支 C (BLASTN):      {cnt_c:,}")
    print(f"  分支 D (genus_len):   {cnt_d:,}")
    print(f"  最终无冗余:           {n_final:,}")
    print(f"  最终输出:             {centroids_final}")
    print(f"  耗时:                 {elapsed / 60:.1f} min")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
