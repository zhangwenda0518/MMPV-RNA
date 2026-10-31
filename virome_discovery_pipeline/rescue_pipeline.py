#!/usr/bin/env python3
"""
rescue_pipeline.py -- 病毒基因组两支路级联拯救 v3.2
==================================================

独立拯救脚本: 接收已聚类的 centroids + clusters + 拆分文件, 直接执行拯救。
不重新聚类, 不重新过滤。聚类由 cluster_pipeline.py 完成, rescue 只做拯救。

流程:
  分支 A: CheckV 并行评估 → completeness>90% 直通; 失败 → 分支 B
  分支 B: Virseqimprover reads 延伸 (cluster 内多样本聚合) → CheckV 复评
          (CheckV=NA/不达标但长度达属均 85% 时回退拯救, 节段属按范围判定)
  A/B 均失败 → 记为未拯救, 不入 HQ 目录

分支 C (BLASTN 参考引导重建) 与分支 D (genus_len 兜底) 已于 2026-10-11 移除归档
(archive/20261011_rescue_branch_cd/): C 的重建路径因自动库缺 -parse_seqids 全部失败,
实际产出退化为"片段被参考包含即直通" (OneKP 9212/9212 直通, <1kb 占 11%), 与高质量
序列目标相悖; D 在全部实测数据集零产出且与 B 的属长回退重复。详见归档目录 README。

依赖: checkv, Virseqimprover.py
"""

import argparse, subprocess, sys, os, threading, shutil, json, re
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根): 目录名随 MMPV_IO_LAYOUT 解析
# (编排器已 normalize 环境变量, 子进程导入本模块时快照即正确布局)
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import layout_dirs as _layout_dirs
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT"))


def _safe_tqdm(iterable, total=None, desc="", unit=""):
    """TTY 下用 tqdm 进度条；非 TTY (日志重定向) 下用 print 定期打印进度。
    解决 tqdm 在非 TTY 环境下 \r 刷新导致日志里进度条冻结的问题。"""
    import sys, time
    if sys.stderr.isatty():
        return tqdm(iterable, total=total, desc=desc, unit=unit, file=sys.stderr)

    if total is None:
        try:
            total = len(iterable)
        except TypeError:
            total = None

    class _LogIter:
        def __init__(self, it, total, desc):
            self._it = iter(it)
            self._total = total
            self._desc = desc
            self._n = 0
            self._last_print = 0.0
            self._last_pct = -1
            self._t0 = time.time()

        def __iter__(self):
            return self

        def __next__(self):
            try:
                item = next(self._it)
            except StopIteration:
                raise
            self._n += 1
            now = time.time()
            pct = int(self._n / self._total * 100) if self._total else 0
            if now - self._last_print >= 30 or pct != self._last_pct:
                elapsed = int(now - self._t0)
                if self._total:
                    print(f"{self._desc}: {self._n}/{self._total} ({pct}%) [{elapsed}s]", flush=True)
                else:
                    print(f"{self._desc}: {self._n} [{elapsed}s]", flush=True)
                self._last_print = now
                self._last_pct = pct
            return item

    return _LogIter(iterable, total, desc)


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
    """解析 completeness.tsv (aai_completeness 列); 返回 (pass_ids, fail_ids, skip_ids, meta)
    skip_ids: aai_completeness=NA 的 contig
    meta: {contig_id: {aai_confidence, aai_expected_length}} 供后续智能判断"""
    if not Path(qs_path).is_file():
        return set(), set(), set(), {}
    df = pd.read_csv(qs_path, sep='\t')
    comp_col = 'aai_completeness' if 'aai_completeness' in df.columns else 'completeness'
    pass_ids = set()
    fail_ids = set()
    skip_ids = set()
    meta = {}
    for _, row in df.iterrows():
        cid = str(row.get('contig_id', ''))
        val = row.get(comp_col, 0)
        meta[cid] = {
            'confidence': str(row.get('aai_confidence', '')).strip(),
            'expected_len': float(row.get('aai_expected_length', 0) or 0),
        }
        if pd.isna(val) or str(val).strip() in ('NA', 'Not-determined', ''):
            skip_ids.add(cid)
        elif float(val) >= threshold:
            pass_ids.add(cid)
        else:
            fail_ids.add(cid)
    return pass_ids, fail_ids, skip_ids, meta


# ══════════════════════════════════════════════════════════════
# Virseqimprover 调用
# ══════════════════════════════════════════════════════════════

def _read_final_length(out_dir):
    """从 VSI 的 output-log.txt 读 final_length (无则 None)"""
    log = Path(out_dir) / "output-log.txt"
    if not log.is_file():
        return None
    try:
        for line in log.read_text(errors="ignore").splitlines():
            m = re.search(r"final_length=(\d+)", line)
            if m:
                return int(m.group(1))
    except Exception:
        pass
    return None


def _read_seq_len(fa):
    """读 fasta 第一条序列长度 (无则 None)"""
    try:
        for rec in SeqIO.parse(fa, "fasta"):
            return len(rec.seq)
    except Exception:
        return None


def _extract_vsi_error(log_path, out_dir):
    """从 VSI 日志提取崩溃原因 (运行时明确输出)"""
    text = ""
    if Path(log_path).is_file():
        text += Path(log_path).read_text(errors="ignore")
    olog = Path(out_dir) / "output-log.txt"
    if olog.is_file():
        text += "\n" + olog.read_text(errors="ignore")
    patterns = [
        (r"libdivsufsort|generalized text length", "salmon index 崩溃 (scaffold 过短)"),
        (r"spades", "SPAdes 组装失败"),
        (r"MemoryError", "内存不足"),
        (r"FileNotFoundError|No such file or directory", "文件缺失"),
        (r"salmon", "salmon 崩溃"),
        (r"CalledProcessError", "子命令返回非0"),
        (r"Traceback", "Python 异常"),
        (r"ERROR", "VSI 报告错误"),
    ]
    for pat, desc in patterns:
        if re.search(pat, text):
            return desc
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    return "未知错误: " + (" | ".join(l[:80] for l in lines[-3:]) if lines else "无日志")


def run_vsi(ref_fa, r1, r2, out_dir, threads, vsi_path, salmon_bin, checkv_db=None, genus_avg_len=0):
    """运行 Virseqimprover。
    返回 (scaffold_fa, status, reason):
      status = "ok" 延伸成功 / "no_extend" 正常结束未延伸 / "error" 崩溃 / "timeout" 超时
    """
    d = Path(out_dir); d.mkdir(parents=True, exist_ok=True)
    vsi = vsi_path if vsi_path and Path(vsi_path).is_file() else None
    if not vsi:
        return None, "error", "Virseqimprover.py 未找到"
    cmd = ["python", str(vsi), "-1", str(r1), "-2", str(r2),
           "-scaffold", str(ref_fa), "-o", str(d),
           "-salmon", str(salmon_bin), "-t", str(threads)]
    if checkv_db:
        cmd += ["-checkv_db", str(checkv_db)]
    if genus_avg_len > 0:
        cmd += ["-genus_avg_len", str(genus_avg_len)]
    ok_marker = d / "run.ok"
    if ok_marker.is_file():
        for name in ["scaffold.fasta", "pilon_out.fasta"]:
            f = d / name
            if f.is_file() and f.stat().st_size > 0:
                return str(f), "ok", "已完成"
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=14400, check=False)
    except subprocess.TimeoutExpired:
        print(f"    [VSI-TIMEOUT] {Path(ref_fa).stem[:50]} | 超过 4h", flush=True)
        return None, "timeout", "VSI 超过 4h 超时"
    vsi_log = d / "vsi_stdout.log"
    with open(vsi_log, 'w') as lf:
        lf.write(result.stdout or "")
        if result.stderr:
            lf.write("\n=== STDERR ===\n")
            lf.write(result.stderr)
    orig_len = _read_seq_len(ref_fa)
    if result.returncode == 0:
        for name in ["scaffold.fasta", "pilon_out.fasta"]:
            f = d / name
            if f.is_file() and f.stat().st_size > 0:
                final_len = _read_final_length(d)
                if final_len is not None and orig_len is not None and final_len > orig_len:
                    ok_marker.write_text("ok")
                    return str(f), "ok", f"延伸 {orig_len}→{final_len}bp"
                return str(f), "no_extend", f"未延伸 (final={final_len or '?'}bp, orig={orig_len or '?'}bp)"
        return None, "error", "VSI 正常结束但无 scaffold 输出"
    reason = _extract_vsi_error(vsi_log, d)
    return None, "error", reason


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
    合并两种来源并去重:
      - flye_sample_map 精确映射 (Flye contig 自身及同簇 Flye 成员)
      - vclust clusters 成员前缀 (正常 contig)
    NC_/ref| 前缀返回 None (参考序列, 无需 rescue)。"""
    seen = set()
    result = []

    def _add(s):
        if s not in seen:
            seen.add(s)
            result.append(s)
        return max_samples > 0 and len(result) >= max_samples

    # 1. Flye 精确映射 (contig 自身是 Flye 产物)
    if flye_sample_map and contig_id in flye_sample_map:
        for s in flye_sample_map[contig_id]:
            if _add(s):
                return result

    # 2. 查 clusters 补充成员样本
    if clusters:
        for cname, info in clusters.items():
            if contig_id in info["members"]:
                if max_samples > 0 and fasta_info and len(info["members"]) > max_samples:
                    ranked = sorted(info["members"],
                                    key=lambda m: fasta_info.get(m, {}).get("length", 0),
                                    reverse=True)
                else:
                    ranked = list(info["members"])
                for m in ranked:
                    if m == contig_id:
                        continue
                    if flye_sample_map and m in flye_sample_map:
                        # Flye 成员 → 精确映射
                        for s in flye_sample_map[m]:
                            if _add(s):
                                return result
                    elif m.startswith('flye_'):
                        # flye_ 前缀但无映射 → 跳过 (无 reads 可追溯)
                        continue
                    else:
                        # 正常成员 → 前缀
                        if _add(m.split('_')[0]):
                            return result
                break

    # 3. 回退: 前缀解析
    if not result:
        parts = contig_id.split('_')
        if len(parts) < 2 or parts[0].startswith('NC_') or parts[0].startswith('ref|'):
            return None
        return [parts[0]]

    return result


# ══════════════════════════════════════════════════════════════
# 分支 A: CheckV 并行评估
# ══════════════════════════════════════════════════════════════

def branch_a(ref_fasta, work_dir, checkv_db, threads, jobs, threshold=90.0, resume=True):
    """CheckV 并行评估 centroids。返回 (pass_fa, fail_fa, pass_count, fail_count)。支持断点续传。"""
    d = Path(work_dir) / "branch_a"; d.mkdir(parents=True, exist_ok=True)

    # 断点续传: 最终输出文件已存在 → 直接恢复 (resume=False 时清除重算)
    pass_fa = d / "branchA_pass.fasta"
    fail_fa = d / "branchA_fail.fasta"
    if resume and pass_fa.is_file() and fail_fa.is_file():
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

        # 断点续传: completeness.tsv 已存在 → 直接解析 (resume=False 时重算)
        if resume and cv_tsv.is_file():
            pids, fids, sids, _ = parse_checkv(cv_tsv, threshold)
        else:
            SeqIO.write(chunk, chunk_fa, "fasta")
            qs = run_checkv(chunk_fa, tmp_dir / "checkv_out", checkv_db, t_per)
            pids, fids, sids, _ = parse_checkv(qs, threshold)
        for r in chunk:
            if r.id in pids:
                lp.append(r)
            elif r.id in sids:
                # NA 但 ≥2000bp → 还有机会, 当 fail 处理进分支 B/C
                # NA 且 <2000bp → short fragment无蛋白基因, 跳过
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

def _cleanup_merged(r1, r2, merged_reads_dir=None):
    """删除单个 contig 的合并 reads (仅限 merged_reads 目录内), 释放磁盘。
    绝不删除 00b 原始 reads / 软链接 (单样本 contig 直接引用 00b 路径)。"""
    if merged_reads_dir is None:
        return
    mr = str(Path(merged_reads_dir).resolve())
    for f in (r1, r2):
        if f and os.path.isfile(f):
            try:
                if str(Path(f).resolve()).startswith(mr):
                    os.unlink(f)
            except OSError:
                pass


def branch_b(fail_fa, fastq_dir, work_dir, checkv_db, threads, jobs, vsi_path, salmon_bin, clusters=None, fasta_info=None, max_vsi_samples=10, min_vsi_len=2000, threshold=90.0, genus_map=None, genus_lens=None, flye_sample_map=None, seg_meta=None, resume=True):
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
        fail_fa_out = d / "branchB_fail.fasta"
        SeqIO.write(list(fail_records), fail_fa_out, "fasta")
        print(f"  分支 B: 全部 <{min_vsi_len}bp, 全部直升分支C")
        return None, str(fail_fa_out), 0, total

    # 所有长度达标的序列都进 VSI 任务池
    # 属分类信息在 _do 内部用于决定 genus_avg_len 参数 (有属传实际值, 无属传 0)
    fail_records = vsi_records
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

        ok_marker = out_dir / "run.ok"
        fail_marker = out_dir / ".vsi_failed"
        err_marker = out_dir / ".vsi_error"
        if ok_marker.is_file() and resume:
            if vsi_fa.is_file() and vsi_fa.stat().st_size > 0:
                with lock: stats["ok"] += 1
                return sid, str(vsi_fa), "已完成"
        elif ok_marker.is_file() and not resume:
            # [R6c] --no-resume: 已完成目录清除后重算 (err/fail 语义不变)
            import shutil as _shutil
            _shutil.rmtree(str(out_dir), ignore_errors=True)
        elif err_marker.is_file():
            # 延伸报错(崩溃/超时/缺reads): 最多重跑 1 次, 之后跳过进 branch C
            retries = 0
            try:
                retries = int(err_marker.read_text().strip().split("|")[0])
            except Exception:
                pass
            if retries <= 1:
                err_marker.unlink()
                print(f"    [RE-TRY] {sid[:50]} | 报错重跑 (第{retries}次)", flush=True)
            else:
                with lock: stats["skip"] += 1
                return sid, None, "报错(已达重跑上限): " + err_marker.read_text()
        elif fail_marker.is_file():
            # 延伸失败(正常结束未延伸, 如低覆盖度) → 跳过, 不重跑
            with lock: stats["skip"] += 1
            return sid, None, "未延伸(跳过): " + fail_marker.read_text()
        elif vsi_fa.is_file() and vsi_fa.stat().st_size > 0:
            # 有输出目录但无 run.ok → 不完整的残留, 删除重新跑
            import shutil as _shutil
            _shutil.rmtree(str(out_dir), ignore_errors=True)
            print(f"    [RE-RUN] {sid[:50]} | 清除残留, 重新运行", flush=True)

        sample_ids = _contig_cluster_samples(sid, clusters, fasta_info, max_vsi_samples, flye_sample_map=flye_sample_map) if clusters else [sample]
        if sample_ids is None:
            with lock: stats["skip"] += 1
            return sid, None, "参考序列(不在cluster中)"
        r1, r2, nsamp = _gather_cluster_reads(fastq_dir, sample_ids, merged_reads_dir, prefix=sid[:80])
        if not r1:
            out_dir.mkdir(parents=True, exist_ok=True)
            err_marker.write_text(f"1|reads 未找到 (samples: {sample_ids[:5]}{'...' if len(sample_ids)>5 else ''})")
            with lock: stats["fail"] += 1
            return sid, None, f"reads 未找到 (samples: {sample_ids[:5]}{'...' if len(sample_ids)>5 else ''})"

        sf = scaffold_dir / f"{sid}.fasta"
        SeqIO.write([rec], sf, "fasta")

        # 查找属平均长度 (CheckV NA 时的备选截止条件)
        gal = 0
        if genus_map and genus_lens:
            g = genus_map.get(sid) or genus_map.get(sample, "")
            g_clean = g.replace("g__", "").replace("G__", "")
            # 多级回退: 原名 → 带 g__ 前缀 → 去 ICTV 前缀 (Acarallexivirus → Allexivirus)
            gal = genus_lens.get(g_clean, 0) or genus_lens.get("g__" + g_clean, 0)
            if gal == 0:
                # ICTV 亚属名通常有前缀(Acar/Alpha/Beta...), 尝试去掉前缀匹配 NCBI 属名
                for prefix in ("Acar", "Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta",
                               "Eta", "Theta", "Iota", "Kappa", "Lambda", "Mu", "Nu", "Xi",
                               "Omicron", "Pi", "Rho", "Sigma", "Tau", "Upsilon", "Phi"):
                    if g_clean.startswith(prefix) and len(g_clean) > len(prefix) + 2:
                        c = g_clean[len(prefix):]
                        for c2 in (c, c.capitalize()):
                            gal = genus_lens.get(c2, 0) or genus_lens.get("g__" + c2, 0)
                            if gal > 0: break
                        if gal > 0: break

        task_log = log_dir / f"{sid}.log"
        with open(task_log, "w") as lf:
            lf.write(f"# VSI multi-sample task: {sid}\n"
                     f"# Cluster samples ({nsamp}): {', '.join(sample_ids)}\n"
                     f"# Merged Reads: {r1}, {r2}\n# Scaffold: {sf}\n"
                     f"# Genus: {g if genus_map else 'N/A'} → avg_len: {gal}bp\n")

        _, vsi_status, vsi_reason = run_vsi(sf, r1, r2, out_dir, threads_per_job, vsi_path, salmon_bin, checkv_db, genus_avg_len=gal)

        vsi_fa = out_dir / "scaffold-truncated" / "scaffold.fasta"
        if not vsi_fa.is_file():
            vsi_fa = out_dir / "scaffold.fasta"
        if not vsi_fa.is_file():
            vsi_fa = out_dir / "pilon_out.fasta"
        if vsi_status == "ok" and vsi_fa.is_file():
            out_dir.mkdir(parents=True, exist_ok=True)
            ok_marker.touch()
            with lock: stats["ok"] += 1
            print(f"    [C-OK] {sid[:60]}  {vsi_reason}  ✓", flush=True)
            _cleanup_merged(r1, r2, merged_reads_dir)
            return sid, str(vsi_fa), vsi_reason
        out_dir.mkdir(parents=True, exist_ok=True)
        if vsi_status == "no_extend":
            # 延伸失败(正常结束未延伸, 如低覆盖度) → 跳过, 不重跑
            fail_marker.write_text(vsi_reason)
            with lock: stats["skip"] += 1
            print(f"    [NO-EXTEND] {sid[:60]}  {vsi_reason}", flush=True)
            _cleanup_merged(r1, r2, merged_reads_dir)
            return sid, None, vsi_reason
        # 报错(崩溃/超时): 标记 .vsi_error (resume 时可重跑 1 次)
        retries = 0
        try:
            retries = int(err_marker.read_text().strip().split("|")[0])
        except Exception:
            pass
        err_marker.write_text(f"{retries + 1}|{vsi_reason}")
        with lock: stats["fail"] += 1
        print(f"    [VSI-{vsi_status.upper()}] {sid[:60]}  {vsi_reason}", flush=True)
        _cleanup_merged(r1, r2, merged_reads_dir)
        return sid, None, f"[{vsi_status}] {vsi_reason}"

    print(f"  并行启动 {jobs} 个 Virseqimprover (各 {threads_per_job} 线程)...")
    results = {}
    if jobs > 1 and total > 1:
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futures = {ex.submit(_do, rec): rec.id for rec in fail_records}
            for f in _safe_tqdm(as_completed(futures), total=total, desc="  分支 B", unit="task"):
                sid, fa, msg = f.result()
                if fa:
                    results[sid] = fa
    else:
        for rec in _safe_tqdm(fail_records, desc="  分支 B", unit="task"):
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
        # VSI 全部失败 → 所有原始序列写进 fail_fa 存档 (2026-10-11 起无分支 C, 记为未拯救)
        fail_fa_out = d / "branchB_fail.fasta"
        SeqIO.write(list(fail_records) + list(short_records), fail_fa_out, "fasta")
        return None, str(fail_fa_out), 0, len(fail_records) + len(short_records)

    extended_fa = d / "branchB_extended.fasta"
    SeqIO.write(vsi_records, extended_fa, "fasta")

    qs = run_checkv(extended_fa, d / "checkv_out", checkv_db, threads)
    pass_ids, fail_ids, skip_ids, cv_meta = parse_checkv(qs, threshold)

    # ── CheckV low-confidence 智能纠偏: 期望长度与属平均差 >20% → 用 genus ──
    if fail_ids and genus_lens:
        for cid in list(fail_ids):
            info = cv_meta.get(cid, {})
            if info.get('confidence', '') != 'low':
                continue
            cv_elen = info.get('expected_len', 0)
            if cv_elen <= 0:
                continue
            gal = 0
            for key in (cid, cid.split('_')[0]):
                g = genus_map.get(key, '')
                g_clean = g.replace('g__', '').replace('G__', '')
                gal = genus_lens.get(g_clean, 0) or genus_lens.get("g__" + g_clean, 0)
                if gal == 0:
                    for prefix in ("Acar", "Alpha", "Beta", "Gamma", "Delta", "Epsilon",
                                   "Zeta", "Eta", "Theta", "Iota", "Kappa", "Lambda", "Mu"):
                        if g_clean.startswith(prefix) and len(g_clean) > len(prefix) + 2:
                            c = g_clean[len(prefix):]
                            for c2 in (c, c.capitalize()):
                                gal = genus_lens.get(c2, 0) or genus_lens.get("g__" + c2, 0)
                                if gal > 0: break
                            if gal > 0: break
                if gal > 0: break
            if gal > 0:
                diff = abs(cv_elen - gal) / gal
                if diff > 0.20:
                    fail_ids.discard(cid)
                    skip_ids.add(cid)

    # ── CheckV=NA 或 <threshold 的 genus 回退拯救 ──
    # VSI 自身已停止延伸, CheckV 数值不可靠时, 用属平均长度兜底
    genus_rescued = set()
    genus_tol = 0.85  # 属长度容忍度 (原分支D已归档, 现仅分支 B 使用)
    rescue_candidates = set()
    if skip_ids: rescue_candidates.update(skip_ids)
    if fail_ids: rescue_candidates.update(fail_ids)
    if rescue_candidates and genus_map and genus_lens:
        for rec in vsi_records:
            if rec.id in rescue_candidates:
                gal = 0
                sinfo = None
                for key in (rec.id, rec.id.split('_')[0]):
                    g = genus_map.get(key, '')
                    g_clean = g.replace('g__', '').replace('G__', '')
                    gal = genus_lens.get(g_clean, genus_lens.get(g, 0))
                    if gal > 0:
                        sinfo = seg_meta.get(g_clean) if seg_meta else None
                        break
                if _genus_len_ok(len(rec.seq), gal, sinfo, genus_tol):
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
    vsi_fail_n = 0
    for rec in fail_records:
        if rec.id in results:
            # VSI 成功 → CheckV 决定 pass/fail
            in_fail = rec.id in fail_ids and rec.id not in genus_rescued
            in_skip = rec.id in skip_ids and rec.id not in genus_rescued
            if in_fail or in_skip:
                vsi_match = [r for r in vsi_records if r.id == rec.id]
                fail_to_write.append(vsi_match[0] if vsi_match else rec)
        else:
            # VSI 失败 (超时/崩溃/未延伸/缺reads) → 原始序列记为未拯救 (仅存档)
            vsi_fail_n += 1
            fail_to_write.append(rec)
    # short_records (< min_vsi_len) 也加入 fail 存档 (不参与 VSI, 记为未拯救)
    if short_records:
        fail_to_write.extend(short_records)
    SeqIO.write(fail_to_write, fail_fa_out, "fasta")

    n_fail_total = len(fail_ids - genus_rescued) + len(skip_ids - genus_rescued) + len(short_records) + vsi_fail_n
    na_fail = len(skip_ids - genus_rescued)
    num_fail = len(fail_ids - genus_rescued)
    print(f"  分支 B: pass={len(pass_ids):,}  未拯救={n_fail_total} (含{len(short_records)}短序列, {na_fail}条NA, {num_fail}条<阈值)")
    return str(pass_fa) if Path(pass_fa).is_file() else None, \
           str(fail_fa_out) if Path(fail_fa_out).is_file() else None, \
           len(pass_ids), n_fail_total


# ══════════════════════════════════════════════════════════════
# 属分类 / 属长工具 (分支 B 回退判定 + 报告共用)
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
    # 2. 补充 ICTV→NCBI 映射: ICTV 亚属名 → NCBI 属名 → genus_lens 中的平均长度
    ictv_map_path = Path(genus_len_path).parent / "ictv_to_ncbi_genus.json"
    if ictv_map_path.is_file():
        try:
            ictv_ncbi = json.load(open(ictv_map_path))
            added = 0
            for ictv_name, ncbi_name in ictv_ncbi.items():
                if ictv_name not in mapping and ncbi_name in mapping:
                    mapping[ictv_name] = mapping[ncbi_name]
                    added += 1
            if added:
                print(f"  ICTV→NCBI 映射补充: {added} 个 genus")
        except Exception:
            pass

    print(f"  genus_len 加载: {len(mapping)} 个 genus")
    return mapping


def _load_genus_seg_meta(genus_len_path):
    """加载 genus_len v2 表的节段属元数据 (第3列起: is_segmented/ref_min/ref_max).

    旧版 2 列表 (genus/total) 无这些列, 返回空 dict, 调用方自动回退到属均值规则。
    返回 {genus(原名与去g__前缀): {"is_segmented": bool, "ref_min": float, "ref_max": float}}
    """
    meta = {}
    if not genus_len_path or not Path(genus_len_path).is_file():
        return meta
    with open(genus_len_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if parts[0].strip().lower() == 'genus':   # 表头
                if len(parts) < 5:
                    return meta                        # 老表: 无节段列
                continue
            if len(parts) < 5:
                continue
            genus = parts[0].strip()
            try:
                is_seg = parts[2].strip().lower() in ('yes', 'true', '1')
                rmin = float(parts[3])
                rmax = float(parts[4])
            except (ValueError, IndexError):
                continue
            if rmin <= 0 or rmax < rmin:
                continue
            for gname in (genus, genus.replace('g__', '').replace('G__', '')):
                meta[gname] = {"is_segmented": is_seg, "ref_min": rmin, "ref_max": rmax}
    if meta:
        print(f"  genus_len 节段元数据: {sum(1 for v in meta.values() if v['is_segmented'])} 节段属 / {len(meta)}")
    return meta


def _genus_len_ok(contig_len, gal, seg_info, genus_tol=0.85):
    """分支 B/D 长度判定 (节段属感知).

    节段属且有范围信息: ref_min*tol <= 长度 <= ref_max*(2-tol)
      -- tol=0.85 即下界允许短 15%, 上界允许长 15%; 范围覆盖该属全部完整单元
         (单段到全基因组), 节段病毒的完整单段 contig 可被拯救。
    非节段属 / 未知 / 无范围信息: 沿用旧规则 长度 >= 属平均长度*tol (单侧)。
    """
    if seg_info and seg_info.get("is_segmented") and seg_info.get("ref_min") and seg_info.get("ref_max"):
        lo = seg_info["ref_min"] * genus_tol
        hi = seg_info["ref_max"] * (2.0 - genus_tol)
        return lo <= contig_len <= hi
    if gal and gal > 0:
        return contig_len >= gal * genus_tol
    return False


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
                         taxonomy_tsv=None, genus_lens_path=None, min_vsi_len=2000,
                         prepass=0):
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

            # v6.6.1: 分支 C 的 Completo 直达 pass (BLASTN qcov≥98%) 此前没有对应分支,
            # 它的 id 就是 contig id 且不进属组, 全被记成 fail, 273/9212 条 C 拯救在报告里凭空消失
            elif cid in pass_c and cid not in contig_to_genus:
                info = pass_c[cid]
                flen = info['length']
                pct2 = f"{flen/gal*100:.1f}%" if gal > 0 else "-"
                desc = info.get('desc', '')
                cv = "-"
                for part in desc.split():
                    if part.startswith('CheckV='):
                        cv = part.replace('CheckV=', '')
                        break
                tf.write(f"{cid}\t{clen}\t{flen}\tC\tblastn_completo\t{genus}\t{gal:.0f}\t{pct2}\t{cv}\tBLASTN qcov≥98% 直达\n")

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
                    "should_rescue_but_skipped": "Supported length but not rescued",
                    "below_all_thresholds": "未通过任何分支",
                    "unknown": "未拯救",
                }.get(reason, reason)
                tf.write(f"{cid}\t{clen}\t-\tfail\t-\t{genus}\t{gal:.0f}\t{pct}\t-\t{note}\n")

    # v6.6.1: cnt_a 含 orchestrator 预通过 (prepass) 的免拯救 contig, 它们不在候选表内,
    # 不能进候选表的「已拯救/未拯救」分母, 否则未拯救数与原因分布会自相矛盾 (原因数 > 未拯救数)
    total = len(centroids_records)
    rescued_rows = (cnt_a - max(0, prepass)) + cnt_b + cnt_c + cnt_d
    rescued = cnt_a + cnt_b + cnt_c + cnt_d  # 含 prepass, 仅用于总量展示
    failed = total - rescued_rows

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
        a_inrun = cnt_a - max(0, prepass)  # 与 TSV 行口径一致: 只数候选表内的 A 分支拯救
        a_note = f"（另含 {prepass} 条免拯救 prepass, 不在候选表内）" if prepass else ""
        mf.write(f"| A | CheckV ≥90% (蛋白完整) | {a_inrun} {a_note} | {a_inrun/total*100:.1f}% |\n")
        mf.write(f"| B | VSI reads延伸 + genus回退 | {cnt_b} | {cnt_b/total*100:.1f}% |\n")
        if cnt_c or cnt_d:   # 历史输出兼容 (2026-10-11 起 C/D 已归档, 新运行恒为 0 不再列出)
            mf.write(f"| C | BLASTN + RGA/ragtag | {cnt_c} | {cnt_c/total*100:.1f}% |\n")
            mf.write(f"| D | genus_len 属水平兜底 | {cnt_d} | {cnt_d/total*100:.1f}% |\n")
        mf.write(f"| **合计 (候选表内)** | | **{rescued_rows}** | **{rescued_rows/total*100:.1f}%** |\n")
        mf.write(f"| 未拯救 (候选表内) | | {failed} | {failed/total*100:.1f}% |\n")
        if prepass:
            mf.write(f"| 免拯救 (prepass, 候选表外) | CheckV ≥90% | {prepass} | - |\n")
        mf.write("\n")

        if fail_reasons:
            mf.write(f"## 未拯救原因分布\n\n")
            mf.write(f"| 原因 | 数量 | 占未拯救 |\n")
            mf.write(f"|------|------|---------|\n")
            _reasons_sum = sum(fail_reasons.values())
            for reason, count in sorted(fail_reasons.items(), key=lambda x: -x[1]):
                label = reason_labels.get(reason, reason)
                _pct = (count / failed * 100) if failed else 0.0
                mf.write(f"| {label} | {count} | {_pct:.1f}% |\n")
            mf.write(f"| **合计** | **{_reasons_sum}** | **{(_reasons_sum / failed * 100) if failed else 0.0:.1f}%** |\n")
            mf.write("\n")

        mf.write(f"最终无冗余 vOTU: **{n_final}**\n\n")
        mf.write(f"详细: [{tsv_path.name}]({tsv_path.name})\n")

    print(f"  报告: {tsv_path}")
    print(f"  摘要: {md_path}")


# ══════════════════════════════════════════════════════════════
# 主流程
# ══════════════════════════════════════════════════════════════

def main():
    p = argparse.ArgumentParser(description="病毒基因组两支路级联拯救 v3.2 (A: CheckV≥90% 直通 → B: VSI reads 延伸; C/D 已归档见 archive/20261011_rescue_branch_cd/)")
    p.add_argument("--centroids", "-c", required=True, help="输入 centroids FASTA")
    p.add_argument("--clusters-tsv", required=True, help="vclust 聚类结果 TSV")
    p.add_argument("--split-dir", required=True, help="per-cluster 拆分文件目录 (split_fastas/)")
    p.add_argument("--output-dir", "-o", required=True)
    p.add_argument("--fastq-dir", "-fq", required=True, help="原始 reads 目录")
    p.add_argument("--checkv-db", "-cv", required=True)
    p.add_argument("--blast-db", "-db", default=None, help="[已废弃 2026-10-11] 分支C已归档, 接受但忽略")
    p.add_argument("--virseqimprover-path", default=None)
    p.add_argument("--salmon-bin", default=os.path.expanduser("~/mambaforge/envs/Virseqimprover/bin/salmon"))
    p.add_argument("--max-vsi-samples", type=int, default=10, help="VSI 最大合并样本数 (0=不限制, 默认10)")
    p.add_argument("--min-vsi-len", type=int, default=2000, help="VSI 最小 contig 长度 bp (短于此值跳过VSI记为未拯救, 默认2000)")
    p.add_argument("--checkv-threshold", type=float, default=90.0, help="CheckV completeness 通过阈值 (默认90, Plant建议80)")
    p.add_argument("--taxonomy-tsv", default=None, help="05_Taxonomy 分类结果 final_integrated_classification.tsv (B 属长回退判定用)")
    p.add_argument("--cdhit-fa", default=None, help="[已废弃 2026-10-11] 分支C已归档, 接受但忽略")
    p.add_argument("--genus-len", default=None, help="genus_len 参考文件 (属平均长度, B 属长回退判定用)")
    p.add_argument("--genus-tolerance", type=float, default=0.85, help="[已废弃 2026-10-11] 分支D已归档, 接受但忽略 (B 内部固定 0.85)")
    p.add_argument("--threads", "-t", type=int, default=64)
    p.add_argument("--jobs", "-j", type=int, default=4, help="Virseqimprover 并行数")
    p.add_argument("--ani", type=float, default=0.95, help="最终 vclust ANI")
    p.add_argument("--qcov", type=float, default=0.85, help="最终 vclust QCOV")
    p.add_argument("--resume", action="store_true", default=True, help="断点续传 (默认开启, --no-resume 关闭)")
    p.add_argument("--no-resume", action="store_false", dest="resume", help="禁用断点续传, 强制重跑")
    p.add_argument("--flye-sample-map", default=None,
                   help="Flye 共组装 contig→样本 映射 JSON (merge 阶段产出, 供 VSI reads 回退)")
    args = p.parse_args()
    for _dep in ("blast_db", "cdhit_fa", "genus_tolerance"):
        if getattr(args, _dep, None) not in (None, 0.85):
            print(f"  [DEPRECATED] --{_dep.replace('_', '-')} 已随分支C/D归档被忽略 (2026-10-11)")

    # 进程守护
    try:
        import process_guard
        process_guard.install()
    except Exception:
        pass

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

    # 统计长短序列分布 (短序列跳过 B, 记为未拯救)
    n_long = sum(1 for r in centroids_records if len(r.seq) >= args.min_vsi_len)
    n_short = len(centroids_records) - n_long
    if n_short > 0:
        print(f"  长度分布: ≥{args.min_vsi_len}bp={n_long} 条, <{args.min_vsi_len}bp={n_short} 条 (A后跳过B, 不拯救)")

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
        fa_a_pass, fa_a_fail, cnt_a, cnt_a_fail = branch_a(str(tmp_centroids), out, args.checkv_db, args.threads, args.jobs, args.checkv_threshold, resume=args.resume)

    # 合并 orchestrator 预通过的 CheckV pass 计数 (免拯救, 不在 rescue centroids 中)
    # v6.6.1: n_pre 单独留存, 报告里必须与候选表内的拯救数分开记账
    n_pre = 0
    pre_pass_file = out / "branch_a_prepass.txt"
    if pre_pass_file.is_file():
        n_pre = sum(1 for _ in open(pre_pass_file))
        cnt_a += n_pre
        print(f"  分支 A: +{n_pre} 条预通过 (CheckV≥90%, 免拯救)")

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
        seg_meta = {}
        if getattr(args, 'taxonomy_tsv', None):
            genus_map = _load_taxonomy(args.taxonomy_tsv)
        if getattr(args, 'genus_len', None):
            genus_lens = _load_genus_len(args.genus_len)
            seg_meta = _load_genus_seg_meta(args.genus_len)
        fa_b_pass, fa_b_fail, cnt_b, cnt_b_fail = branch_b(fa_a_fail, args.fastq_dir, out, args.checkv_db, args.threads, args.jobs, args.virseqimprover_path, args.salmon_bin, clusters, fasta_info, args.max_vsi_samples, args.min_vsi_len, args.checkv_threshold, genus_map=genus_map, genus_lens=genus_lens, flye_sample_map=flye_sample_map, seg_meta=seg_meta, resume=args.resume)

    # 3. 合并
    print("\n── Step 3: 合并 ──")
    d4 = out / "merged"; d4.mkdir(exist_ok=True)
    merged = d4 / "all_HQ.fasta"
    centroids_final = out / "centroids" / "final_centroids.fasta"

    if args.resume and centroids_final.is_file() and centroids_final.stat().st_size > 0:
        print("  [RESUME] Step 3+4 已有最终结果, 跳过")
        total_m = sum(1 for _ in SeqIO.parse(merged, "fasta")) if merged.is_file() else cnt_a + cnt_b
        n_final = sum(1 for _ in SeqIO.parse(centroids_final, "fasta"))
    else:
        with open(merged, "w") as mf:
            for fp in [fa_a_pass, fa_b_pass]:
                if fp and os.path.isfile(fp):
                    with open(fp) as inf:
                        mf.write(inf.read())
        total_m = sum(1 for _ in SeqIO.parse(merged, "fasta"))
        print(f"  A:{cnt_a}  B:{cnt_b}  →  {total_m} 条  →  {merged}")

        if total_m < 2:
            print("  [SKIP] <2 条, 跳过最终 vclust")
            centroids_final.parent.mkdir(parents=True, exist_ok=True)
            with open(centroids_final, "w") as cf, open(merged) as mf2:
                cf.write(mf2.read())
            n_final = total_m
        else:
            # 4. 最终 vclust 去重
            print("\n── Step 4: vclust 最终去重 ──")
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
    print("\n── Step 5: 生成报告 ──")
    _write_rescue_report(out, centroids_records, clusters,
                          fa_a_pass, fa_b_pass, None, None,
                          cnt_a, cnt_b, 0, 0, n_final,
                          taxonomy_tsv=getattr(args, 'taxonomy_tsv', None),
                          genus_lens_path=getattr(args, 'genus_len', None),
                          min_vsi_len=args.min_vsi_len,
                          prepass=n_pre)

    elapsed = (datetime.now() - start).total_seconds()
    print(f"\n{'=' * 60}")
    print(f"  分支 A (CheckV):      {cnt_a:,}")
    print(f"  分支 B (VSI):         {cnt_b:,}")
    print(f"  最终无冗余:           {n_final:,}")
    print(f"  最终输出:             {centroids_final}")
    print(f"  耗时:                 {elapsed / 60:.1f} min")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
