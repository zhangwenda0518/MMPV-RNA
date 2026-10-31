#!/usr/bin/env python3
"""
batch_extract_reads.py — 靶向病毒 reads 提取 (known-virus 管线 Stage 4 extract 独立脚本)
=========================================================================================
[阶段重排] 自 batch_virus_variants.py 独立: 只做"从比对 BAM 提取每病毒×样本的 reads",
不做比对/summary 重算/consensus/变异检测。

输入:
  --summary      detect/filter 的 summary TSV (含 Sample/Virus 列; 编排器传 get_summary())
  --bam          variants 阶段产出的 BAM 目录 ({output_dir}/virus-bam/); 缺 .bai 自动补索引
  --output_dir   03_Variants 根目录 (产物写 {output_dir}/virus_reads/, 日志写 {output_dir}/logs/)

产物:
  virus_reads/{sample}.{virus}_R1/_R2/_single.fastq.gz   (full 阶段 de novo 组装主料)

用法 (编排器 Stage 4):
  python batch_extract_reads.py --summary <tsv> --bam <variants_dir>/virus-bam \
      --output_dir <variants_dir> --threads 15 --jobs 15 --resume
"""
import argparse
import os
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

# 提取 worker 单一来源 (batch_virus_variants 为 guarded import, 无副作用)
from batch_virus_variants import worker_extract_reads, safe_name

import pandas as pd


def build_bam_map(bam_dir, samples, threads, log_file):
    """glob *.bam → {sample: {"bam": path, "is_fasta": False}}, 与 variants --bam 分支同语义:
    stem 剥 .sorted/.sort/.bam, 同名优先取带 sorted 的, 缺/旧 .bai 自动补索引。"""
    import subprocess
    bam_map = {}
    for f in sorted(Path(bam_dir).glob("*.bam")):
        stem = f.name.replace(".sorted.bam", "").replace(".sort.bam", "").replace(".bam", "")
        if stem in samples and (stem not in bam_map or "sorted" in f.name):
            bam_map[stem] = {"bam": str(f), "is_fasta": False}
            bai = Path(str(f) + ".bai")
            if not bai.exists() or bai.stat().st_mtime < f.stat().st_mtime:
                subprocess.run(f"samtools index -@ {threads} '{f}'",
                               shell=True, executable='/bin/bash')
    return bam_map


def main():
    p = argparse.ArgumentParser(
        description="靶向病毒 reads 提取 (Stage 4 extract; 从比对 BAM 按病毒×样本抽取 R1/R2/single)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--summary", required=True,
                   help="detect/filter 的 summary TSV (需含 Sample/Virus 列)")
    p.add_argument("--bam", required=True,
                   help="BAM 目录 (variants 阶段产物 virus-bam/)")
    p.add_argument("--output_dir", required=True,
                   help="输出根目录 (03_Variants; 产物写 virus_reads/, 日志写 logs/)")
    p.add_argument("--threads", type=int, default=4, help="单任务 IO 线程数")
    p.add_argument("--jobs", type=int, default=4, help="并行任务数")
    p.add_argument("--resume", action="store_true",
                   help="跳过 R1 已存在 (>50B) 的任务")
    args = p.parse_args()

    bam_dir = Path(args.bam)
    if not bam_dir.is_dir():
        sys.exit(f"ERROR: BAM 目录不存在: {bam_dir} (先运行 --stage variants 产出比对 BAM)")

    summary = Path(args.summary)
    if not summary.is_file():
        sys.exit(f"ERROR: summary 不存在: {summary} (先运行 --stage detect/filter)")
    df = pd.read_csv(summary, sep="\t")
    for col in ("Sample", "Virus"):
        if col not in df.columns:
            sys.exit(f"ERROR: summary 缺少 {col} 列: {summary}")

    samples = set(df["Sample"].dropna().astype(str).unique())
    bam_map = build_bam_map(bam_dir, samples, args.threads,
                            str(Path(args.output_dir) / "logs" / "samtools_index.log"))
    missing = sorted(samples - set(bam_map))
    if not bam_map:
        sys.exit("ERROR: BAM 目录中没有任何匹配样本的 BAM — 先运行 --stage variants")
    if missing:
        print(f"[WARN] {len(missing)} 个样本无对应 BAM, 跳过 (示例: {', '.join(missing[:5])})",
              file=sys.stderr)

    out_root = Path(args.output_dir)
    d_reads = out_root / "virus_reads"
    d_logs = out_root / "logs"
    d_reads.mkdir(parents=True, exist_ok=True)
    d_logs.mkdir(parents=True, exist_ok=True)

    tasks = []
    for _, r in df[["Sample", "Virus"]].drop_duplicates().iterrows():
        s, v = str(r["Sample"]), str(r["Virus"])
        if s not in bam_map:
            continue
        prefix = f"{safe_name(s)}.{safe_name(v)}"
        ext = 'fasta.gz' if bam_map[s]["is_fasta"] else 'fastq.gz'
        tasks.append((
            s, v, bam_map[s]["bam"], bam_map[s]["is_fasta"],
            str(d_reads / f"{prefix}_R1.{ext}"),
            str(d_reads / f"{prefix}_R2.{ext}"),
            str(d_reads / f"{prefix}_single.{ext}"),
            args.threads, args.resume,
            str(d_logs / f"{prefix}_extract_reads.log"),
            str(d_logs / "extract_reads_master.log"),
        ))

    if not tasks:
        sys.exit("ERROR: 0 个提取任务 (summary 病毒×样本与 BAM 均不匹配)")

    n_done = 0
    print(f"🧾 提取任务: {len(tasks)} 个 (病毒×样本) | BAM: {len(bam_map)} 样本 | jobs={args.jobs}",
          file=sys.stderr, flush=True)
    with ProcessPoolExecutor(max_workers=min(len(tasks), max(1, args.jobs))) as ex:
        futures = [ex.submit(worker_extract_reads, t) for t in tasks]
        if tqdm:
            for _ in tqdm(as_completed(futures), total=len(futures),
                          desc="extract_reads", unit="task", dynamic_ncols=True):
                n_done += 1
        else:
            for _ in as_completed(futures):
                n_done += 1
    print(f"✅ 靶向 reads 提取完成: {n_done}/{len(tasks)} 任务 -> {d_reads}",
          file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
