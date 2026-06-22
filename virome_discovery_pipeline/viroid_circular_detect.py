#!/usr/bin/env python3
"""
Viroid 环状基因组检测脚本 (独立运行)
======================================

参照 Virseqimprover checkCircularity 逻辑:
  对 viroid 候选 (200-1000bp) 进行自 BLASTN 环状检测。
  提取 5' 端为 subject, 3' 端为 query, 自比对判定是否环状。

用法:
  python viroid_circular_detect.py -i viroids.candidate.fasta -o circular_output -t 64
"""

import argparse, subprocess, sys, os
from pathlib import Path
from datetime import datetime
from Bio import SeqIO


def run(cmd, desc="", check=True):
    print(f"  [RUN] {desc}", flush=True)
    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=check)
    except subprocess.CalledProcessError as e:
        print(f"  [FAIL] 退出码 {e.returncode}")
        if e.stderr:
            print(f"        {e.stderr.strip()[:500]}")
        if check:
            sys.exit(1)
        return None


def detect_circular_genomes(fasta_in, out_dir, threads=4, min_id=95, min_qcov=90):
    """自 BLASTN 环状检测。返回 (circular_ids: set, log_dir: Path)"""
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    log_dir = out / "logs"; log_dir.mkdir(exist_ok=True)
    tmp_dir = out / "tmp"; tmp_dir.mkdir(exist_ok=True)

    circular = set()
    seqs = list(SeqIO.parse(fasta_in, "fasta"))
    total = len(seqs)
    print(f"  输入: {total} 条 viroid 候选 (200-1000bp)")

    for i, rec in enumerate(seqs):
        slen = len(rec.seq)
        check_len = max(40, slen // 3)
        if slen < check_len * 2:
            continue

        sid = rec.id.replace('/', '_').replace('|', '_')[:60]
        prefix = tmp_dir / sid

        # 5' 端 → subject, 3' 端 → query
        subj_fa = tmp_dir / f"{sid}_subject.fa"
        qry_fa = tmp_dir / f"{sid}_query.fa"
        bt_out = tmp_dir / f"{sid}_blastn.txt"

        with open(subj_fa, "w") as sf:
            sf.write(f">{sid}_subject\n{str(rec.seq[:slen - check_len])}\n")
        with open(qry_fa, "w") as qf:
            qf.write(f">{sid}_query\n{str(rec.seq[slen - check_len:])}\n")

        task_log = log_dir / f"{sid}.log"
        with open(task_log, "w") as lf:
            lf.write(f"# Viroid circularity check: {sid} ({slen}bp)\n")
            lf.write(f"# check_len={check_len}bp\n")

        run(["makeblastdb", "-in", str(subj_fa), "-dbtype", "nucl"], f"makeblastdb: {sid}", check=False)
        run(["blastn", "-query", str(qry_fa), "-db", str(subj_fa),
             "-num_threads", str(threads), "-outfmt", "7", "-out", str(bt_out)],
            f"blastn self: {sid}", check=False)

        if bt_out.is_file() and bt_out.stat().st_size > 0:
            with open(bt_out) as bf:
                for line in bf:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 12:
                        continue
                    try:
                        ident = float(parts[2])
                        aln_len = int(parts[3])
                        if ident >= min_id and aln_len >= int(check_len * min_qcov / 100):
                            circular.add(rec.id)
                            print(f"    [CIRCULAR] {rec.id} ({slen}bp, id={ident:.1f}%, aln={aln_len}bp)")
                            with open(task_log, "a") as lf:
                                lf.write(f"# RESULT: CIRCULAR (ident={ident:.1f}%, aln={aln_len}bp)\n")
                            break
                    except (ValueError, IndexError):
                        continue

        if rec.id not in circular and (i + 1) % 50 == 0:
            print(f"  [{i+1}/{total}] 已检测 {len(circular)} 个环状...", flush=True)

    # 清理临时文件
    for f in tmp_dir.glob("*"):
        try:
            f.unlink()
        except OSError:
            pass
    try:
        tmp_dir.rmdir()
    except OSError:
        pass

    return circular, log_dir


def main():
    p = argparse.ArgumentParser(description="Viroid 环状基因组检测")
    p.add_argument("--input", "-i", required=True, help="输入 FASTA (viroid 候选, 200-1000bp)")
    p.add_argument("--output-dir", "-o", required=True, help="输出目录")
    p.add_argument("--threads", "-t", type=int, default=64)
    p.add_argument("--min-id", type=float, default=95, help="最小 identity %% (默认: 95)")
    p.add_argument("--min-qcov", type=float, default=90, help="最小覆盖 %% (默认: 90)")
    args = p.parse_args()

    start = datetime.now()
    print(f"开始: {start.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"输入: {args.input}")
    print(f"参数: min_id={args.min_id}% min_qcov={args.min_qcov}% threads={args.threads}")

    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
    circular_ids, log_dir = detect_circular_genomes(
        args.input, out, args.threads, args.min_id, args.min_qcov
    )

    # 写出环状 viroid
    circular_fa = out / "viroids.circular.fasta"
    non_circular_fa = out / "viroids.non_circular.fasta"
    circ_recs, non_recs = [], []
    for rec in SeqIO.parse(args.input, "fasta"):
        (circ_recs if rec.id in circular_ids else non_recs).append(rec)

    SeqIO.write(circ_recs, circular_fa, "fasta")
    SeqIO.write(non_recs, non_circular_fa, "fasta")

    elapsed = (datetime.now() - start).total_seconds()
    print(f"\n{'=' * 60}")
    print(f"  总输入:       {circ_recs + non_recs:,} 条")
    print(f"  环状 (HQ):    {len(circ_recs):,} → {circular_fa}")
    print(f"  非环状:       {len(non_recs):,} → {non_circular_fa}")
    print(f"  日志:         {log_dir}/")
    print(f"  耗时:         {elapsed:.1f}s")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
