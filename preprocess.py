#!/usr/bin/env python3
"""
preprocess.py — 数据预处理独立脚本 (清洗 + 去宿主)

  1. clean-data.py: Fastp → Seqkit → Clumpify
  2. host_depletion.py: Kraken2 → Align → Ribodetector

用法:
  python preprocess.py \
    -i /data/raw_fastq/ \
    -o /data/clean_out/ \
    --kraken2_db /db/kraken2_host/ \
    --host_align_db /db/bowtie2_host/host \
    --threads 20 --jobs 4 --rrna

依赖: clean-data.py, host_depletion.py (同目录)
"""

import argparse, subprocess, sys, os, logging, shutil
from pathlib import Path


def setup_logger(out_dir):
    logger = logging.getLogger("Preprocess")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    os.makedirs(out_dir, exist_ok=True)
    ch = logging.StreamHandler(); ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(ch)
    fh = logging.FileHandler(os.path.join(out_dir, 'preprocess.log'))
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s %(message)s'))
    logger.addHandler(fh)
    return logger


def run(cmd, logger, name):
    logger.info("[%s] %s", name, cmd)
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=None, stderr=subprocess.PIPE)
        logger.info("[%s] ✓", name)
        return True
    except subprocess.CalledProcessError as e:
        logger.error("[%s] ✗ (exit=%d)", name, e.returncode)
        if e.stderr:
            logger.error("  %s", e.stderr.decode()[-500:] if isinstance(e.stderr, bytes) else e.stderr[-500:])
        return False


def main():
    p = argparse.ArgumentParser(description="数据预处理 (清洗 + 去宿主)")
    p.add_argument("-i", "--input", required=True, help="原始 FASTQ 目录")
    p.add_argument("-o", "--output-dir", required=True, help="输出根目录 (生成 clean/ 和 host/ 子目录)")

    # 清洗
    g = p.add_argument_group("清洗 (clean-data.py)")
    g.add_argument("--skip-clumpify", action="store_true", help="跳过 Clumpify")
    g.add_argument("--dedup", action="store_true", help="启用 fastp 去重")

    # 去宿主
    g = p.add_argument_group("去宿主 (host_depletion.py)")
    g.add_argument("--kraken2_db", help="Kraken2 宿主库 (去宿主必需)")
    g.add_argument("--host_align_db", help="宿主比对索引 (去宿主必需)")
    g.add_argument("--tool", default="bowtie2", choices=["bowtie2", "hisat2", "minimap2"], help="比对工具")
    g.add_argument("--seq-type", default="rna-short", choices=["dna-short", "rna-short", "nanopore", "pacbio"])
    g.add_argument("--rrna", action="store_true", help="开启 Ribodetector 去 rRNA")
    g.add_argument("--confidence", type=float, default=0.4, help="Kraken2 置信度阈值")
    g.add_argument("--steps", default="kraken2,align,rrna", help="去宿主步骤")

    # 跳过
    g = p.add_argument_group("流程控制")
    g.add_argument("--skip-clean", action="store_true", help="跳过清洗")
    g.add_argument("--skip-depletion", action="store_true", help="跳过去宿主")
    g.add_argument("--force", action="store_true", help="强制重跑")
    g.add_argument("--threads", type=int, default=20, help="线程数")
    g.add_argument("--jobs", type=int, default=4, help="并行数")

    args = p.parse_args()
    script_dir = Path(__file__).parent.resolve()
    out = Path(args.output_dir).resolve()
    in_dir = Path(args.input).resolve()

    logger = setup_logger(str(out))
    logger.info("=" * 50)
    logger.info("数据预处理: %s → %s", in_dir, out)
    logger.info("  清洗: %s  去宿主: %s", not args.skip_clean, not args.skip_depletion)
    logger.info("=" * 50)

    reads_dir = str(in_dir)

    # ── Step 1: 清洗 ──
    if not args.skip_clean:
        clean_dir = out / "clean"
        parts = [
            f"python {script_dir / 'clean-data.py'}",
            f"-i {reads_dir}",
            f"-o {clean_dir}",
            f"-j {args.jobs}",
            f"-t {args.threads}",
        ]
        if args.skip_clumpify:
            parts.append("--skip-clumpify")
        if args.dedup:
            parts.append("--dedup")
        if args.force:
            parts.append("--force")
        if not run(' '.join(parts), logger, "clean-data"):
            sys.exit(1)

        # 更新 reads 指针: 优先 clumpify
        cl = clean_dir / "3.clumpify"
        fa = clean_dir / "2.fasta"
        reads_dir = str(cl if (cl.exists() and any(cl.iterdir())) else fa)
        logger.info("  清洗完成 → %s", reads_dir)
    else:
        logger.info("  [SKIP] 清洗")

    # ── Step 2: 去宿主 ──
    if not args.skip_depletion:
        if not args.kraken2_db or not args.host_align_db:
            logger.error("去宿主需要 --kraken2_db 和 --host_align_db")
            sys.exit(1)

        host_dir = out / "host"
        parts = [
            f"python {script_dir / 'host_depletion.py'}",
            f"--tool {args.tool}",
            f"--seq-type {args.seq_type}",
            f"--kraken2_index {args.kraken2_db}",
            f"--step2_index {args.host_align_db}",
            f"--input-dir {reads_dir}",
            f"--outdir {host_dir}",
            f"--jobs {args.jobs}",
            f"--threads {args.threads}",
            f"--logs_dir {host_dir}/logs",
            f"--steps {args.steps}",
            f"--confidence {args.confidence}",
            "--filter true",
        ]
        if args.rrna:
            parts.append("--rrna")
        if not run(' '.join(parts), logger, "host_depletion"):
            sys.exit(1)

        logger.info("  去宿主完成 → %s", host_dir)
    else:
        logger.info("  [SKIP] 去宿主")

    logger.info("=" * 50)
    logger.info("预处理完成!")
    if not args.skip_depletion:
        logger.info("  清洁 reads: %s", out / "host")
    logger.info("=" * 50)


if __name__ == "__main__":
    main()
