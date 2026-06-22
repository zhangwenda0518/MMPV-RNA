#!/usr/bin/env python3
"""
data_preprocessing.py — 数据预处理独立脚本 (清洗 + 去宿主)
============================================================
  Stage 1: Clean  → Fastp质控 → Seqkit格式转换 → Clumpify光学去重
  Stage 2: Deplete → Kraken2分类 → Bowtie2/Minimap2精准去宿主 → rRNA剔除

Usage:
  # 只清洗
  python data_preprocessing.py --stage clean --input_reads raw/ --output_dir out/ -t 40 -j 10

  # 只去宿主 (自动从 clean 输出读取)
  python data_preprocessing.py --stage deplete --output_dir out/ --host_db host_db/ -t 40 -j 10

  # 全跑
  python data_preprocessing.py --stage all --input_reads raw/ --output_dir out/ --host_db host_db/ -t 40 -j 10
"""

import argparse, logging, os, shutil, subprocess, sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

def setup_logger(output_dir, level='INFO'):
    logger = logging.getLogger("DataPrep")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    os.makedirs(output_dir, exist_ok=True)
    ch = logging.StreamHandler()
    ch.setLevel(getattr(logging, level.upper(), logging.INFO))
    ch.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s %(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(ch)
    fh = logging.FileHandler(os.path.join(output_dir, 'preprocessing.log'), encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s %(message)s'))
    logger.addHandler(fh)
    return logger

def run_cmd(cmd, log, step):
    log.info("[%s] %s", step, cmd[:200])
    try:
        subprocess.run(cmd, shell=True, check=True)
        return True
    except subprocess.CalledProcessError as e:
        log.error("[%s] 失败 (exit=%d)", step, e.returncode)
        return False

def find_host_db(args):
    """从 --host_db 自动推导 kraken2 和比对索引路径"""
    kraken2_db = args.kraken2_db
    host_align_db = args.host_align_db
    if args.host_db:
        hdb = Path(args.host_db)
        if not kraken2_db:
            for sub in ['kraken2', 'kraken2_db', 'kraken']:
                if (hdb / sub).is_dir(): kraken2_db = str(hdb / sub); break
        if not host_align_db:
            aligner = args.aligner
            for sub in [aligner, f'{aligner}_index', 'align']:
                d = hdb / sub
                if not d.is_dir(): continue
                for prefix in ['host', 'index', 'genome']:
                    test = d / prefix
                    if aligner == 'bowtie2' and ((test.parent / (prefix + '.1.bt2')).is_file() or (test.parent / (prefix + '.1.bt2l')).is_file()):
                        host_align_db = str(test); break
                    if aligner == 'hisat2' and (test.parent / (prefix + '.1.ht2')).is_file():
                        host_align_db = str(test); break
                    if aligner == 'minimap2':
                        mmi = test.parent / f'{prefix}_{args.seq_type}.mmi'
                        if mmi.is_file(): host_align_db = str(mmi); break
                if host_align_db: break
    return kraken2_db, host_align_db

def run_clean(args, log, reads_dir):
    """Step 1: Fastp + Seqkit + Clumpify"""
    clean_dir = Path(args.output_dir) / '00a_CleanData'
    clean_dir.mkdir(parents=True, exist_ok=True)

    log.info("=" * 50)
    log.info("[1/2] 数据清洗: Fastp → Seqkit → Clumpify")

    clean_script = SCRIPT_DIR / 'clean-data.py'
    parts = [
        f"python {clean_script}",
        f"--input {reads_dir}",
        f"--output {clean_dir}",
        f"--fastp-threads {args.threads}",
        f"--jobs {args.jobs}",
    ]
    if args.skip_clumpify: parts.append("--skip-clumpify")
    if args.force: parts.append("--force")
    if args.dedup: parts.append("--dedup")
    if args.clumpify_memory: parts.append(f"--clumpify-memory {args.clumpify_memory}")
    if args.no_compress: parts.append("--no-compress")

    if not run_cmd(' '.join(parts), log, "Clean"): sys.exit(1)

    # 更新 reads 指针
    cl = clean_dir / '3.clumpify'
    fa = clean_dir / '2.fasta'
    reads_dir = cl if (cl.exists() and any(cl.iterdir())) else fa
    log.info("  Reads → %s", reads_dir)
    return str(reads_dir)

def run_deplete(args, log, reads_dir):
    """Step 2: Kraken2 + Align + rRNA removal"""
    hostdep_dir = Path(args.output_dir) / '00b_HostDepletion'
    hostdep_dir.mkdir(parents=True, exist_ok=True)

    kraken2_db, host_align_db = find_host_db(args)
    if not kraken2_db:
        log.error("致命: 需要 --kraken2_db 或 --host_db"); sys.exit(1)
    if not host_align_db:
        log.error("致命: 需要 --host_align_db 或 --host_db"); sys.exit(1)

    log.info("=" * 50)
    log.info("[2/2] 去宿主: Kraken2 → %s → rRNA", args.aligner)

    deplete_script = SCRIPT_DIR / 'host_depletion.py'
    parts = [
        f"python {deplete_script}",
        f"--tool {args.aligner}",
        f"--seq-type {args.seq_type}",
        f"--kraken2_index {kraken2_db}",
        f"--step2_index {host_align_db}",
        f"--input-dir {reads_dir}",
        f"--outdir {hostdep_dir}",
        f"--jobs {args.jobs}",
        f"--threads {args.threads}",
        f"--logs_dir {hostdep_dir}/logs",
        "--filter true",
    ]
    if args.rrna:
        parts.append("--rrna")
        parts.append(f"--rrna_tool {args.rrna_tool}")
        if args.silva_index: parts.append(f"--silva_index {args.silva_index}")
    if args.force: parts.append("--force")
    if args.tmp_dir: parts.append(f"--tmp {args.tmp_dir}")
    parts.append(f"--confidence {args.kraken2_confidence}")
    if args.keep_rrna: parts.append("--keep_rrna")

    if not run_cmd(' '.join(parts), log, "Deplete"): sys.exit(1)
    log.info("  Depleted reads → %s", hostdep_dir)
    return str(hostdep_dir)

def main():
    p = argparse.ArgumentParser(description="数据预处理独立脚本 — 清洗 + 去宿主",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--stage', default=['all'], nargs='+',
                   choices=['all','clean','deplete'],
                   help='运行阶段 (可多个: --stage clean deplete)')
    p.add_argument('--input_reads', help='输入 FASTQ/FASTA 目录 (clean 阶段必需)')
    p.add_argument('--output_dir', '-o', required=True, help='输出根目录')

    g = p.add_argument_group('数据库 (去宿主)')
    g.add_argument('--host_db', help='宿主数据库根目录 (自动查找 kraken2/ bowtie2/ hisat2/ minimap2/)')
    g.add_argument('--kraken2_db', help='Kraken2 宿主库')
    g.add_argument('--host_align_db', help='宿主比对索引前缀 (bowtie2/minimap2)')

    g = p.add_argument_group('工具与算法')
    g.add_argument('--aligner', default='bowtie2', choices=['bowtie2','hisat2','minimap2'], help='比对工具 (默认: bowtie2)')
    g.add_argument('--seq_type', default='rna-short', choices=['dna-short','rna-short','nanopore','pacbio'])
    g.add_argument('--rrna', action='store_true', help='开启 rRNA 剔除')
    g.add_argument('--rrna_tool', default='ribodetector', choices=['ribodetector','silva'], help='rRNA 工具')
    g.add_argument('--silva_index', help='SILVA Bowtie2 索引前缀 (--rrna_tool silva 时必需)')
    g.add_argument('--kraken2_confidence', type=float, default=0.0, help='Kraken2 置信度阈值 (默认: 0.0)')

    g = p.add_argument_group('清洗参数')
    g.add_argument('--skip_clumpify', action='store_true', help='跳过 Clumpify 光学去重')
    g.add_argument('--dedup', action='store_true', help='fastp 自带去重')
    g.add_argument('--clumpify_memory', default='10g', help='clumpify 内存 (默认: 10g)')
    g.add_argument('--no_compress', action='store_true', help='输出不压缩')

    g = p.add_argument_group('计算资源')
    g.add_argument('--threads', '-t', type=int, default=20, help='线程数 (默认: 20)')
    g.add_argument('--jobs', '-j', type=int, default=2, help='并行样本数 (默认: 2)')

    g = p.add_argument_group('流程控制')
    g.add_argument('--force', action='store_true', help='强制重跑')
    g.add_argument('--keep_rrna', action='store_true', help='保留 rRNA reads')
    g.add_argument('--tmp_dir', help='临时目录')

    args = p.parse_args()
    stages = set(args.stage)
    _all = 'all' in stages

    log = setup_logger(args.output_dir)

    log.info("=" * 50)
    log.info("Data Preprocessing Pipeline")
    log.info("  Stage:  %s", ','.join(sorted(stages)))
    log.info("  Output: %s", args.output_dir)
    log.info("=" * 50)

    # 确定 reads 目录
    if args.input_reads:
        reads_dir = args.input_reads
    elif 'deplete' in stages and not _all:
        # deplete standalone: 从 clean 输出读取
        cl = Path(args.output_dir) / '00a_CleanData' / '3.clumpify'
        fa = Path(args.output_dir) / '00a_CleanData' / '2.fasta'
        reads_dir = str(cl) if cl.exists() and any(cl.iterdir()) else str(fa)
        if not Path(reads_dir).exists():
            log.error("未找到 clean 输出, 请先 --stage clean 或指定 --input_reads")
            sys.exit(1)
        log.info("  Reads → %s (auto)", reads_dir)
    else:
        log.error("需要 --input_reads")
        sys.exit(1)

    # 执行
    if _all or 'clean' in stages:
        reads_dir = run_clean(args, log, reads_dir)
    if _all or 'deplete' in stages:
        run_deplete(args, log, reads_dir)

    log.info("=" * 50)
    log.info("预处理完成!")
    log.info("  清洗:   %s/00a_CleanData", args.output_dir)
    log.info("  去宿主: %s/00b_HostDepletion", args.output_dir)
    log.info("=" * 50)

if __name__ == '__main__':
    main()
