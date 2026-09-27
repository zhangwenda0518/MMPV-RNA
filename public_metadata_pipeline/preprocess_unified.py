#!/usr/bin/env python3
"""
preprocess_unified.py — 统一数据预处理 (从 .sra 一条命令到 clean reads)
=======================================================================
convert → clean → [hostref] → deplete → [bbnorm]，
一条命令从 .sra 到宿主剔除后 clean reads。

阶段:
  convert   SRA → FASTQ.GZ (sra2fastx.py, 可选)
  clean     FASTQ → QC'd FASTA (clean-data.py)
  hostref   宿主参考基因组下载 + 四索引构建 (可选, download_host_genome.py + build_hostbase.py)
  deplete   FASTA → 宿主剔除 + rRNA 去除 (host_depletion.py)
  bbnorm    k-mer 覆盖度归一化 (可选, --bbnorm; run_bbnorm.py, target=70 mindepth=2)

用法:
  # 从 SRA 开始全流程
  python preprocess_unified.py --sra-dir ./sra_files/ \
      --kraken2-db ~/database/host_db/kraken2 \
      --step2-index ~/database/host_db/bowtie2/host \
      --outdir ./preprocessed/

  # 已有 FASTQ, 跳过 convert
  python preprocess_unified.py --fastq-dir ./fastq_files/ \
      --kraken2-db ~/database/host_db/kraken2 \
      --step2-index ~/database/host_db/bowtie2/host \
      --outdir ./preprocessed/ --skip-convert

  # 加大深度样本共组装前的覆盖度归一化 (可选, 默认不跑)
  python preprocess_unified.py --fastq-dir ./fastq_files/ --outdir ./preprocessed/ \
      --skip-convert --kraken2-db ... --step2-index ... --bbnorm

输出:
  outdir/
    00_Converted/      ← .sra → .fastq.gz
    00a_CleanData/     ← fastp + clumpify + FASTA
    00b_HostDepletion/ ← host-depleted + rRNA-removed reads
    00c_BBnorm/        ← 覆盖度归一化 reads (仅 --bbnorm 时产出)
    preprocess.log     ← 全流程日志

目录命名与 data_preprocessing_pipeline / virome_discovery_pipeline 一致
(00a_CleanData → 00b_HostDepletion → [00c_BBnorm])，
因此产出可直接作为 `virome_pipeline.py --input_reads` 或
`auto_known_virus.py --reads_dir` 的输入。
"""

import argparse
import os
import sys
import subprocess
import logging
from pathlib import Path
from datetime import datetime

SCRIPT_DIR = Path(__file__).resolve().parent
PREPROC_DIR = SCRIPT_DIR.parent / "data_preprocessing_pipeline"

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根); 环境变量会传给子进程 ②
_REPO_ROOT = SCRIPT_DIR.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
from mmpv_common.io_layout import dir_name as layout_dir_name, normalize_layout_env


def setup_logger(outdir):
    logger = logging.getLogger("Preprocess")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    os.makedirs(outdir, exist_ok=True)
    log_file = os.path.join(outdir, "preprocess.log")

    fmt = logging.Formatter("[%(asctime)s] %(levelname)s %(message)s", datefmt="%H:%M:%S")
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger


def run_cmd(cmd, logger, desc, timeout=None):
    logger.info("[%s] 执行中...", desc)
    logger.debug("  CMD: %s", cmd)
    try:
        result = subprocess.run(cmd, shell=True, text=True,
                                capture_output=True, timeout=timeout)
        for line in (result.stdout or "").splitlines():
            logger.debug("  %s", line)
        if result.returncode != 0:
            logger.error("[%s] 失败 (exit=%d)\n%s", desc, result.returncode,
                          (result.stderr or "")[-500:])
            return False
        logger.info("[%s] 完成", desc)
        return True
    except subprocess.TimeoutExpired:
        logger.error("[%s] 超时", desc)
        return False
    except Exception as e:
        logger.error("[%s] 异常: %s", desc, e)
        return False


def has_reads(d):
    """目录存在且含可直接归一化的 reads (bbnorm 输入检查)"""
    d = Path(d)
    if not d.is_dir():
        return False
    pats = ("*.fastq.gz", "*.fq.gz", "*.fasta.gz", "*.fa.gz", "*.fasta", "*.fa")
    return any(next(d.glob(p), None) is not None for p in pats)


def stage_convert(sra_dir, outdir, logger, jobs=2, threads=6):
    """SRA → FASTQ.GZ"""
    from sra2fastx import convert_all
    logger.info("=== [convert] SRA → FASTQ.GZ ===")
    if not any(Path(sra_dir).rglob("*.sra")):
        logger.info("  未找到 .sra 文件, 跳过")
        return True

    success, failures = convert_all(sra_dir, jobs=jobs, threads=threads, progress=True)
    if failures:
        logger.error("  %d 个转换失败", len(failures))
        return False
    logger.info("  %d 个 .sra 转换完成 → %s", success, outdir)
    return True


def stage_clean(fastq_dir, outdir, logger,
                fastp_threads=4, jobs=2, clumpify_memory="10g",
                dedup=False, skip_clumpify=False):
    """FASTQ → QC'd FASTA (fastp + clumpify + seqkit)"""
    logger.info("=== [clean] FASTQ → QC'd FASTA ===")
    cmd_parts = [
        f"python {PREPROC_DIR / 'clean-data.py'}",
        f"-i {fastq_dir}",
        f"-o {outdir}",
        f"--fastp-threads {fastp_threads}",
        f"--clumpify-memory {clumpify_memory}",
        f"-j {jobs}",
    ]
    if dedup:
        cmd_parts.append("--dedup")
    if skip_clumpify:
        cmd_parts.append("--skip-clumpify")
    return run_cmd(" ".join(cmd_parts), logger, "clean")


def stage_hostref(logger, host_species, host_taxid, host_fasta,
                  host_work_dir, ncbi_api, tool, seq_type, threads,
                  host_download_source="auto"):
    """宿主参考构建: 直接调 download_host_genome.py (四通道) + build_hostbase.py, 不嵌套子管道。

    返回 (kraken2_db, step2_index_prefix)；失败返回 None。
    """
    logger.info("=== [hostref] 宿主参考基因组获取与索引构建 ===")
    work = Path(host_work_dir)
    hostdb = work / "hostdb"

    # genome: 用户 FASTA 优先, 否则四通道下载
    if host_fasta:
        fasta = Path(host_fasta)
        if not fasta.is_file():
            logger.error("[hostref] 宿主基因组不存在: %s", fasta)
            return None
    else:
        fasta = work / "genome" / "all.genome.uniq.fasta"
        if not fasta.is_file():
            genome_dir = work / "genome"
            genome_dir.mkdir(parents=True, exist_ok=True)
            cmd = (f"python {SCRIPT_DIR / 'download_host_genome.py'} "
                   f"--species {host_species} --taxid {host_taxid} "
                   f"--outdir {genome_dir} --source {host_download_source}")
            if ncbi_api:
                cmd += f" --ncbi-api {ncbi_api}"
            if not run_cmd(cmd, logger, "hostref-genome"):
                return None
    if not fasta.is_file():
        logger.error("[hostref] 未找到宿主基因组 FASTA: %s", fasta)
        return None

    # hostdb: Kraken2/Bowtie2/HISAT2/Minimap2 四索引
    cmd = (f"python {SCRIPT_DIR / 'build_hostbase.py'} "
           f"--tool kraken2,bowtie2,hisat2,minimap2 "
           f"--input {fasta} --output {hostdb} "
           f"--threads {threads} "
           f"--seq-type dna-short,rna-short,nanopore,pacbio "
           f"--k2-libs archaea,bacteria,plasmid,fungi,protozoa,UniVec "
           f"--taxid {host_taxid}")
    if not run_cmd(cmd, logger, "hostdb"):
        return None
    kraken2_db = hostdb / "kraken2"
    step2_index = hostdb / tool / "host"
    if not (kraken2_db / "hash.k2d").is_file():
        logger.error("[hostref] Kraken2 库缺失: %s", kraken2_db)
        return None
    logger.info("  宿主参考就绪: kraken2=%s, %s 索引前缀=%s", kraken2_db, tool, step2_index)
    return str(kraken2_db), str(step2_index)


def stage_deplete(fasta_dir, outdir, logger,
                  kraken2_db, step2_index,
                  tool="bowtie2", seq_type="rna-short",
                  threads=4, jobs=4, rrna=True,
                  tmp=None):
    """FASTA → host-depleted + rRNA-removed reads"""
    logger.info("=== [deplete] 宿主剔除 + rRNA 去除 ===")
    cmd_parts = [
        f"python {PREPROC_DIR / 'host_depletion.py'}",
        f"-I {fasta_dir}",
        f"-O {outdir}",
        f"--tool {tool}",
        f"--seq-type {seq_type}",
        f"-k {kraken2_db}",
        f"-x {step2_index}",
        f"--threads {threads}",
        f"--jobs {jobs}",
    ]
    if rrna:
        cmd_parts.append("--rrna")
    if tmp:
        cmd_parts.extend(["-T", tmp])
    return run_cmd(" ".join(cmd_parts), logger, "deplete")


def main():
    parser = argparse.ArgumentParser(
        description="统一数据预处理: convert → clean → [hostref] → deplete [→ bbnorm 可选]",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 从 SRA 到 clean reads
  python preprocess_unified.py --sra-dir ./sra/ --outdir ./out/ \\
      --kraken2-db ~/database/host_db/kraken2 \\
      --step2-index ~/database/host_db/bowtie2/host

  # 已有 FASTQ, 跳过 convert
  python preprocess_unified.py --fastq-dir ./fq/ --outdir ./out/ \\
      --skip-convert --kraken2-db ~/database/host_db/kraken2 \\
      --step2-index ~/database/host_db/bowtie2/host
        """)

    g = parser.add_argument_group("输入 (二选一)")
    g.add_argument("--sra-dir", help="包含 .sra 文件的目录")
    g.add_argument("--fastq-dir", help="包含 .fastq.gz 文件的目录 (跳过 convert)")

    g = parser.add_argument_group("输出")
    g.add_argument("--outdir", "-o", default="./preprocessed", help="输出根目录 (默认: ./preprocessed)")
    g.add_argument("--io-layout", choices=["legacy", "standard"], default=None,
                   help="I/O 目录布局: legacy=v3.0 现行名 (默认); standard=统一编号"
                        "(doc/IO_LAYOUT_DESIGN.md); 亦可 MMPV_IO_LAYOUT 环境变量")

    g = parser.add_argument_group("阶段控制")
    g.add_argument("--skip-convert", action="store_true", help="跳过 SRA 转换")
    g.add_argument("--skip-clean", action="store_true", help="跳过清洗 (fastp+clumpify)")
    g.add_argument("--skip-deplete", action="store_true", help="跳过宿主剔除")
    g.add_argument("--bbnorm", action="store_true",
                   help="可选: 去宿主后 BBNorm 覆盖度归一化 → 00c_BBnorm/ (默认不运行; "
                        "需 deplete 产物, 与 --skip-deplete 同用时要求 00b_HostDepletion/ 已存在)")

    g = parser.add_argument_group("宿主数据库 (已有库时指定)")
    g.add_argument("--kraken2-db", help="Kraken2 宿主数据库路径 (已有库时指定; 否则可用 --host-* 自动构建)")
    g.add_argument("--step2-index", help="第二步比对索引前缀 (bowtie2/hisat2/minimap2; 已有索引时指定)")

    g = parser.add_argument_group("宿主参考构建 (可选: 自动构建索引后接续 deplete)")
    g.add_argument("--host-species", help="宿主物种拉丁学名 (与 --host-taxid 同用)")
    g.add_argument("--host-taxid", type=int, help="宿主 NCBI Taxonomy ID")
    g.add_argument("--host-fasta", help="已有宿主基因组 FASTA (跳过 NCBI 下载)")
    g.add_argument("--host-work-dir", help="宿主参考构建输出目录 (默认: <outdir>/host_reference)")
    g.add_argument("--ncbi-api", help="NCBI API Key")
    g.add_argument("--host-download-source", default="auto",
                   choices=["auto", "datasets", "ngd", "ftp", "gget"],
                   help="宿主基因组下载通道: auto=datasets→ngd→FTP 自动回退; gget 走 Ensembl")
    g.add_argument("--tool", default="bowtie2", choices=["bowtie2", "hisat2", "minimap2"], help="比对工具 (默认: bowtie2)")
    g.add_argument("--seq-type", default="rna-short", choices=["dna-short", "rna-short", "nanopore", "pacbio"])

    g = parser.add_argument_group("clean 阶段")
    g.add_argument("--dedup", action="store_true", help="fastp 内置去重")
    g.add_argument("--skip-clumpify", action="store_true", help="跳过 clumpify 聚类重排")
    g.add_argument("--clumpify-memory", default="10g", help="clumpify 堆内存 (默认: 10g)")

    g = parser.add_argument_group("计算资源")
    g.add_argument("--threads", "-t", type=int, default=40, help="总线程 (默认: 40)")
    g.add_argument("--jobs", "-j", type=int, default=4, help="并行任务数 (默认: 4)")
    g.add_argument("--tmp", help="临时目录")

    args = parser.parse_args()

    # 布局解析 (CLI > MMPV_IO_LAYOUT > legacy) 并写回环境供子进程继承
    normalize_layout_env(getattr(args, "io_layout", None))

    # 输入校验
    if not args.sra_dir and not args.fastq_dir:
        parser.error("必须指定 --sra-dir 或 --fastq-dir")
    if args.sra_dir and not os.path.isdir(args.sra_dir):
        parser.error(f"SRA 目录不存在: {args.sra_dir}")
    if args.fastq_dir and not os.path.isdir(args.fastq_dir):
        parser.error(f"FASTQ 目录不存在: {args.fastq_dir}")

    wants_hostref = bool(args.host_species or args.host_fasta)
    has_host_db = bool(args.kraken2_db and args.step2_index)
    if not args.skip_deplete and not has_host_db and not wants_hostref:
        parser.error("需要 --kraken2-db + --step2-index（已有库），"
                     "或提供 --host-species/--host-taxid、--host-fasta 以自动构建宿主参考")
    if wants_hostref and not args.host_taxid:
        parser.error("--host-taxid 必填（与 --host-species 或 --host-fasta 配合）")

    outdir = Path(args.outdir)
    convert_dir = outdir / "00_Converted"
    clean_dir = outdir / layout_dir_name('d_clean')
    deplete_dir = outdir / layout_dir_name('d_hostdep')

    # bbnorm 的输入是去宿主后的 reads: 跳过 deplete 时必须已有现成产物, 否则提前报错
    # (不要等到跑完 convert/clean 才发现 bbnorm 无输入)
    if args.bbnorm and args.skip_deplete and not has_reads(deplete_dir):
        parser.error(f"--bbnorm 需要去宿主后的 reads, 但 {deplete_dir} 不存在或为空; "
                     "请去掉 --skip-deplete 让本条命令生成, 或先单独跑完 deplete 阶段")

    logger = setup_logger(str(outdir))
    logger.info("=" * 50)
    logger.info("统一数据预处理")
    logger.info("  模式: %s", "SRA→" if args.sra_dir else "FASTQ→")
    logger.info("  流程: %s",
                " → ".join([s for s, skip in [
                    ("convert", args.skip_convert or not args.sra_dir),
                    ("clean", args.skip_clean),
                    ("hostref", not wants_hostref or args.skip_deplete),
                    ("deplete", args.skip_deplete),
                    ("bbnorm", not args.bbnorm),
                ] if not skip]))
    logger.info("  输出: %s", outdir)
    logger.info("=" * 50)

    start_time = datetime.now()
    failed = []

    # ── Stage 1: convert ──
    if not args.skip_convert and args.sra_dir:
        convert_dir.mkdir(parents=True, exist_ok=True)
        if not stage_convert(args.sra_dir, str(convert_dir), logger,
                             jobs=args.jobs, threads=min(args.threads, 12)):
            failed.append("convert")
    elif args.skip_convert:
        logger.info("[convert] 跳过")

    # ── 确定 clean 输入 ──
    if args.fastq_dir:
        clean_input = args.fastq_dir
    elif not args.skip_convert and args.sra_dir:
        clean_input = str(convert_dir)
    else:
        clean_input = str(args.sra_dir)

    # ── Stage 2: clean ──
    if not args.skip_clean:
        clean_dir.mkdir(parents=True, exist_ok=True)
        # clean-data.py 输出到 1.fastp/, 2.fasta/, 3.clumpify/ 三个子目录
        if not stage_clean(clean_input, str(clean_dir), logger,
                           fastp_threads=min(args.threads, 8),
                           jobs=args.jobs,
                           clumpify_memory=args.clumpify_memory,
                           dedup=args.dedup,
                           skip_clumpify=args.skip_clumpify):
            failed.append("clean")
    else:
        logger.info("[clean] 跳过")

    # ── Stage 3: hostref (可选, 自动构建宿主参考) ──
    if wants_hostref and not args.skip_deplete:
        host_work = args.host_work_dir or str(
            outdir / os.path.split(layout_dir_name("h_genome"))[0])  # legacy=host_reference / standard=05_HostRef
        built = stage_hostref(logger, args.host_species, args.host_taxid,
                              args.host_fasta, host_work, args.ncbi_api,
                              args.tool, args.seq_type, threads=args.threads,
                              host_download_source=args.host_download_source)
        if built:
            args.kraken2_db, args.step2_index = built
        else:
            failed.append("hostref")
            logger.error("[hostref] 宿主参考构建失败, deplete 阶段终止")
            args.skip_deplete = True

    # ── Stage 4: deplete ──
    if not args.skip_deplete:
        deplete_dir.mkdir(parents=True, exist_ok=True)
        # 优先使用 clumpify 产出, 回退到 fasta
        fasta_input = clean_dir / "3.clumpify"
        if not fasta_input.is_dir() or not any(fasta_input.iterdir()):
            fasta_input = clean_dir / "2.fasta"
        if not str(fasta_input):
            logger.error("[deplete] clean 阶段无有效输出, 终止")
            sys.exit(1)
        logger.info("  deplete 输入: %s", fasta_input)

        if not stage_deplete(str(fasta_input), str(deplete_dir), logger,
                             kraken2_db=args.kraken2_db,
                             step2_index=args.step2_index,
                             tool=args.tool,
                             seq_type=args.seq_type,
                             threads=args.threads,
                             jobs=args.jobs,
                             rrna=True,
                             tmp=args.tmp):
            failed.append("deplete")
    else:
        logger.info("[deplete] 跳过")

    # ── Stage: bbnorm (可选, 仅 --bbnorm 显式启用; 默认流程不含) ──
    if args.bbnorm:
        if "deplete" in failed:
            logger.warning("[bbnorm] 因 deplete 失败而跳过")
        elif not has_reads(deplete_dir):
            logger.error("[bbnorm] 输入目录无 reads: %s (需先完成 deplete 阶段)", deplete_dir)
            failed.append("bbnorm")
        else:
            bbnorm_dir = outdir / layout_dir_name('d_bbnorm')
            bbnorm_dir.mkdir(parents=True, exist_ok=True)
            cmd = (f"python {PREPROC_DIR / 'run_bbnorm.py'} -i {deplete_dir} "
                   f"-o {bbnorm_dir} -t {args.threads} -j {args.jobs}")
            if not run_cmd(cmd, logger, "bbnorm"):
                failed.append("bbnorm")
            else:
                logger.info("  Normalized reads → %s", bbnorm_dir)

    elapsed = datetime.now() - start_time
    logger.info("=" * 50)
    if failed:
        logger.error("失败阶段: %s", ", ".join(failed))
        sys.exit(1)
    logger.info("全流程完成 | 耗时: %s", elapsed)
    logger.info("  convert → %s", convert_dir)
    logger.info("  clean   → %s", clean_dir)
    logger.info("  deplete → %s (下游 --input_reads)", deplete_dir)
    logger.info("=" * 50)


if __name__ == "__main__":
    main()
