#!/usr/bin/env python3
"""
data_preprocessing.py — 数据预处理独立脚本 (清洗 + 去宿主)
============================================================
  Stage 1: Clean  → Fastp质控 → Seqkit格式转换 → Clumpify聚类重排
  Stage 2: Deplete → Kraken2分类 → Bowtie2/Minimap2精准去宿主 → rRNA剔除
  依赖: 同目录 clean-data.py ＋ host_depletion.py

Usage:
  # 只清洗
  python data_preprocessing.py --stage clean --input_reads raw/ --output_dir out/ -t 40 -j 10

  # 只去宿主 (自动从 clean 输出读取)
  python data_preprocessing.py --stage deplete --output_dir out/ --host_db host_db/ -t 40 -j 10

  # 全跑
  python data_preprocessing.py --stage all --input_reads raw/ --output_dir out/ --host_db host_db/ -t 40 -j 10
"""

import argparse, logging, os, re, shutil, subprocess, sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根); 布局经环境变量向子进程传播
_REPO_ROOT = SCRIPT_DIR.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
from mmpv_common.io_layout import dir_name as layout_dir_name, normalize_layout_env

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
    """从 --host_db 自动推导 kraken2 和比对索引路径.

    --host_db 可指向: hostdb/ 目录本身、host_reference 根 (自动下钻 hostdb/)、
    或直接的 kraken2 库目录 (含 hash.k2d)。
    """
    kraken2_db = args.kraken2_db
    host_align_db = args.host_align_db
    if args.host_db:
        hdb = Path(args.host_db)
        if not kraken2_db and (hdb / 'hash.k2d').is_file():
            kraken2_db = str(hdb)
        if (hdb / 'hostdb').is_dir():
            hdb = hdb / 'hostdb'
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

def ensure_hostdb(args, log):
    """无现成库时现场构建宿主参考 (--host_fasta 或 --host_species+--host_taxid).

    直接调用 download_host_genome.py (四通道下载) + build_hostbase.py (只建 kraken2 +
    当前比对工具), 结果写入 {output_dir}/host_reference/hostdb/, 并回填 args.host_db。
    """
    host_reference = Path(args.output_dir) / 'host_reference'
    hostdb = host_reference / 'hostdb'
    if hostdb.is_dir() and any(hostdb.iterdir()) and not args.force:
        log.info("  复用已有宿主参考: %s (--force 重建)", hostdb)
        args.host_db = str(hostdb)
        return

    bin_dir = SCRIPT_DIR.parent / 'public_metadata_pipeline'
    genome_dir = host_reference / 'genome'
    fasta = Path(args.host_fasta) if args.host_fasta else genome_dir / 'all.genome.uniq.fasta'

    if not fasta.is_file():
        if args.host_fasta:
            log.error("宿主基因组不存在: %s", args.host_fasta)
            sys.exit(1)
        log.info("[前置] 宿主基因组下载 (download_host_genome.py 四通道)")
        parts = [
            f"python {bin_dir / 'download_host_genome.py'}",
            f"--species {args.host_species}",
            f"--taxid {args.host_taxid}",
            f"--outdir {genome_dir}",
            "--source auto",
        ]
        if not run_cmd(' '.join(parts), log, "GenomeDown"):
            sys.exit(1)
    if not fasta.is_file():
        log.error("未找到宿主基因组 FASTA: %s", fasta)
        sys.exit(1)

    log.info("[前置] 宿主索引构建 (build_hostbase.py: kraken2 + %s)", args.aligner)
    parts = [
        f"python {bin_dir / 'build_hostbase.py'}",
        f"--tool kraken2,{args.aligner}",
        f"--input {fasta}",
        f"--output {hostdb}",
        f"--threads {args.threads}",
        f"--seq-type {args.seq_type}",
        "--k2-libs archaea,bacteria,plasmid,fungi,protozoa,UniVec",
    ]
    if args.host_taxid:
        parts.append(f"--taxid {args.host_taxid}")
    if not run_cmd(' '.join(parts), log, "HostDB"):
        sys.exit(1)
    args.host_db = str(hostdb)


def run_clean(args, log, reads_dir):
    """Step 1: Fastp + Seqkit + Clumpify"""
    clean_dir = Path(args.output_dir) / layout_dir_name('d_clean')
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
    hostdep_dir = Path(args.output_dir) / layout_dir_name('d_hostdep')
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

def run_report(args, log):
    """Step 3: 汇总各步统计 → summary TSV + assembly_ready.list + HTML 报告"""
    log.info("=" * 50)
    log.info("[3/3] Report: 汇总 + 交接清单 (preprocess_report.py)")
    # 显式传 clean/deplete 目录 (preprocess_report 的默认值是 legacy 名, 布局感知须显式)
    parts = [
        f"python {SCRIPT_DIR / 'preprocess_report.py'}",
        f"--output-dir {args.output_dir}",
        f"--clean-dir {Path(args.output_dir) / layout_dir_name('d_clean')}",
        f"--deplete-dir {Path(args.output_dir) / layout_dir_name('d_hostdep')}",
        f"--warn-retained {args.warn_retained}",
    ]
    ok = run_cmd(' '.join(parts), log, "Report")
    if ok:
        log.info("  汇总表:   %s/preprocessing_summary.tsv", args.output_dir)
        log.info("  交接清单: %s/assembly_ready.list", args.output_dir)
        log.info("  HTML:     %s/preprocessing_report.html", args.output_dir)
    return ok


def cleanup_00a_after_deplete(output_dir, log):
    """
    00b 完成后按样本删除 00a 中对应的测序数据文件。
    保留目录结构、日志 (clean.log / logs/) 和 fastp JSON 报告。
    不改动 00b 的任何数据，失败不影响管道继续。
    """
    clean_dir = Path(output_dir) / layout_dir_name('d_clean')
    if not clean_dir.is_dir():
        return

    hostdep_dir = Path(output_dir) / layout_dir_name('d_hostdep')
    if not hostdep_dir.is_dir() or not any(hostdep_dir.iterdir()):
        log.warning("  00b 目录不存在或为空, 跳过 00a 清理 (安全保护)")
        return

    # 从 00b 提取已完成的样本名: ERR2040117_clean_1.fa.gz → ERR2040117
    done_samples = set()
    for f in hostdep_dir.iterdir():
        if not f.is_file():
            continue
        name = f.name
        for pat in [r'^(.+?)_clean_[12]\.fa\.gz$', r'^(.+?)_clean_[12]\.fq\.gz$',
                    r'^(.+?)_[12]\.fa\.gz$', r'^(.+?)_[12]\.fq\.gz$']:
            m = re.match(pat, name)
            if m:
                done_samples.add(m.group(1))
                break

    if not done_samples:
        log.warning("  00b 中无法识别样本名, 跳过 00a 清理")
        return

    # 00a 的数据在子目录里: 1.fastp_tmp/, 2.fasta/, 3.clumpify/
    data_subdirs = ['1.fastp_tmp', '2.fasta', '3.clumpify']
    deleted_count = 0
    total_freed = 0

    for sub in data_subdirs:
        sub_path = clean_dir / sub
        if not sub_path.is_dir():
            continue
        for f in sub_path.iterdir():
            if not f.is_file():
                continue
            # 提取样本名并匹配
            name = f.name
            sample = None
            for pat in [r'^(.+?)_clean', r'^(.+?)_[12]\.', r'^(.+?)\.']:
                m = re.match(pat, name)
                if m:
                    sample = m.group(1)
                    break
            if sample and sample in done_samples:
                try:
                    size = f.stat().st_size
                    f.unlink()
                    total_freed += size
                    deleted_count += 1
                except Exception as e:
                    log.warning("  删除 %s 失败: %s", f.name, e)

    if deleted_count > 0:
        log.info("  00a 已清理 %d 个样本文件 (释放 %.1fG), 目录结构和日志保留",
                 deleted_count, total_freed / (1024**3))
    else:
        log.info("  00a 无需清理 (无可匹配的样本文件)")


def main():
    p = argparse.ArgumentParser(description="数据预处理独立脚本 — 清洗 + 去宿主",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--stage', default=['all'], nargs='+',
                   choices=['all','clean','deplete','bbnorm','report'],
                   help='运行阶段 (可多个: --stage clean deplete report; bbnorm 仅显式指定时运行, 不含于 all)')
    p.add_argument('--input_reads', help='输入 FASTQ/FASTA 目录 (clean 阶段必需)')
    p.add_argument('--output_dir', '-o', required=True, help='输出根目录')
    p.add_argument('--io-layout', choices=['legacy', 'standard'], default=None,
                   help='I/O 目录布局: legacy=v3.0 现行名 (默认); standard=统一编号'
                        '(doc/IO_LAYOUT_DESIGN.md); 亦可 MMPV_IO_LAYOUT 环境变量')

    g = p.add_argument_group('数据库 (去宿主)')
    g.add_argument('--host_db', help='宿主数据库目录 (支持 hostdb/ 本身、host_reference 根、kraken2 库目录, 自动查找子索引)')
    g.add_argument('--kraken2_db', help='Kraken2 宿主库 (与 --host_align_db 同用, 精确指定)')
    g.add_argument('--host_align_db', help='宿主比对索引前缀 (bowtie2/minimap2)')
    g.add_argument('--host_fasta', help='宿主基因组 FASTA (无现成库时现场构建到 {output_dir}/host_reference/)')
    g.add_argument('--host_species', help='宿主物种拉丁学名 (与 --host_taxid 同用, 现场构建)')
    g.add_argument('--host_taxid', type=int, help='宿主 NCBI Taxonomy ID')

    g = p.add_argument_group('工具与算法')
    g.add_argument('--aligner', default='bowtie2', choices=['bowtie2','hisat2','minimap2'], help='比对工具 (默认: bowtie2)')
    g.add_argument('--seq_type', default='rna-short', choices=['dna-short','rna-short','nanopore','pacbio'])
    g.add_argument('--rrna', action='store_true', help='开启 rRNA 剔除')
    g.add_argument('--rrna_tool', default='ribodetector', choices=['ribodetector','silva'], help='rRNA 工具')
    g.add_argument('--silva_index', help='SILVA Bowtie2 索引前缀 (--rrna_tool silva 时必需)')
    g.add_argument('--kraken2_confidence', type=float, default=0.2, help='Kraken2 置信度阈值 (默认: 0.2)')

    g = p.add_argument_group('清洗参数')
    g.add_argument('--skip_clumpify', action='store_true', help='跳过 Clumpify 聚类重排')
    g.add_argument('--dedup', action='store_true', help='fastp 自带去重 (默认未开启)')
    g.add_argument('--clumpify_memory', default='10g', help='clumpify 内存 (默认: 10g)')
    g.add_argument('--no_compress', action='store_true', help='输出不压缩')

    g = p.add_argument_group('计算资源')
    g.add_argument('--threads', '-t', type=int, default=20, help='线程数 (默认: 20)')
    g.add_argument('--jobs', '-j', type=int, default=2, help='并行样本数 (默认: 2)')

    g = p.add_argument_group('流程控制')
    g.add_argument('--force', action='store_true', help='强制重跑')
    g.add_argument('--keep_rrna', action='store_true', help='保留 rRNA reads')
    g.add_argument('--tmp_dir', help='临时目录')
    g.add_argument('--warn_retained', type=float, default=20.0,
                   help='report 阶段: 留存率低于该百分比标记 LOW_RETAINED (默认: 20)')

    args = p.parse_args()
    stages = set(args.stage)
    _all = 'all' in stages

    # 布局解析 (CLI > MMPV_IO_LAYOUT > legacy) 并写回环境供子进程继承
    normalize_layout_env(getattr(args, 'io_layout', None))

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
        cl = Path(args.output_dir) / layout_dir_name('d_clean') / '3.clumpify'
        fa = Path(args.output_dir) / layout_dir_name('d_clean') / '2.fasta'
        reads_dir = str(cl) if cl.exists() and any(cl.iterdir()) else str(fa)
        if not Path(reads_dir).exists():
            log.error("未找到 clean 输出, 请先 --stage clean 或指定 --input_reads")
            sys.exit(1)
        log.info("  Reads → %s (auto)", reads_dir)
    elif _all or 'clean' in stages:
        log.error("需要 --input_reads")
        sys.exit(1)
    else:
        reads_dir = None  # 仅 report 阶段无需输入

    # 执行
    if _all or 'clean' in stages:
        reads_dir = run_clean(args, log, reads_dir)
    if _all or 'deplete' in stages:
        has_db = args.host_db or (args.kraken2_db and args.host_align_db)
        has_src = args.host_fasta or (args.host_species and args.host_taxid)
        if not has_db and not has_src:
            log.error("去宿主需要宿主参考, 三选一:")
            log.error("  ① --host_db <宿主库目录>                  (hostdb/ 或 host_reference 根, 自动查找子索引)")
            log.error("  ② --kraken2_db <dir> --host_align_db <prefix>   (分开精确指定)")
            log.error("  ③ --host_fasta <FASTA> 或 --host_species X --host_taxid N   (现场构建到 {output_dir}/host_reference/)")
            sys.exit(1)
        if not has_db:
            ensure_hostdb(args, log)
        run_deplete(args, log, reads_dir)
        cleanup_00a_after_deplete(args.output_dir, log)
    # ── BBNorm (可选: 仅 --stage bbnorm 显式运行; --stage all 不包含) ──
    if 'bbnorm' in stages:
        bbnorm_out = Path(args.output_dir) / layout_dir_name('d_bbnorm')
        if not reads_dir:
            dep = Path(args.output_dir) / layout_dir_name('d_hostdep')
            reads_dir = str(dep) if dep.exists() and any(dep.iterdir()) else None
        if not reads_dir:
            log.error("bbnorm 需要 --input_reads 或先运行 deplete/clean")
            sys.exit(1)
        parts = [
            f"python {SCRIPT_DIR / 'run_bbnorm.py'}",
            f"-i {reads_dir}",
            f"-o {bbnorm_out}",
            f"-t {args.threads}",
            f"-j {args.jobs}",
        ]
        log.info("[bbnorm] 覆盖度归一化 (target=70 mindepth=2), 输入: %s", reads_dir)
        if not run_cmd(' '.join(parts), log, "BBNorm"):
            sys.exit(1)
        log.info("  Normalized reads → %s", bbnorm_out)

    if _all or 'report' in stages:
        run_report(args, log)

    log.info("=" * 50)
    log.info("预处理完成!")
    log.info("  清洗:   %s/%s", args.output_dir, layout_dir_name('d_clean'))
    log.info("  去宿主: %s/%s", args.output_dir, layout_dir_name('d_hostdep'))
    log.info("  报告:   %s/preprocessing_report.html", args.output_dir)
    log.info("=" * 50)

if __name__ == '__main__':
    main()
