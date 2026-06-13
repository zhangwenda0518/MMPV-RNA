#!/usr/bin/env python3
"""
auto_known_virus.py — 已知病毒分析主控流水线

流程:
  1. batch_virus_depth40.py     → 已知病毒快速检测
  2. batch_virus_variants.py    → 变异分析 (FreeBayes+SnpEff+SnpGenie)
  3. batch_virus_full.py        → 单倍型全长组装

用法:
  python auto_known_virus.py \\
    --reads_dir host/ --output_dir ./virus_analysis/ \\
    --ref_info ~/db/ref_info.tsv --reference ~/db/ref.fasta \\
    --tool salmon --snpeff --snpgenie \\
    --threads 40 --jobs 4 --resume
"""

import argparse, subprocess, sys, os, logging
from pathlib import Path
from datetime import datetime


def logger(out_dir):
    l = logging.getLogger("KnownVirus"); l.setLevel(logging.DEBUG); l.handlers.clear()
    os.makedirs(out_dir, exist_ok=True)
    for h in [logging.StreamHandler(), logging.FileHandler(os.path.join(out_dir, 'known_virus.log'))]:
        h.setLevel(logging.INFO)
        h.setFormatter(logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S'))
        l.addHandler(h)
    return l


def run(cmd, log, name):
    log.info("[%s] %s", name, cmd)
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=None, stderr=subprocess.PIPE)
        log.info("[%s] ✓", name)
        return True
    except subprocess.CalledProcessError as e:
        log.error("[%s] ✗ (exit=%d)", name, e.returncode)
        if e.stderr:
            tail = e.stderr.decode()[-500:] if isinstance(e.stderr, bytes) else e.stderr[-500:]
            log.error("  %s", tail)
        return False


STAGE_HELP = {
    'detect': """
  --stage detect — 已知病毒快速检测

  运行: batch_virus_depth40.py
    伪比对 (Salmon/Kallisto) 或传统比对 (Bowtie2/BWA/...) → Poisson Ratio 去假阳性

  输出: 1_FastViromeExplorer/summary/summary.tsv

  关键参数: --tool, --coverage, --ratio, --sp_thresh, --genes_cov
""",
    'variants': """
  --stage variants — 变异分析

  运行: batch_virus_variants.py
    提取 reads → 共识序列 → FreeBayes 变异检测 → SnpEff/SnpGenie 注释

  输入: 1_FastViromeExplorer/summary/summary.tsv
  输出: 2_Virus_variants_Results/{virus}/

  关键参数: --variant_caller, --snpeff, --snpgenie
""",
    'full': """
  --stage full — 单倍型全长组装

  运行: batch_virus_full.py → virus-full.py
    多工具全长组装 (SPAdes/IVA/...), 迭代 3 轮

  输入: 2_Virus_variants_Results/summary/all_summary.tsv
  输出: 3_Virus_assemblies_final/

  关键参数: --assembly_tools, --min_covered, --extra_args
""",
}


def main():
    # --stage <name> --help 处理
    for i, a in enumerate(sys.argv[1:], 1):
        if a in ('-h', '--help'):
            for j, b in enumerate(sys.argv[1:], 1):
                if b == '--stage' and j < len(sys.argv) - 1:
                    s = sys.argv[j + 1]
                    if s in STAGE_HELP:
                        print(STAGE_HELP[s])
                        sys.exit(0)
            # 默认 help
            p = argparse.ArgumentParser(description="已知病毒分析主控流水线")
            p.parse_args(['--help'])
            sys.exit(0)
    p = argparse.ArgumentParser(description="已知病毒分析主控流水线",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reads_dir", required=True, help="清洁 reads 目录 (preprocess 输出)")
    p.add_argument("--output_dir", "-o", required=True, help="输出根目录")

    # ── 参考数据库 ──
    g = p.add_argument_group("参考数据库")
    g.add_argument("--ref_info", required=True)
    g.add_argument("--reference", required=True)

    # ── Step 1: 快速检测 ──
    g = p.add_argument_group("Step 1: 快速检测 (batch_virus_depth40)")
    g.add_argument("--tool", default="bowtie2",
                   choices=["salmon","kallisto","bowtie2","bwa","minimap2","strobealign","bwa-mem2","hisat2"])
    g.add_argument("--batch_size", type=int, default=20)
    g.add_argument("--coverage", type=float, default=10.0, help="全长覆盖度下限 %%")
    g.add_argument("--ratio", type=float, default=0.3, help="泊松覆盖度比值下限")
    g.add_argument("--meandepth", type=float, default=0.0, help="最小平均深度")
    g.add_argument("--min_tpm", type=float, default=0.0)
    g.add_argument("--min_uniq_reads", type=int, default=1)
    g.add_argument("--sp_thresh", type=float, default=95.0, help="物种ANI阈值 %%")
    g.add_argument("--taxid_clusters", help="同义 TaxID 合并映射文件")
    g.add_argument("--genes_cov", help="转录覆盖率文件 (双轨过滤)")
    g.add_argument("--min_gene_total_cov", type=float, default=80.0)
    g.add_argument("--min_gene_avr_cov", type=float, default=5.0)
    g.add_argument("--use_coverm", action="store_true", help="启用 CoverM 严格清洗 (仅传统比对)")
    g.add_argument("--min_aln_len", type=int, default=80)
    g.add_argument("--min_aln_prop", type=float, default=0.85)
    g.add_argument("--min_pid", type=float, default=0.90)
    g.add_argument("--single_end", action="store_true", help="强制单端模式")
    g.add_argument("--keep_tmp", action="store_true", help="保留中间 BAM 文件")
    g.add_argument("--verbose", action="store_true", help="详细底层日志")

    # ── Step 2: 变异分析 ──
    g = p.add_argument_group("Step 2: 变异分析 (batch_virus_variants)")
    g.add_argument("--variant_caller", default="freebayes", choices=["freebayes","ivar","lofreq"])
    g.add_argument("--snpeff", action="store_true")
    g.add_argument("--snpeff_jar", default="~/biosoft/snpEff/snpEff.jar")
    g.add_argument("--snpeff_config", default="~/biosoft/snpEff/snpEff.config")
    g.add_argument("--snpeff_mem", default="4g")
    g.add_argument("--snpgenie", action="store_true")
    g.add_argument("--no_extract_reads", action="store_true", help="不提取 reads")
    g.add_argument("--no_consensus", action="store_true", help="不生成共识序列")
    g.add_argument("--no_call_variants", action="store_true", help="不检出变异")
    g.add_argument("--bam", help="已有 BAM 文件夹 (替代 --reads_dir)")
    g.add_argument("--disable_dynamic_vcf", action="store_true", help="禁用动态 VCF")
    g.add_argument("--vc_qual", type=int, default=20, help="共识序列最低质量 (默认20)")
    g.add_argument("--vc_depth", type=int, default=10, help="共识序列最低深度 (默认10)")
    g.add_argument("--vc_freq", type=float, default=0.5, help="共识序列最低频率 (默认0.5)")
    g.add_argument("--vc_ambig", type=float, default=0.25, help="共识序列 IUPAC 模糊阈值 (默认0.25)")

    # ── Step 3: 全长组装 ──
    g = p.add_argument_group("Step 3: 全长组装 (batch_virus_full)")
    g.add_argument("--assembly_tools", default="all")
    g.add_argument("--min_covered", type=float, default=10.0)
    g.add_argument("--extra_args", default="--iter 3 --vc-min-depth 1", help="virus-full 额外参数")
    g.add_argument("--virus_full_script", default=None, help="virus-full.py 路径")

    # ── 流程控制 ──
    g = p.add_argument_group("流程控制")
    g.add_argument("--stage", default="all", choices=["all","detect","variants","full"],
                   help="运行阶段 (all/detect/variants/full). --help-stage detect 查看详情")
    g.add_argument("--resume", action="store_true", help="断点续传 (Step1+2)")
    g.add_argument("--threads", type=int, default=40)
    g.add_argument("--jobs", type=int, default=4)
    g.add_argument("--align_threads", type=int, default=8, help="单样本比对线程")

    args = p.parse_args()
    script_dir = Path(__file__).parent.resolve()
    out = Path(args.output_dir).resolve()
    reads = Path(args.reads_dir).resolve()

    if not reads.exists():
        sys.exit(f"ERROR: reads 目录不存在: {reads}")

    log = logger(str(out))
    log.info("=" * 50)
    log.info("已知病毒分析 | Stage=%s | Threads=%d Jobs=%d", args.stage, args.threads, args.jobs)
    log.info("  Reads:  %s", reads)
    log.info("  Output: %s", out)
    log.info("=" * 50)

    # ═══ Step 1: 快速检测 ═══
    detect_dir = out / "1_FastViromeExplorer"
    if args.stage in ("all", "detect"):
        log.info("=" * 50)
        log.info("[1] 已知病毒快速检测")
        parts = [
            f"python {script_dir / 'batch_virus_depth40.py'}",
            f"--input_dir {reads}",
            f"--output_dir {detect_dir}",
            f"--ref_info {args.ref_info}",
            f"--reference {args.reference}",
            f"--tool {args.tool}",
            f"--threads {args.threads}",
            f"--align_threads {args.align_threads}",
            f"--batch_size {args.batch_size}",
            f"--coverage {args.coverage}",
            f"--ratio {args.ratio}",
            f"--meandepth {args.meandepth}",
            f"--min_tpm {args.min_tpm}",
            f"--min_uniq_reads {args.min_uniq_reads}",
            f"--sp_thresh {args.sp_thresh}",
        ]
        if args.genes_cov:
            parts += [
                f"--genes_cov {args.genes_cov}",
                f"--min_gene_total_cov {args.min_gene_total_cov}",
                f"--min_gene_avr_cov {args.min_gene_avr_cov}",
            ]
        if args.taxid_clusters:
            parts.append(f"--taxid_clusters {args.taxid_clusters}")
        if args.use_coverm:
            parts.append("--use_coverm")
            parts += [f"--min_aln_len {args.min_aln_len}", f"--min_aln_prop {args.min_aln_prop}", f"--min_pid {args.min_pid}"]
        if args.single_end:
            parts.append("--single_end")
        if args.keep_tmp:
            parts.append("--keep_tmp")
        if args.verbose:
            parts.append("--verbose")
        if args.resume:
            parts.append("--resume")
        if not run(' '.join(parts), log, "batch_virus_depth40"):
            sys.exit(1)
        log.info("  检测完成 → %s", detect_dir)

    # ═══ Step 2: 变异分析 ═══
    summary_in = detect_dir / "summary" / "summary.tsv"
    variants_dir = out / "2_Virus_variants_Results"
    if args.stage in ("all", "variants"):
        log.info("=" * 50)
        log.info("[2] 变异分析")
        if not summary_in.exists():
            log.error("summary 不存在: %s (请先 --stage detect)", summary_in)
            sys.exit(1)
        parts = [
            f"python {script_dir / 'batch_virus_variants.py'}",
            f"--summary {summary_in}",
            f"--info {args.ref_info}",
            f"--reference {args.reference}",
            f"--variant_caller {args.variant_caller}",
            f"--output_dir {variants_dir}",
            f"--threads {args.threads}",
            f"--jobs {args.jobs}",
        ]
        if not args.no_extract_reads:
            parts.append("--extract_reads")
        if not args.no_consensus:
            parts.append("--consensus")
        if not args.no_call_variants:
            parts.append("--call_variants")
        if args.snpeff:
            parts += ["--snpeff", f"--snpeff_jar {args.snpeff_jar}",
                      f"--snpeff_config {args.snpeff_config}", f"--snpeff_mem {args.snpeff_mem}"]
        if args.snpgenie:
            parts.append("--snpgenie")
        if args.bam:
            parts += ["--bam", args.bam]
        else:
            parts += ["--fastq", str(reads)]
        if args.disable_dynamic_vcf:
            parts.append("--disable_dynamic_vcf")
        parts += [f"-q {args.vc_qual}", f"-d {args.vc_depth}", f"-f {args.vc_freq}", f"-a {args.vc_ambig}"]
        if args.resume:
            parts.append("--resume")
        if not run(' '.join(parts), log, "batch_virus_variants"):
            log.warning("  变异分析部分失败, 检查日志")
        else:
            log.info("  变异分析完成 → %s", variants_dir)

    # ═══ Step 3: 全长组装 ═══
    full_dir = out / "3_Virus_assemblies_final"
    if args.stage in ("all", "full"):
        log.info("=" * 50)
        log.info("[3] 单倍型全长组装")
        var_summary = variants_dir / "summary" / "all_summary.tsv"
        if not var_summary.exists():
            log.error("变异 summary 不存在: %s (请先 --stage variants)", var_summary)
            sys.exit(1)
        vsi = args.virus_full_script or str(script_dir.parent / "virus-full.py")
        parts = [
            f"python {script_dir / 'batch_virus_full.py'}",
            f"--downstream_dir {variants_dir}",
            f"--summary {var_summary}",
            f"--clean_data {reads}",
            f"--virus_full_script {vsi}",
            f"--outdir {full_dir}",
            f"--assembly_tools {args.assembly_tools}",
            f"--jobs {args.jobs}",
            f"--threads {args.threads}",
            f"--extra_args '{args.extra_args}'",
            f"--min_covered {args.min_covered}",
        ]
        if args.resume:
            parts.append("--resume")
        if not run(' '.join(parts), log, "batch_virus_full"):
            log.warning("  全长组装部分失败, 检查日志")
        else:
            log.info("  全长组装完成 → %s", full_dir)

    log.info("=" * 50)
    log.info("已知病毒分析完成! | %s", datetime.now().strftime('%H:%M:%S'))
    log.info("=" * 50)


if __name__ == "__main__":
    main()
