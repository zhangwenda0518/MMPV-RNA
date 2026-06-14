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


def logger(out_dir, level='INFO'):
    l = logging.getLogger("KnownVirus"); l.setLevel(logging.DEBUG); l.handlers.clear()
    os.makedirs(out_dir, exist_ok=True)
    console_level = getattr(logging, level.upper(), logging.INFO)
    for h in [logging.StreamHandler(), logging.FileHandler(os.path.join(out_dir, 'known_virus.log'))]:
        h.setLevel(console_level)
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
  --stage detect — 已知病毒快速检测 (伪比对/传统比对 + Poisson打假 + 双轨过滤)

  子脚本: batch_virus_depth40.py  (支持 Parquet 批次断点续传)
    双引擎       : Salmon/Kallisto 极速伪比对, Bowtie2/BWA/StrobeAlign 传统比对
    Poisson打假  : Poisson_Ratio 建模, 精准剔除局部堆叠假阳性
    双轨过滤     : A轨(全基因组≥coverage+泊松≥ratio) + B轨(基因转录覆盖)
                   RNA病毒→A+B双轨, DNA病毒→仅A轨
    输出指标     : EM_Reads, CPM, FPKM, TPM, Avg_Read_ANI, Poisson_Ratio, Pi

  编排器透传参数 (★ 必需 | · 可选):
    ★ --reads_dir         清洁 reads 目录 (preprocess/去宿主后)
    ★ --output_dir        输出根目录 (子目录: 1_FastViromeExplorer/)
    ★ --ref_info          病毒参考信息 TSV (Accession/Taxid/Species/Segment 列)
    ★ --reference         参考基因组 FASTA (建索引自动复用)
    · --tool              {salmon,kallisto,bowtie2,bwa,bwa-mem2,hisat2,minimap2,strobealign}
                           (默认 bowtie2)
    · --batch_size        批次刷盘保护数量 (默认 20)
    · --coverage          全长覆盖度下限 %% (默认 10.0)
    · --ratio             泊松覆盖度比值下限 (默认 0.3)
    · --meandepth         最小平均深度 (默认 0.0)
    · --min_tpm           最小 TPM 值 (默认 0.0)
    · --min_uniq_reads    最少独特比对 reads 数 (默认 1)
    · --sp_thresh         物种 ANI 阈值 %% (仅传统比对, 默认 95.0)
    · --taxid_clusters    同义 TaxID 合并映射文件
    · --genes_cov         转录覆盖率文件 (启动双轨B轨, RNA病毒全基因组+基因区过滤)
    · --min_gene_total_cov 最低转录区总覆盖 %% (默认 80.0)
    · --min_gene_avr_cov  最低转录区平均覆盖 (默认 5.0)
    · --use_coverm        启用 CoverM 严格清洗 (仅传统比对)
    · --min_aln_len       CoverM 最小比对长度 (默认 80)
    · --min_aln_prop      CoverM 最小比对比例 (默认 0.85)
    · --min_pid           CoverM 最小序列相似度 (默认 0.90)
    · --single_end        强制单端模式
    · --keep_tmp          保留中间 BAM 文件 (调试用)
    · --verbose           输出详细底层日志
    · --resume            批次断点续传 (跳过已完成 .parquet)
    · -t/--threads        并发进程数 (默认 40)
    · --align_threads     单样本内部比对线程 (默认 8)

  输出: 1_FastViromeExplorer/summary/summary.tsv (含 best.summary.tsv)
""",
    'variants': """
  --stage variants — 变异分析 (reads提取 → 共识 → 变异检出 → SnpEff → SnpGenie)

  子脚本: batch_virus_variants.py  (支持任务级断点续传)
    三池解耦       : 变异检出池 / SnpEff注释池 / SNPGenie进化池 独立并行
    ─────────────────────────────────────────────────
    [1] 提取 reads  : 从 BAM 提取目标病毒 reads (或复用已有 --bam)
    [2] 共识序列    : iVar/viral_consensus 生成共识 (--vc_qual/depth/freq 质量控制)
    [3] 变异检出    : FreeBayes/LoFreq/iVar → 动态 VCF 过滤 (QUAL/DP/AF自适应)
    [4] SnpEff 注释 : 自动 NCBI 下载 GenBank 构建本地 DB → 变异功能注释
    [5] SnpGenie     : dN/dS 选择压力分析 (种群遗传学)

  编排器透传参数 (★ 必需 | · 可选):
    ★ --reads_dir         清洁 reads 目录
    ★ --output_dir        输出根目录 (子目录: 2_Virus_variants_Results/)
    ★ --ref_info          病毒参考信息 TSV
    ★ --reference         参考基因组 FASTA
    · --variant_caller    {freebayes,ivar,lofreq} (默认 freebayes)
    · --snpeff            启用 SnpEff 注释
    · --snpeff_jar        snpEff.jar 路径 (默认 ~/biosoft/snpEff/snpEff.jar)
    · --snpeff_config     snpEff 配置文件 (默认 ~/biosoft/snpEff/snpEff.config)
    · --snpeff_mem        SnpEff 内存限制 (默认 4g)
    · --snpgenie          启用 SnpGenie dN/dS 分析
    · --no_extract_reads  跳过 reads 提取 (已有提取结果)
    · --no_consensus      跳过共识序列生成
    · --no_call_variants  跳过变异检出 (仅做注释/进化分析)
    · --bam               已有 BAM 目录 (替代 --reads_dir re-extraction)
    · --disable_dynamic_vcf 禁用动态 VCF 质量过滤 (使用固定阈值)
    · -q/--vc_qual        共识碱基最低质量 Phred (默认 20)
    · -d/--vc_depth       共识最低深度 (默认 5)
    · -f/--vc_freq        变异最低频率阈值 (默认 0.5)
    · -a/--vc_ambig       低覆盖碱基填充字符 (默认 N)
    · -t/--threads        单任务线程数 (默认 40)
    · -j/--jobs           并行任务数 (默认 4)
    · --resume            断点续传 (跳过已完成任务)

  输出:
    2_Virus_variants_Results/
    ├── virus-variants/       (VCF + 变异频率 TSV)
    ├── virus-SnpEff/         (SnpEff 注释 VCF + 摘要 TSV)
    ├── virus-SNPGenie/       (SNPGenie dN/dS 结果)
    ├── virus-consensus/      (共识序列 FASTA)
    ├── virus_reads/          (靶向提取 reads)
    └── summary/all_summary.tsv
""",
    'full': """
  --stage full — 单倍型全长组装 (批量调度 virus-full.py)

  子脚本: batch_virus_full.py  (支持断点续传 + 目录归档)
    流程      : 靶向提取 reads + 原始 Clean reads → virus-full.py 多工具组装
    目录结构  : {Taxonomy}_{Accession}/{Sample}_{Accession}/  按分类归档
    组装工具  : SPAdes/IVA/... (via --assembly_tools)

  编排器透传参数 (★ 必需 | · 可选):
    ★ --reads_dir         清洁 reads 目录 (原始 reads)
    ★ --output_dir        输出根目录 (子目录: 3_Virus_assemblies_final/)
    · --assembly_tools    启用的组装工具 (默认 all)
    · --min_covered       跳过 Covered% 小于等于该值的记录 (默认 10.0)
    · --extra_args        传递给 virus-full.py 的额外参数 (默认: --iter 3 --vc-min-depth 1)
    · --virus_full_script  virus-full.py 路径 (默认: 上级目录)
    · --gb                全局 GenBank (.gb) 文件 (病毒注释信息)
    · -t/--threads        单任务线程数 (默认 40)
    · -j/--jobs           并行任务数 (默认 4)
    · --resume            断点续传 (检测 final.fasta 存在则跳过)

  输出:
    3_Virus_assemblies_final/{Taxonomy}_{Accession}/{Sample}_{Accession}/
    ├── final.fasta        (最终组装)
    ├── scaffolds.fasta    (scaffolds)
    └── contigs.fasta      (contigs)
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
    g.add_argument("--vc_depth", type=int, default=5, help="共识序列最低深度 (默认5)")
    g.add_argument("--vc_freq", type=float, default=0.5, help="共识序列最低频率 (默认0.5)")
    g.add_argument("--vc_ambig", type=str, default="N", help="低于阈值的碱基用此字符代替 (默认: N)")

    # ── Step 3: 全长组装 ──
    g = p.add_argument_group("Step 3: 全长组装 (batch_virus_full)")
    g.add_argument("--assembly_tools", default="all")
    g.add_argument("--min_covered", type=float, default=10.0)
    g.add_argument("--extra_args", default="--iter 3 --vc-min-depth 1", help="virus-full 额外参数")
    g.add_argument("--virus_full_script", default=None, help="virus-full.py 路径")
    g.add_argument("--gb", help="全局 GenBank (.gb) 文件 (virus-full.py 注释用)")

    # ── 流程控制 ──
    g = p.add_argument_group("流程控制")
    g.add_argument("--stage", default="all", choices=["all","detect","variants","full"],
                   help="运行阶段 (all/detect/variants/full). --help-stage detect 查看详情")
    g.add_argument("--resume", action="store_true", help="断点续传 (Step1+2)")
    g.add_argument("--force", action="store_true", help="强制重跑 (忽略断点, 覆盖已有结果)")
    g.add_argument("--dry-run", action="store_true", help="仅显示配置和任务概览, 不实际执行")
    g.add_argument("--log-level", default="INFO", choices=["DEBUG","INFO","WARNING","ERROR"],
                   help="日志级别 (默认: INFO)")
    g.add_argument("--threads", type=int, default=40)
    g.add_argument("--jobs", type=int, default=4)
    g.add_argument("--align_threads", type=int, default=8, help="单样本比对线程")

    args = p.parse_args()
    script_dir = Path(__file__).parent.resolve()
    out = Path(args.output_dir).resolve()
    reads = Path(args.reads_dir).resolve()

    if not reads.exists():
        sys.exit(f"ERROR: reads 目录不存在: {reads}")

    log = logger(str(out), level=args.log_level)
    log.info("=" * 50)
    log.info("已知病毒分析 | Stage=%s | Threads=%d Jobs=%d", args.stage, args.threads, args.jobs)
    log.info("  Reads:  %s", reads)
    log.info("  Output: %s", out)
    if args.force:
        log.info("  模式:    强制重跑 (--force)")
    elif args.resume:
        log.info("  模式:    断点续传 (--resume)")
    log.info("  Log:     %s", args.log_level)
    log.info("=" * 50)

    # ── --dry-run: 仅显示配置 ──
    if args.dry_run:
        log.info("")
        log.info("═══ DRY-RUN 模式 — 不执行任何计算 ═══")
        log.info("  阶段:    %s", args.stage)
        log.info("  Reads:   %s", reads)
        log.info("  Output:  %s", out)
        log.info("  检测工具: %s", args.tool)
        if args.stage in ("all", "detect"):
            log.info("  [1/3] 快速检测:  batch_virus_depth40.py")
        if args.stage in ("all", "variants"):
            log.info("  [2/3] 变异分析:  batch_virus_variants.py")
            log.info("         caller=%s snpeff=%s snpgenie=%s", args.variant_caller, args.snpeff, args.snpgenie)
        if args.stage in ("all", "full"):
            log.info("  [3/3] 全长组装:  batch_virus_full.py")
        log.info("═══ DRY-RUN 结束 ═══")
        return

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
        if args.resume and not args.force:
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
        if args.resume and not args.force:
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
        if args.resume and not args.force:
            parts.append("--resume")
        if args.gb:
            parts.append(f"--gb {args.gb}")
        if not run(' '.join(parts), log, "batch_virus_full"):
            log.warning("  全长组装部分失败, 检查日志")
        else:
            log.info("  全长组装完成 → %s", full_dir)

    log.info("=" * 50)
    log.info("已知病毒分析完成! | %s", datetime.now().strftime('%H:%M:%S'))
    log.info("=" * 50)


if __name__ == "__main__":
    main()
