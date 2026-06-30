#!/usr/bin/env python3
"""
auto_known_virus.py - Known Virus Analysis Pipeline
====================================================
7-stage automated pipeline for known virus detection, variant analysis,
full-length assembly, and post-hoc characterization.

Stages:
  1. detect    - Rapid quantification (batch_virus_depth.py)
  2. filter    - High-confidence filtering (filter_summary.py) [auto with --filter]
  3. variants  - Variant calling + SnpEff + SnpGenie (batch_virus_variants.py)
  4. full      - De novo full-length assembly (virus-full.py)
  5. extract   - Extract longest contigs (extract_full_fasta.py)
  6. post      - VCF visualization + SnpEff macro + MAF + SnpGenie
  7. capheine   - Positive selection analysis (capheine_pipeline.py)
  8. similarity - Full-length similarity panorama (virus_auto_pipeline.py)
  9. dvg        - DVG & recombination detection (batch_virema_dvg.py)
 10. report     - Generate summary report + AI interpretation prompts

Output structure:
  output_dir/
    1_FastViromeExplorer/     Stage 1: detection results
    2_Virus_variants_Results/ Stage 3: variant analysis
    3_Virus_assemblies_final/ Stage 4: full assemblies
    4_assemblies_clean/       Stage 5: extracted contigs
    5_post_analysis/          Stage 6: post-hoc viz
    6_capheine/               Stage 7: selection analysis
    7_similarity/             Stage 8: similarity panorama
    logs/                     Pipeline logs
"""

import argparse
import subprocess
import sys
import os
import logging
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

SCRIPT_DIR = Path(__file__).resolve().parent

def setup_logger(out_dir, level="INFO"):
    logger = logging.getLogger("KnownVirus")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    os.makedirs(out_dir, exist_ok=True)
    console_level = getattr(logging, level.upper(), logging.INFO)
    for handler in [
        logging.StreamHandler(),
        logging.FileHandler(os.path.join(out_dir, "known_virus.log")),
    ]:
        handler.setLevel(console_level)
        handler.setFormatter(
            logging.Formatter("[%(asctime)s] %(message)s", datefmt="%H:%M:%S")
        )
        logger.addHandler(handler)
    return logger


def run(cmd, log, name):
    """Execute a shell command with error handling.

    stdout 和 stderr 均继承父进程（由 shell 重定向到日志文件），
    确保 tqdm 进度条等 stderr 输出可见。
    """
    log.info("[%s] %s", name, cmd)
    try:
        subprocess.run(
            cmd, shell=True, check=True, stdout=None, stderr=None
        )
        log.info("[%s] OK", name)
        return True
    except subprocess.CalledProcessError as e:
        log.error("[%s] FAILED (exit=%d)", name, e.returncode)
        return False


def find_virus_dir(base_dir, sub_dir, acc):
    """Match a virus directory by accession substring."""
    d = base_dir / sub_dir
    if not d.exists():
        return None
    if (d / acc).exists():
        return d / acc
    for child in d.iterdir():
        if child.is_dir() and acc in child.name:
            return child
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Known Virus Analysis Pipeline (7-stage)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ---- I/O ----
    g = parser.add_argument_group("Input/Output")
    g.add_argument("--reads_dir", default=None, help="Clean reads directory")
    g.add_argument("--output_dir", "-o", default=None, help="Output root directory")
    g.add_argument("--ref_info", default=None, help="Reference info TSV")
    g.add_argument("--reference", default=None, help="Reference genome FASTA")

    # ---- Stage control ----
    g = parser.add_argument_group("Stage Control")
    g.add_argument(
        "--stage",
        default="all",
        choices=["all", "detect", "filter", "variants", "full", "extract", "post", "capheine", "similarity", "dvg", "report"],
        help="Which stage to run (default: all)",
    )
    g.add_argument("--no-resume", action="store_true", help="Disable checkpoint resume (always re-run)")
    g.add_argument("--force", action="store_true", help="Force re-run, ignore all checkpoints")
    g.add_argument("--dry_run", action="store_true", help="Preview only, no execution")
    g.add_argument("--log_level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    # ---- Concurrency ----
    g = parser.add_argument_group("Concurrency")
    g.add_argument("--threads", type=int, default=60)
    g.add_argument("--jobs", type=int, default=4)
    g.add_argument("--align_threads", type=int, default=8)

    # ---- Stage 1: Detect ----
    g = parser.add_argument_group("Stage 1: Detection (batch_virus_depth)")
    g.add_argument("--tool", default="salmon", choices=["salmon", "kallisto", "bowtie2", "bwa", "minimap2", "strobealign", "bwa-mem2", "hisat2"])
    g.add_argument("--batch_size", type=int, default=20)
    g.add_argument("--coverage", type=float, default=10.0, help="Min coverage %%")
    g.add_argument("--ratio", type=float, default=0.3, help="Min Poisson ratio")
    g.add_argument("--meandepth", type=float, default=0.5, help="Min mean depth (default 0.5X)")
    g.add_argument("--min_tpm", type=float, default=1.0, help="Min TPM (default 1.0)")
    g.add_argument("--min_uniq_reads", type=int, default=10, help="Min unique reads (default 10)")
    g.add_argument("--sp_thresh", type=float, default=95.0, help="Species ANI threshold %%")
    g.add_argument("--genes_cov", help="Gene coverage file (dual-track filter)")
    g.add_argument("--min_gene_total_cov", type=float, default=80.0)
    g.add_argument("--min_gene_avr_cov", type=float, default=5.0)
    g.add_argument("--taxid_clusters", help="Synonymous TaxID mapping file")
    g.add_argument("--use_coverm", action="store_true", help="Enable CoverM (traditional only)")
    g.add_argument("--min_aln_len", type=int, default=80)
    g.add_argument("--min_aln_prop", type=float, default=0.85)
    g.add_argument("--min_pid", type=float, default=0.90)
    g.add_argument("--single_end", action="store_true", help="Force single-end mode")
    g.add_argument("--keep_tmp", action="store_true", help="Keep intermediate BAM files")
    g.add_argument("--verbose", action="store_true", help="Verbose logging")

    # ---- Stage 2: Filter ----
    g = parser.add_argument_group("Stage 2: Filter (filter_summary)")
    g.add_argument("--filter", action="store_true", help="Enable high-confidence filtering")
    g.add_argument("--filter_cov", type=float, default=50.0, help="Min coverage %% for filter")
    g.add_argument("--filter_depth", type=float, default=5.0, help="Min depth for filter")
    g.add_argument("--filter_reads", type=float, default=100.0, help="Min reads for filter")
    g.add_argument("--filter_keyword", type=str, help="Keyword filter (e.g. Cytorhabdovirus)")
    g.add_argument("--filter_tpm", type=float, default=0.0, help="Min TPM for filter")
    g.add_argument("--filter_poisson", type=float, default=0.0, help="Min Poisson_Ratio for filter")

    # ---- Stage 3: Variants ----
    g = parser.add_argument_group("Stage 3: Variants (batch_virus_variants)")
    g.add_argument("--variant_caller", default="ivar", choices=["freebayes", "ivar", "lofreq"])
    g.add_argument("--snpeff", action="store_true", help="Enable SnpEff annotation")
    g.add_argument("--snpeff_jar", default=str(SCRIPT_DIR / "../biosoft/snpEff/snpEff.jar"))
    g.add_argument("--snpeff_config", default=str(SCRIPT_DIR / "../biosoft/snpEff/snpEff.config"))
    g.add_argument("--snpeff_mem", default="4g")
    g.add_argument("--snpgenie", action="store_true", help="Enable SnpGenie analysis")
    g.add_argument("--no_extract_reads", action="store_true")
    g.add_argument("--no_consensus", action="store_true")
    g.add_argument("--no_call_variants", action="store_true")
    g.add_argument("--bam", help="Existing BAM directory (alternative to --reads_dir)")
    g.add_argument("--disable_dynamic_vcf", action="store_true")
    g.add_argument("--vc_qual", type=int, default=20, help="Consensus min quality")
    g.add_argument("--vc_depth", type=int, default=5, help="Consensus min depth")
    g.add_argument("--vc_freq", type=float, default=0.5, help="Consensus min frequency")
    g.add_argument("--vc_ambig", type=str, default="N", help="Low-coverage base fill char")

    # ---- Stage 4: Full Assembly ----
    g = parser.add_argument_group("Stage 4: Full Assembly (virus-full)")
    g.add_argument("--assembly_tools", default="all")
    g.add_argument("--min_covered", type=float, default=10.0)
    g.add_argument("--extra_args", default="--iter 3 --vc-min-depth 1")
    g.add_argument("--virus_full_script", default=None, help="Path to virus-full.py")
    g.add_argument("--gb", help="GenBank file for annotation")

    # ---- Stage 5: Extract ----
    g = parser.add_argument_group("Stage 5: Extract Assemblies")
    g.add_argument("--extract_target", default="11.Ultimate_Circular_Result.fasta")
    g.add_argument("--max_n_genome", type=float, default=5.0, help="N content threshold%% for extract N-fill (default: 5)")
    g.add_argument("--min_length", type=int, default=150, help="Min contig length for extract (default: 150)")

    # ---- Stage 6: Post-hoc ----
    g = parser.add_argument_group("Stage 6: Post-hoc Visualization")
    g.add_argument("--post_min_dp", type=int, default=50, help="VCF min depth for post-hoc")
    g.add_argument("--post_min_af", type=float, default=0.05, help="VCF min allele freq")
    g.add_argument("--skip_vcf_viz", action="store_true")
    g.add_argument("--skip_vcf_merge", action="store_true", help="Skip VCF merge + PCA + distance matrix")
    g.add_argument("--skip_snpeff_macro", action="store_true")
    g.add_argument("--skip_maftools", action="store_true")
    g.add_argument("--skip_snpgenie", action="store_true")

    # ---- Stage 7: Capheine ----
    g = parser.add_argument_group("Stage 7: Capheine (Positive Selection)")
    g.add_argument("--capheine_ref", help="Reference CDS FASTA for capheine")
    g.add_argument("--capheine_unaligned", help="Unaligned sequences FASTA")
    g.add_argument("--capheine_fg", help="Foreground taxa list")
    g.add_argument("--capheine_code", default="1", help="Genetic code (default: 1=Universal)")

    # ---- Stage 8: Similarity ----
    g = parser.add_argument_group("Stage 8: Similarity Panorama (virus_auto_pipeline)")
    g.add_argument("--sim_ref", help="GenBank accession or .gb file for similarity analysis")
    g.add_argument("--sim_mode", default="filter", choices=["strict", "filter", "fill", "all"])
    g.add_argument("--sim_cdhit", action="store_true", help="Enable CD-HIT dedup")

    # ---- Stage 9: DVG ----
    g = parser.add_argument_group("Stage 9: DVG & Recombination (batch_virema_dvg)")
    g.add_argument("--virema_script", default=str(SCRIPT_DIR / "../biosoft/virema/ViReMa.py"), help="Path to ViReMa.py")
    g.add_argument("--dvg_seed", type=int, default=25, help="ViReMa seed length (default: 25)")
    g.add_argument("--dvg_mindel", type=int, default=15, help="Microdeletion threshold (default: 15)")
    g.add_argument("--dvg_min_cov", type=float, default=80.0, help="Min coverage%% for DVG analysis (default: 80)")
    g.add_argument("--dvg_shm", action="store_true", help="Use /dev/shm RAM disk for ViReMa")
    g.add_argument("--dvg_reads", default=None, help="FASTQ reads dir for DVG (default: same as --reads_dir)")

    # ---- Stage 10: Report ----
    g = parser.add_argument_group("Stage 10: Report Generation (generate_pipeline_report)")
    g.add_argument("--report_ai", action="store_true", help="Include AI interpretation prompts")
    g.add_argument("--ai_api_key", default=None, help="DeepSeek/OpenAI API key for AI interpretation")
    g.add_argument("--ai_model", default="deepseek-chat", help="LLM model")

    # ---- Profile support ----
    g = parser.add_argument_group("Profile")
    g.add_argument("--profile", default=None, help="YAML profile with default parameters")

    # Pre-parse only --profile from raw argv
    profile_file = None
    for i, a in enumerate(sys.argv[1:], 1):
        if a == '--profile' and i < len(sys.argv):
            profile_file = sys.argv[i + 1]
        elif a.startswith('--profile='):
            profile_file = a.split('=', 1)[1]

    if profile_file:
        try:
            import yaml
            with open(profile_file, 'r') as f:
                profile = yaml.safe_load(f)
            parser.set_defaults(**{k: v for k, v in profile.items() if v is not None})
        except ImportError:
            print(f"[WARNING] pyyaml 未安装, 跳过 profile 加载: {profile_file}", file=sys.stderr)
        except Exception as e:
            print(f"[WARNING] profile 加载失败 ({profile_file}): {e}", file=sys.stderr)

    args = parser.parse_args()

    # Manual validation (required args can come from profile)
    missing = []
    for param in ['reads_dir', 'output_dir', 'ref_info', 'reference']:
        if getattr(args, param, None) is None:
            missing.append(f'--{param}')
    if missing:
        parser.error(f"missing required arguments: {', '.join(missing)} (use --profile or pass explicitly)")

    # ---- Setup ----
    script_dir = Path(__file__).parent.resolve()
    out = Path(args.output_dir).resolve()

    # Redirect all temp files to pipeline's own tmp dir (avoid /tmp overflow)
    _pipeline_tmp = out / "tmp"
    _pipeline_tmp.mkdir(parents=True, exist_ok=True)
    os.environ["TMPDIR"] = str(_pipeline_tmp)
    os.environ["TMP"] = str(_pipeline_tmp)
    os.environ["TEMP"] = str(_pipeline_tmp)
    reads = Path(args.reads_dir).resolve()

    if not reads.exists():
        sys.exit(f"ERROR: reads directory not found: {reads}")

    log = setup_logger(str(out), level=args.log_level)
    log.info("=" * 55)
    log.info("Known Virus Pipeline | Stage=%s | Threads=%d Jobs=%d", args.stage, args.threads, args.jobs)
    log.info("  Reads:  %s", reads)
    log.info("  Output: %s", out)
    if args.force:
        log.info("  Mode:   FORCE (full re-run)")
    elif args.no_resume:
        log.info("  Mode:   NO-RESUME (always re-run)")
    else:
        log.info("  Mode:   RESUME (skip completed, default)")
    log.info("=" * 55)

    # ---- Dry-run ----
    if args.dry_run:
        s = args.stage
        log.info("")
        log.info("=== DRY-RUN ===")
        log.info("  Tool:    %s", args.tool)
        if s in ("all", "detect"):   log.info("  [1/10] Detect:   batch_virus_depth.py")
        if args.filter or s == "filter": log.info("  [2/10] Filter:   filter_summary.py")
        if s in ("all", "variants"):  log.info("  [3/10] Variants: batch_virus_variants.py (caller=%s snpeff=%s snpgenie=%s)", args.variant_caller, args.snpeff, args.snpgenie)
        if s in ("all", "full"):      log.info("  [4/10] Full:     batch_virus_full.py")
        if s in ("all", "extract"):   log.info("  [5/10] Extract:  extract_full_fasta.py")
        if s in ("all", "post"):      log.info("  [6/10] Post-hoc: VCF viz + SnpEff + MAF + SnpGenie")
        if s in ("all", "capheine"):  log.info("  [7/10] Capheine: positive selection analysis")
        if s in ("all", "similarity"): log.info("  [8/10] Similarity: virus_auto_pipeline.py")
        if s in ("all", "dvg"):       log.info("  [9/10] DVG:    batch_virema_dvg.py")
        if s in ("all", "report"):    log.info("  [10/10] Report: generate_pipeline_report.py")
        log.info("=== DRY-RUN END ===")
        return

    # ---- Shared paths ----
    detect_dir = out / "1_FastViromeExplorer"
    filter_dir = out / "2_Virus_result_filter"
    variants_dir = out / "3_Virus_variants_Results"
    full_dir = out / "4_Virus_assemblies_final"
    extract_dir = out / "5_assemblies_clean"
    post_dir = out / "6_post_analysis"
    capheine_dir = out / "7_capheine"
    similarity_dir = out / "8_similarity"
    dvg_dir = out / "9_virema_dvg"
    report_dir = out / "10_Reports"

    best_summary = detect_dir / "summary" / "all_viruses.best.summary.tsv"
    high_conf = filter_dir / "high_conf.summary.tsv"

    def get_summary():
        """动态获取当前最优 summary（filter 运行后自动切换到 high_conf）"""
        return high_conf if high_conf.exists() else best_summary

    def add_stage_log(stage_dir, stage_name):
        """Attach a per-stage FileHandler so logs go to both console and stage dir."""
        stage_dir.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(str(stage_dir / f"run_{stage_name}.log"), mode='a')
        fh.setLevel(logging.INFO)
        fh.setFormatter(logging.Formatter("[%(asctime)s] %(message)s", datefmt="%H:%M:%S"))
        log.addHandler(fh)
        return fh

    def remove_stage_log(handler):
        """Remove stage-specific handler after completion."""
        if handler:
            log.removeHandler(handler)
            handler.close()

    # ═══════════════════════════════════════════
    # Stage 1: Detection
    # ═══════════════════════════════════════════
    if args.stage in ("all", "detect"):
        _sh = add_stage_log(detect_dir, "1_detect")
        if not args.force and best_summary.exists():
            log.info("[1/10] Detection: checkpoint OK, skip")
        else:
            log.info("-" * 40)
            log.info("[1/10] Rapid Virus Detection")
            # checkpoint check
            parts = [
            f"python {script_dir / 'batch_virus_depth.py'}",
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
        parts += [
            f"--min_aln_len {args.min_aln_len}",
            f"--min_aln_prop {args.min_aln_prop}",
            f"--min_pid {args.min_pid}",
            ]
        if args.single_end:
            parts.append("--single_end")
        if args.keep_tmp:
            parts.append("--keep_tmp")
        if args.verbose:
            parts.append("--verbose")
        if not args.no_resume and not args.force:
            parts.append("--resume")
        if not run(" ".join(parts), log, "batch_virus_depth"):
            sys.exit(1)
        # ── Stage 1 viz: batch_plot_virus_depth.py ──
        batch_plot = script_dir / 'batch_plot_virus_depth.py'
        if batch_plot.is_file() and best_summary.exists():
            log.info("  [1/10 viz] Generating Stage 1 visualization plots...")
            run(f"python {batch_plot} --mode sample "
                f"-m {best_summary} -o {detect_dir / 'sample_distribution'}",
                log, "plot_sample_dist")
            run(f"python {batch_plot} --mode coabundance "
                f"-m {best_summary} -o {detect_dir / 'coabundance'}",
                log, "plot_coabundance")
            depth_stat_dir = detect_dir / "stat"
            if depth_stat_dir.exists():
                run(f"python {batch_plot} --mode depth "
                    f"-d {depth_stat_dir} -m {best_summary} "
                    f"-o {detect_dir / 'depth_plots'} -t {args.threads} -g",
                    log, "plot_depth_all")
            log.info("  [1/10 viz] Stage 1 plots complete -> %s", detect_dir)
        log.info("  Detection complete -> %s", detect_dir)

    # ═══════════════════════════════════════════
    # Stage 2: Filter (optional auto-filter)
    # ═══════════════════════════════════════════
    if args.filter or args.stage == "filter":
        _sh2 = add_stage_log(filter_dir, "2_filter")
        if not args.force and high_conf.exists():
            log.info("[2/10] Filter: checkpoint OK, skip")
        elif best_summary.exists():
            log.info("-" * 40)
            log.info("[2/10] High-Confidence Filtering")
            filter_parts = [
                f"python {script_dir / 'utils/filter_summary.py'}",
                f"-i {best_summary}",
                f"-o {high_conf}",
                f"-c {args.filter_cov}",
                f"-d {args.filter_depth}",
                f"-r {args.filter_reads}",
                f"--summary {filter_dir / 'filter_stats'}",
                "--plot",
            ]
            if args.filter_keyword:
                filter_parts.append(f"-k {args.filter_keyword}")
            if args.filter_tpm > 0:
                filter_parts.append(f"--min_tpm {args.filter_tpm}")
            if args.filter_poisson > 0:
                filter_parts.append(f"--min_poisson {args.filter_poisson}")
            if run(" ".join(filter_parts), log, "filter_summary"):
                log.info("  Filter complete -> %s", high_conf)
        else:
            log.warning("  best.summary not found, skipping filter")

    # ═══════════════════════════════════════════
    # Stage 3: Variant Analysis
    # ═══════════════════════════════════════════
    if args.stage in ("all", "variants"):
        _sh3 = add_stage_log(variants_dir, "3_variants")
        if not args.force and (variants_dir / "summary" / "all_summary.tsv").exists():
            log.info("[3/10] Variants: checkpoint OK, skip")
        else:
            log.info("-" * 40)
            log.info("[3/10] Variant Analysis")
            summary_in = get_summary()
            if not summary_in.exists():
                log.error("Summary not found: %s (run Stage 1 first)", summary_in)
                sys.exit(1)

            parts = [
                f"python {script_dir / 'batch_virus_variants.py'}",
                f"--summary {get_summary()}",
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
                parts += [
                    "--snpeff",
                    f"--snpeff_jar {args.snpeff_jar}",
                    f"--snpeff_config {args.snpeff_config}",
                    f"--snpeff_mem {args.snpeff_mem}",
                ]
            if args.snpgenie:
                parts.append("--snpgenie")
            if args.bam:
                parts += ["--bam", args.bam]
            else:
                parts += ["--fastq", str(reads)]
            if args.disable_dynamic_vcf:
                parts.append("--disable_dynamic_vcf")
            parts += [
                f"-q {args.vc_qual}",
                f"-d {args.vc_depth}",
                f"-f {args.vc_freq}",
                f"-a {args.vc_ambig}",
            ]
            if not args.no_resume and not args.force:
                parts.append("--resume")
            if not run(" ".join(parts), log, "batch_virus_variants"):
                log.warning("  Variant analysis partially failed, check logs")
            else:
                log.info("  Variants complete -> %s", variants_dir)

    # ═══════════════════════════════════════════
    # Stage 4: Full-length Assembly
    # ═══════════════════════════════════════════
    if args.stage in ("all", "full"):
        _sh4 = add_stage_log(full_dir, "4_full")
        if not args.force and full_dir.exists() and any(p.is_dir() for p in full_dir.iterdir()):
            log.info("[4/10] Assembly: checkpoint OK, skip")
        else:
            log.info("-" * 40)
            log.info("[4/10] Full-length Assembly")
            var_summary = variants_dir / "summary" / "all_summary.tsv"
            if not var_summary.exists():
                log.error("Variant summary not found: %s (run Stage 3 first)", var_summary)
                sys.exit(1)

            vsi = args.virus_full_script or str(script_dir / "virus-full.py")
            parts = [
                f"python {script_dir / 'batch_virus_full.py'}",
                f"--downstream_dir {variants_dir}",
                f"--summary {get_summary()}",
                f"--clean_data {reads}",
                f"--virus_full_script {vsi}",
                f"--outdir {full_dir}",
                f"--assembly_tools {args.assembly_tools}",
                f"--jobs {args.jobs}",
                f"--threads {args.threads}",
                f'--extra_args "{args.extra_args}"',
                f"--min_covered {args.min_covered}",
            ]
            if args.gb:
                parts.append(f"--gb {args.gb}")
            if not run(" ".join(parts), log, "batch_virus_full"):
                log.warning("  Assembly partially failed, check logs")
            else:
                log.info("  Assemblies complete -> %s", full_dir)

    # ═══════════════════════════════════════════
    # Stage 5: Extract Clean Assemblies
    # ═══════════════════════════════════════════
    if args.stage in ("all", "extract"):
        _sh5 = add_stage_log(extract_dir, "5_extract")
        if not args.force and extract_dir.exists() and any(p.is_dir() for p in extract_dir.iterdir()):
            log.info("[5/10] Extract: checkpoint OK, skip")
        else:
            log.info("-" * 40)
            log.info("[5/10] Extract Longest Contigs")
            if full_dir.exists():
                parts = [
                    f"python {script_dir / 'utils/extract_full_fasta.py'}",
                    f"--dir {full_dir}",
                    f"--outdir {extract_dir}",
                    f"--target_file {args.extract_target}",
                    f"--fill",
                    f"--ref_info {args.ref_info}",
                    f"--ref_dir {variants_dir}",
                    f"--max_n_genome {args.max_n_genome}",
                    f"--min_len {args.min_length}",
                    "--plot",
                ]
                if run(" ".join(parts), log, "extract_full_fasta"):
                    log.info("  Extraction complete -> %s", extract_dir)
            else:
                log.warning("  Assembly dir not found, skipping extract")

    # ═══════════════════════════════════════════
    # Stage 6: Post-hoc Visualization
    # ═══════════════════════════════════════════
    if args.stage in ("all", "post"):
        _sh6 = add_stage_log(post_dir, "6_post")
        if not args.force and post_dir.exists() and any(p.is_dir() for p in post_dir.iterdir()):
            log.info("[6/10] Post-hoc: checkpoint OK, skip")
        else:
            log.info("-" * 40)
            log.info("[6/10] Post-hoc Visualization")

            summary_for_post = get_summary()
            if not summary_for_post.exists():
                log.warning("  No summary found, skipping post-hoc")
            else:
                import pandas as pd

                df = pd.read_csv(summary_for_post, sep="\t")
                acc_col = next(
                    (
                        c
                        for c in ["Rep_Accession", "Accession", "Virus"]
                        if c in df.columns
                    ),
                    df.columns[0],
                )
                sp_col = next((c for c in ["Adjusted_Species", "Species_NCBI", "Species_ICTV"] if c in df.columns), None)
                virus_map = {}
                for _, row in df.iterrows():
                    acc = str(row.get(acc_col, ""))
                    if not acc: continue
                    sp = str(row.get(sp_col, acc)) if sp_col else acc
                    safe_name = sp.replace(" ", "_").replace("/", "_").replace("'", "")
                    virus_map[acc] = f"{safe_name}_{acc}"

                if len(virus_map) == 0:
                    log.warning("  No viruses to analyze")
                else:
                    log.info("  Viruses to process: %d", len(virus_map))
                    post_dir.mkdir(parents=True, exist_ok=True)

                    def process_one_virus(vname):
                        """Worker for parallel post-hoc analysis of a single virus."""
                        vname = str(vname)
                        vout = post_dir / virus_map.get(vname, vname)
                        vout.mkdir(parents=True, exist_ok=True)

                        vcf_in = find_virus_dir(variants_dir, "virus-variants", vname)
                        snpeff_in = find_virus_dir(variants_dir, "virus-SnpEff", vname)
                        sg_in = find_virus_dir(variants_dir, "virus-SNPGenie", vname)
                        acc = vname.split("_")[-1] if "_" in vname else vname
                        # Fix NC_ prefix loss: vname like "...tuber_viroid_NC_002030.1"
                        # split("_")[-1] gives "002030.1" (missing "NC_")
                        _parts = vname.split("_")
                        if len(_parts) >= 2 and _parts[-2].isalpha() and _parts[-2].isupper() and len(_parts[-2]) == 2:
                            acc = f"{_parts[-2]}_{_parts[-1]}"

                        tasks_done, tasks_total = 0, 0
                        if not args.skip_vcf_viz and vcf_in:
                            tasks_total += 1
                            if run(f"python {script_dir / 'virus_variants_analyzer.py'} "
                                   f"-i {vcf_in} -o {vout / 'vcf_viz'} -d {args.post_min_dp} "
                                   f"-f {args.post_min_af} -a {acc} -v {vname}", log, f"vcf_{vname}"):
                                tasks_done += 1

                        if not args.skip_vcf_merge and vcf_in:
                            tasks_total += 1
                            merge_flags = f"-d {vcf_in} -o {vout / 'vcf_merge'} --prefix {vname} --visualize"
                            if args.variant_caller == "ivar":
                                merge_flags += " --ivar"
                            if run(f"python {script_dir / 'virus_vcf_pipeline.py'} {merge_flags}", log, f"merge_{vname}"):
                                tasks_done += 1

                        if not args.skip_snpeff_macro and snpeff_in:
                            tasks_total += 1
                            if run(f"python {script_dir / 'snpeff_analysis.py'} "
                                   f"--miner {snpeff_in} --outdir {vout / 'snpeff_macro'}", log, f"eff_{vname}"):
                                tasks_done += 1

                        if not args.skip_maftools and snpeff_in:
                            tasks_total += 2
                            if run(f"python {script_dir / 'snpeff2maf.py'} "
                                   f"-i {snpeff_in} -minDP {args.post_min_dp} "
                                   f"-minAF {args.post_min_af} --filter-pass", log, f"maf_{vname}"):
                                tasks_done += 1
                            if run(f"Rscript {script_dir / 'viral_maftools.R'} "
                                   f"-i {snpeff_in} -o {vout / 'maftools'}", log, f"maftools_{vname}"):
                                tasks_done += 1

                        if not args.skip_snpgenie and sg_in:
                            tasks_total += 1
                            if run(f"python {script_dir / 'snpgenie_master.py'} "
                                   f"-i {sg_in} -o {vout / 'snpgenie'} -r {acc}", log, f"sg_{vname}"):
                                tasks_done += 1

                        return vname, tasks_done, tasks_total

                    with ThreadPoolExecutor(max_workers=min(len(virus_map), args.jobs)) as ex:
                        futures = {ex.submit(process_one_virus, v): v for v in virus_map}
                        for f in as_completed(futures):
                            name, done, total = f.result()
                            log.info("  %s: %d/%d analyses OK", name, done, total)

                    # Virus vs metadata association
                    meta_script = script_dir / 'utils' / 'virus_metadata_plot.py'
                    if meta_script.is_file():
                        run(f"python {meta_script} "
                            f"-m {summary_for_post} "
                            f"-o {post_dir / 'metadata_association'}",
                            log, "meta_association")

                log.info("  Post-hoc complete -> %s", post_dir)

    # ═══════════════════════════════════════════
    # Stage 7: Capheine Positive Selection
    # ═══════════════════════════════════════════
    if args.stage in ("all", "capheine"):
        _sh7 = add_stage_log(capheine_dir, "7_capheine")
        if not args.force and capheine_dir.exists() and any(p.is_dir() for p in capheine_dir.iterdir()):
            log.info("[7/10] Capheine: checkpoint OK, skip")
        else:
            log.info("-" * 40)
            log.info("[7/10] Positive Selection Analysis (Capheine)")

            cap_ref = args.capheine_ref
            cap_unaligned = args.capheine_unaligned

            # Auto-extract CDS from virus-annotations if not provided
            gb_dir = variants_dir / "virus-annotations"
            cap_input_dir = capheine_dir.parent / ".capheine_input"

            if (not cap_ref or not cap_unaligned) and gb_dir.exists():
                log.info("  Auto-extracting CDS from virus-annotations/...")
                cap_input_dir.mkdir(parents=True, exist_ok=True)

                # Find the virus to analyze from summary
                import pandas as _pd
                _df = _pd.read_csv(get_summary(), sep="\t")
                _acc_col = next((c for c in ["Rep_Accession", "Accession"] if c in _df.columns), _df.columns[0])
                _sp_col = next((c for c in ["Adjusted_Species", "Species", "Species_NCBI"] if c in _df.columns), None)
                _targets = _df[_acc_col].dropna().unique()

                # Build virus name mapping {acc: Species_Acc}, same as post-hoc
                _virus_map = {}
                for _, row in _df.drop_duplicates(subset=[_acc_col]).iterrows():
                    _a = str(row[_acc_col])
                    if _sp_col:
                        _sp = str(row[_sp_col]).replace(" ", "_").replace("/", "_").replace("'", "")
                        _virus_map[_a] = f"{_sp}_{_a}"
                    else:
                        _virus_map[_a] = _a

                for _acc in _targets:
                    _gb_file = gb_dir / f"{_acc}.gb"
                    if not _gb_file.exists():
                        _gb_file = gb_dir / f"{_acc.split('.')[0]}.gb"
                    if not _gb_file.exists():
                        continue

                    _vname = _virus_map.get(str(_acc), str(_acc))
                    _vout = cap_input_dir / _vname
                    _vout.mkdir(parents=True, exist_ok=True)
                    # Extract CDS from GB
                    run(f"python {script_dir / 'utils/gbk_extractor.py'} "
                        f"-i {_gb_file} -n {_vout / 'ref_cds.fasta'}", log, f"cds_{_acc}")

                    _ref_cds = _vout / "ref_cds.fasta"
                    _cap_ref = str(_ref_cds) if _ref_cds.exists() else None

                    # Collect assemblies as unaligned
                    import glob as _glob
                    _asm_matches = _glob.glob(str(extract_dir / f"*{_acc.split('.')[0]}*"))
                    _cap_unaligned = None
                    if _asm_matches:
                        _unaligned = _vout / "unaligned.fasta"
                        _all_seqs = []
                        for _fa in Path(_asm_matches[0]).rglob("*.full.fasta"):
                            _all_seqs.append(_fa.read_text())
                        if _all_seqs:
                            with open(_unaligned, "w") as _uf:
                                _uf.write("".join(_all_seqs))
                        if _unaligned.exists():
                            _cap_unaligned = str(_unaligned)

                    # Skip non-coding viruses (e.g. viroids)
                    if _ref_cds.exists() and _ref_cds.stat().st_size < 100:
                        log.info("  %s: no CDS (likely non-coding virus), skipped", _vname)
                        continue

                    if not _cap_ref or not _cap_unaligned:
                        continue

                    # Run capheine
                    _cap_out = capheine_dir / _vname
                    _cap_out.mkdir(parents=True, exist_ok=True)
                    _parts = [
                        f"python {script_dir / 'capheine_pipeline.py'}",
                        f"-r {_cap_ref}", f"-u {_cap_unaligned}",
                        f"-o {_cap_out}", f"--code {args.capheine_code}",
                        f"--workers {args.jobs}",
                        f"--cpus_iqtree {min(args.threads, 16)}",
                        f"--cpus_hyphy {min(args.threads, 32)}",
                    ]
                    if args.capheine_fg:
                        _parts.append(f"--foreground_list {args.capheine_fg}")
                    if run(" ".join(_parts), log, f"capheine_{_acc}"):
                        log.info("  %s: capheine OK", _acc)
                        # Visualize positive selection sites
                        drhip_csv = _cap_out / "drhip" / "combined_sites.csv"
                        cln_dir = _cap_out / "hyphy" / "CLN"
                        if drhip_csv.exists() and cln_dir.exists():
                            run(f"python {script_dir / 'utils/visual_codon_miner.py'} "
                                f"--drhip {drhip_csv} --clndir {cln_dir} "
                                f"-o {_cap_out / 'codon_plots'}", log, f"codon_{_acc}")

                        # Per-gene selection summary bar chart
                        if drhip_csv.exists():
                            _gen_bar = _cap_out / "selection_per_gene.pdf"
                            try:
                                import matplotlib
                                matplotlib.use('Agg')
                                import matplotlib.pyplot as __plt
                                _dr = __pd.read_csv(drhip_csv)
                                _gene_col = next((c for c in ['gene', 'Gene', 'gene_name'] if c in _dr.columns), None)
                                if _gene_col and len(_dr) > 0:
                                    _dr['gene_short'] = _dr[_gene_col].str.replace(
                                        r'.*\.part_', '', regex=True)
                                    _cnts = _dr['gene_short'].value_counts()
                                    _fig, _ax = __plt.subplots(figsize=(max(6, len(_cnts)*0.4), 5))
                                    _colors = __plt.cm.Set2(__plt.Normalize(0, max(len(_cnts)-1, 1))(range(len(_cnts))))
                                    _ax.barh(range(len(_cnts)), _cnts.values, color=_colors, edgecolor='#333')
                                    _ax.set_yticks(range(len(_cnts)))
                                    _ax.set_yticklabels(_cnts.index, fontsize=10, fontweight='bold')
                                    _ax.set_xlabel('Positive Selection Sites', fontweight='bold')
                                    _ax.set_title(f'Positive Selection Sites per Gene ({_acc})',
                                                 fontweight='bold', fontsize=13)
                                    for _j, _v in enumerate(_cnts.values):
                                        _ax.text(_v + max(_cnts.values)*0.02, _j, str(_v),
                                                va='center', fontweight='bold')
                                    __plt.tight_layout()
                                    _fig.savefig(_gen_bar, dpi=300, bbox_inches='tight')
                                    _gen_bar_png = str(_gen_bar).replace('.pdf', '.png')
                                    _fig.savefig(_gen_bar_png, dpi=300, bbox_inches='tight')
                                    __plt.close()
                                    log.info("  %s: gene selection plot -> %s", _acc, _gen_bar)
                            except Exception as _ex:
                                log.warning("  %s: gene selection plot failed: %s", _acc, _ex)
                    else:
                        log.warning("  %s: capheine failed", _acc)

                log.info("  Capheine complete -> %s", capheine_dir)

            elif cap_ref and cap_unaligned:
                capheine_dir.mkdir(parents=True, exist_ok=True)
                parts = [
                    f"python {script_dir / 'capheine_pipeline.py'}",
                    f"-r {cap_ref}", f"-u {cap_unaligned}",
                    f"-o {capheine_dir}", f"--code {args.capheine_code}",
                    f"--workers {args.jobs}",
                    f"--cpus_iqtree {min(args.threads, 16)}",
                    f"--cpus_hyphy {min(args.threads, 32)}",
                ]
                if args.capheine_fg:
                    parts.append(f"--foreground_list {args.capheine_fg}")
                if not run(" ".join(parts), log, "capheine_pipeline"):
                    log.warning("  Capheine analysis failed, check logs")
                else:
                    log.info("  Capheine complete -> %s", capheine_dir)
            else:
                log.warning("  capheine requires CDS input. Provide --capheine_ref/--capheine_unaligned")
                log.info("  or ensure virus-annotations/ and assemblies exist from Stage 3+5.")

    # ═══════════════════════════════════════════
    # Stage 8: Full-length Similarity Panorama
    # ═══════════════════════════════════════════
    if args.stage in ("all", "similarity"):
        _sh8 = add_stage_log(similarity_dir, "8_similarity")
        if not args.force and similarity_dir.exists() and any(p.is_dir() for p in similarity_dir.iterdir()):
            log.info("[8/10] Similarity: checkpoint OK, skip")
        else:
            log.info("-" * 40)
            log.info("[8/10] Full-length Similarity Panorama")

            consensus_base = variants_dir / "virus-consensus"
            if not consensus_base.exists():
                log.warning("  Consensus dir not found: %s", consensus_base)
            else:
                similarity_dir.mkdir(parents=True, exist_ok=True)
                sim_ref = args.sim_ref

                for vdir in consensus_base.iterdir():
                    if not vdir.is_dir(): continue
                    vname = vdir.name
                    vout = similarity_dir / vname
                    vout.mkdir(parents=True, exist_ok=True)

                    # Flatten nested consensus FASTA files into temp dir
                    import tempfile, shutil as _shutil
                    flat_dir = Path(tempfile.mkdtemp(prefix=f"sim_{vname}_", dir=str(_pipeline_tmp)))
                    fasta_files = list(vdir.rglob("*.consensus.fasta")) or list(vdir.rglob("*.fasta"))
                    for ff in fasta_files:
                        _shutil.copy(ff, flat_dir / f"{ff.parent.name}_{ff.name}")
                    if not list(flat_dir.glob("*")):
                        log.warning("  %s: no consensus FASTA found", vname)
                        _shutil.rmtree(flat_dir, ignore_errors=True)
                        continue

                    parts = [
                        f"python {script_dir / 'virus_auto_pipeline.py'}",
                        f"-i {flat_dir}",
                        f"-o {vout}",
                        f"--mode {args.sim_mode}",
                        f"--threads {args.threads}",
                    ]
                    if sim_ref:
                        parts.append(f"-g {sim_ref}")
                    if args.sim_cdhit:
                        parts.append("--cdhit")
                    if not args.no_resume and not args.force:
                        parts.append("--resume")

                    if run(" ".join(parts), log, f"similarity_{vname}"):
                        log.info("  %s: similarity OK", vname)
                    else:
                        log.warning("  %s: similarity failed", vname)
                    _shutil.rmtree(flat_dir, ignore_errors=True)

                log.info("  Similarity complete -> %s", similarity_dir)

    # ═══════════════════════════════════════════
    # Stage 9: DVG & Recombination Analysis
    # ═══════════════════════════════════════════
    if args.stage in ("all", "dvg"):
        _sh9 = add_stage_log(dvg_dir, "9_dvg")
        if not args.force and dvg_dir.exists() and any(p.is_dir() for p in dvg_dir.iterdir()):
            log.info("[9/10] DVG: checkpoint OK, skip")
        else:
            log.info("-" * 40)
            log.info("[9/10] DVG & Recombination Analysis")
            summary_in = get_summary()
            if not summary_in.exists():
                log.warning("  Summary not found, skipping DVG analysis")
            else:
                dvg_dir.mkdir(parents=True, exist_ok=True)
                dvg_reads = Path(args.dvg_reads).resolve() if args.dvg_reads else reads
                parts = [
                    f"python {script_dir / 'batch_virema_dvg.py'}",
                    f"-s {summary_in}",
                    f"-r {args.reference}",
                    f"-d {dvg_reads}",
                    f"-v {args.virema_script}",
                    f"--ref_info {args.ref_info}",
                    f"-o {dvg_dir}",
                    f"--seed {args.dvg_seed}",
                    f"--mindel {args.dvg_mindel}",
                    f"--min_cov {args.dvg_min_cov}",
                    f"-j {args.jobs}",
                    f"-t {args.threads}",
                ]
                if not args.no_resume and not args.force:
                    parts.append("--resume")
                if not run(" ".join(parts), log, "batch_virema_dvg"):
                    log.warning("  DVG analysis failed, check logs")
                else:
                    log.info("  DVG complete -> %s", dvg_dir)

    # ═══════════════════════════════════════════
    # Stage 10: Generate Summary Report
    # ═══════════════════════════════════════════
    if args.stage in ("all", "report"):
        _sh10 = add_stage_log(report_dir, "10_report")
        log.info("-" * 40)
        log.info("[10/10] Generate Pipeline Summary Report")
        parts = [
            f"python {script_dir / 'generate_pipeline_report.py'}",
            f"-d {out}",
            f"-o {report_dir / 'Pipeline_Summary_Report.html'}",
        ]
        if args.ai_api_key:
            parts.append(f"--ai-api-key {args.ai_api_key}")
            parts.append(f"--ai-model {args.ai_model}")
        if run(" ".join(parts), log, "generate_report"):
            log.info("  Report generated -> %s", report_dir / "Pipeline_Summary_Report.html")
        else:
            log.warning("  Report generation failed")

    # ---- Done ----
    log.info("=" * 55)
    log.info("Pipeline complete! | %s", datetime.now().strftime("%H:%M:%S"))
    log.info("=" * 55)


if __name__ == "__main__":
    main()
