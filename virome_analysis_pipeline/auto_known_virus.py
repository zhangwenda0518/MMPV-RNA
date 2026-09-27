#!/usr/bin/env python3
"""
auto_known_virus.py - Known Virus Analysis Pipeline
====================================================
9-stage automated pipeline for known virus detection, variant analysis,
full-length assembly, and post-hoc characterization.

Stages:
  1. detect    - Rapid quantification (batch_virus_depth.py)
  2. filter    - High-confidence filtering (filter_summary.py) [auto with --filter]
  3. variants  - Variant calling + SnpEff + SnpGenie (batch_virus_variants.py)
  4. full      - De novo full-length assembly (virus-full.py)
  5. extract   - Extract longest contigs (extract_full_fasta.py)
  6. post      - VCF visualization + SnpEff macro + MAF + SnpGenie
  7. similarity - Full-length similarity panorama (virus_auto_pipeline.py)
  8. dvg        - DVG & recombination detection (batch_virema_dvg.py)
  9. report     - Generate summary report + AI interpretation prompts

Note: positive selection (capheine) analysis migrated to virome_phylo_pipeline/.

Output structure (legacy 布局; standard 布局见 doc/IO_LAYOUT_DESIGN.md):
  output_dir/
    01_detection/       Stage 1: detection results
    02_filtering/       Stage 2: high-confidence filtering
    03_variants/        Stage 3: variant analysis
    04_post_analysis/   Stage 4: post-hoc viz (VCF merge, PCA, heatmap)
    05_assembly/        Stage 5: full assemblies
    06_extraction/      Stage 6: extracted contigs
    07_similarity/      Stage 7: similarity panorama
    08_dvg/             Stage 8: DVG & recombination
    09_report/          Stage 9: HTML report
    logs/               Pipeline logs
"""

import argparse
import subprocess
import sys
import os
import logging
import shutil
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

SCRIPT_DIR = Path(__file__).resolve().parent

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根)
_REPO_ROOT = SCRIPT_DIR.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
from mmpv_common.io_layout import build_analysis_dirs, normalize_layout_env

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


# ── Schema constants: required columns for inter-stage TSV files ──
STAGE1_SCHEMA = ["Sample", "Rep_Accession", "Rep_Coverage(%)", "Rep_MeanDepth",
                  "Asm_TPM", "Poisson_Ratio", "taxid", "Adjusted_Species"]
STAGE3_SCHEMA = ["Sample", "Accession", "Covered%", "MeanDepth", "Consensus"]


def check_tsv_schema(tsv_path, required_cols, stage_label, logger):
    """Validate that a TSV file has all required columns.
    Raises SystemExit with a clear message if columns are missing.
    """
    if not tsv_path.exists():
        logger.error("[%s] Schema check FAILED: file not found: %s", stage_label, tsv_path)
        logger.error("  → Has the upstream stage completed successfully?")
        sys.exit(1)
    try:
        with open(tsv_path, encoding="utf-8") as f:
            actual = set(f.readline().strip().split("\t"))
    except Exception as e:
        logger.error("[%s] Schema check FAILED: cannot read %s: %s", stage_label, tsv_path, e)
        sys.exit(1)
    missing = [c for c in required_cols if c not in actual]
    if missing:
        logger.error("[%s] Schema check FAILED: missing columns in %s", stage_label, tsv_path.name)
        logger.error("  Required: %s", ", ".join(required_cols))
        logger.error("  Missing:  %s", ", ".join(missing))
        logger.error("  Actual:   %s", ", ".join(sorted(actual)))
        logger.error("  → Column names may have changed between pipeline versions. Check upstream script output.")
        sys.exit(1)
    logger.info("[%s] Schema check PASSED: %d columns (%d required present)",
                stage_label, len(actual), len(required_cols))


def main():
    parser = argparse.ArgumentParser(
        description="Known Virus Analysis Pipeline (10-stage)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ---- I/O ----
    g = parser.add_argument_group("Input/Output")
    g.add_argument("--reads_dir", default=None, help="Clean reads directory")
    g.add_argument("--output_dir", "-o", default=None, help="Output root directory")
    g.add_argument("--io-layout", choices=["legacy", "standard"], default=None,
                   help="I/O directory layout: legacy=v3.0 names (default); "
                        "standard=unified numbering (doc/IO_LAYOUT_DESIGN.md); "
                        "or set env MMPV_IO_LAYOUT")
    g.add_argument("--ref_info", default=None, help="Reference info TSV")
    g.add_argument("--reference", default=None, help="Reference genome FASTA")

    # ---- Stage control ----
    g = parser.add_argument_group("Stage Control")
    g.add_argument(
        "--stage",
        default="all",
        choices=["all", "detect", "filter", "variants", "full", "extract", "post", "similarity", "dvg", "report"],
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
    g.add_argument("--filter", action=argparse.BooleanOptionalAction, default=True, help="Enable high-confidence filtering (default: on, use --no-filter to disable)")
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
    g.add_argument("--meta", type=str, default=None,
                   help="SRA metadata TSV with Lat/Lon cols for Mantel test (passed to virus_vcf_pipeline)")
    g.add_argument("--metadata_tsv", type=str, default=None,
                   help="队列元数据表 (Global_Unified_Metadata_Core14.tsv). 提供后 post 阶段改用发表级 v2 病毒×元数据分析引擎")
    g.add_argument("--skip_metadata_assoc", action="store_true",
                   help="Skip virus × metadata association plots in post stage")
    g.add_argument("--skip_snpeff_macro", action="store_true")
    g.add_argument("--skip_maftools", action="store_true")
    g.add_argument("--skip_snpgenie", action="store_true")

    # ---- Stage 7: Similarity ----
    g = parser.add_argument_group("Stage 8: Similarity Panorama (virus_auto_pipeline)")
    g.add_argument("--sim_ref", help="GenBank accession or .gb file for similarity analysis")
    g.add_argument("--sim_mode", default="filter", choices=["strict", "filter", "fill", "all"])
    g.add_argument("--sim_cdhit", action="store_true", help="Enable CD-HIT dedup")

    # ---- Stage 8: DVG ----
    g = parser.add_argument_group("Stage 8: DVG & Recombination (batch_virema_dvg)")
    g.add_argument("--virema_script", default=str(SCRIPT_DIR / "../biosoft/virema/ViReMa.py"), help="Path to ViReMa.py")
    g.add_argument("--dvg_seed", type=int, default=25, help="ViReMa seed length (default: 25)")
    g.add_argument("--dvg_mindel", type=int, default=15, help="Microdeletion threshold (default: 15)")
    g.add_argument("--dvg_min_cov", type=float, default=80.0, help="Min coverage%% for DVG analysis (default: 80)")
    g.add_argument("--dvg_shm", action="store_true", help="Use /dev/shm RAM disk for ViReMa")
    g.add_argument("--dvg_reads", default=None, help="FASTQ reads dir for DVG (default: same as --reads_dir)")

    # ---- Stage 9: Report ----
    g = parser.add_argument_group("Stage 9: Report Generation (generate_pipeline_report)")
    g.add_argument("--report_ai", action="store_true", help="Include AI interpretation prompts")
    g.add_argument("--ai_api_key", default=os.environ.get("DEEPSEEK_API_KEY") or os.environ.get("AI_API_KEY"),
                   help="DeepSeek/OpenAI API key for AI interpretation (default: $DEEPSEEK_API_KEY/$AI_API_KEY)")
    g.add_argument("--ai_model", default="deepseek-chat", help="LLM model")

    # ---- Profile support ----
    g = parser.add_argument_group("Profile")
    g.add_argument("--profile", default=None, help="Profile name within pipeline_config.yaml (e.g. 'analysis')")
    g.add_argument("--config", default=None, help="Path to pipeline_config.yaml (auto-detected if omitted)")

    # Pre-parse only --profile from raw argv
    # ── Config / Profile loading (YAML) ──
    config_path = None
    profile_name = None
    for i, a in enumerate(sys.argv[1:], 1):
        if a == '--config' and i < len(sys.argv):
            config_path = sys.argv[i + 1]
        elif a.startswith('--config='):
            config_path = a.split('=', 1)[1]
        elif a == '--profile' and i < len(sys.argv):
            profile_name = sys.argv[i + 1]
        elif a.startswith('--profile='):
            profile_name = a.split('=', 1)[1]

    if not config_path:
        candidate = SCRIPT_DIR.parent / "pipeline_config.yaml"
        if candidate.exists():
            config_path = str(candidate)

    if config_path and profile_name:
        try:
            import yaml
            with open(config_path, 'r') as f:
                full_config = yaml.safe_load(f)
            profiles = full_config.get("profiles", {})
            if profile_name not in profiles:
                print(f"[WARNING] profile '{profile_name}' not found in {config_path}. Available: {list(profiles.keys())}", file=sys.stderr)
            else:
                profile_data = profiles[profile_name]
                def _flatten(d, prefix=""):
                    items = {}
                    for k, v in d.items():
                        key = f"{prefix}{k}" if prefix else k
                        if isinstance(v, dict) and not k.startswith("_"):
                            items.update(_flatten(v, f"{key}_"))
                        else:
                            # Set both prefixed and bare key: argparse args don't have prefixes
                            items[key] = v
                            if prefix and k not in items:
                                items[k] = v
                    return items
                flat = {k: v for k, v in _flatten(profile_data).items()
                        if v is not None and not isinstance(v, (dict, list))}
                parser.set_defaults(**flat)
                print(f"[INFO] Loaded profile '{profile_name}' from {config_path} ({len(flat)} keys)", file=sys.stderr)
        except ImportError:
            print(f"[WARNING] pyyaml not installed, skipping config loading", file=sys.stderr)
        except Exception as e:
            print(f"[WARNING] Config loading failed: {e}", file=sys.stderr)

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

    # 布局解析 (CLI > MMPV_IO_LAYOUT > legacy) 并写回环境供子进程继承
    normalize_layout_env(getattr(args, "io_layout", None))

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

    # ── Process Guard: kill entire process tree on Ctrl+C ──
    try:
        from process_guard import install_global_guard
        _guard = install_global_guard(logger=log)
    except ImportError:
        _guard = None
        log.warning("[Guard] process_guard module not found — Ctrl+C may leave orphans")

    # ---- Dry-run ----
    if args.dry_run:
        s = args.stage
        log.info("")
        log.info("=== DRY-RUN ===")
        log.info("  Tool:    %s", args.tool)
        if s in ("all", "detect"):   log.info("  [1/9] Detect:   batch_virus_depth.py")
        if args.filter or s == "filter": log.info("  [2/9] Filter:   filter_summary.py")
        if s in ("all", "variants"):  log.info("  [3/9] Variants: batch_virus_variants.py (caller=%s snpeff=%s snpgenie=%s)", args.variant_caller, args.snpeff, args.snpgenie)
        if s in ("all", "post"):      log.info("  [4/9] Post-hoc: VCF viz + SnpEff + MAF + SnpGenie")
        if s in ("all", "full"):      log.info("  [5/9] Full:     batch_virus_full.py")
        if s in ("all", "extract"):   log.info("  [6/9] Extract:  extract_full_fasta.py")
        if s in ("all", "similarity"): log.info("  [7/9] Similarity: virus_auto_pipeline.py")
        if s in ("all", "dvg"):       log.info("  [8/9] DVG:    batch_virema_dvg.py")
        if s in ("all", "report"):    log.info("  [9/9] Report: generate_pipeline_report.py")
        log.info("=== DRY-RUN END ===")
        return

    try:
        _run_pipeline_stages(args, out, reads, log, script_dir, _pipeline_tmp)
    finally:
        if _guard:
            try:
                from process_guard import uninstall_global_guard
                uninstall_global_guard()
                log.info("[Guard] Process tree cleanup complete.")
            except Exception:
                pass


def _run_pipeline_stages(args, out, reads, log, script_dir, _pipeline_tmp):
    """Execute all pipeline stages. Extracted to enable try/finally guard cleanup."""
    # ---- Shared paths (目录名由 mmpv_common.io_layout 决定, legacy=v3.0 现行名) ----
    _dirs = build_analysis_dirs(out, getattr(args, 'io_layout', None))
    detect_dir = _dirs['detect']
    filter_dir = _dirs['filter']
    variants_dir = _dirs['variants']
    post_dir = _dirs['post']
    assembly_dir = _dirs['assembly']
    extract_dir = _dirs['extract']
    similarity_dir = _dirs['similarity']
    dvg_dir = _dirs['dvg']
    report_dir = _dirs['report']

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

    def _stage_ok(stage_dir):
        return (stage_dir / "run.ok").is_file()

    def _mark_stage_ok(stage_dir):
        stage_dir.mkdir(parents=True, exist_ok=True)
        (stage_dir / "run.ok").write_text("ok")

    # ═══════════════════════════════════════════
    # Stage 1: Detection
    # ═══════════════════════════════════════════
    if args.stage in ("all", "detect"):
        _sh = add_stage_log(detect_dir, "1_detect")
        parts = []
        if not args.force and _stage_ok(detect_dir):
            log.info("[1/9] Detection: checkpoint OK, skip")
        else:
            log.info("-" * 40)
            log.info("[1/9] Rapid Virus Detection")
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
        if parts:
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
        _mark_stage_ok(detect_dir)

    # ═══════════════════════════════════════════
    # Stage 2: Filter (optional auto-filter)
    # ═══════════════════════════════════════════
    if args.filter or args.stage == "filter":
        _sh2 = add_stage_log(filter_dir, "2_filter")
        if not args.force and _stage_ok(filter_dir):
            log.info("[2/9] Filter: checkpoint OK, skip")
        elif best_summary.exists():
            check_tsv_schema(best_summary, STAGE1_SCHEMA, "Stage2→Stage1", log)
            log.info("-" * 40)
            log.info("[2/9] High-Confidence Filtering")
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
                _mark_stage_ok(filter_dir)
        else:
            log.warning("  best.summary not found, skipping filter")

    # ═══════════════════════════════════════════
    # Stage 3: Variant Analysis
    # ═══════════════════════════════════════════
    if args.stage in ("all", "variants"):
        _sh3 = add_stage_log(variants_dir, "3_variants")
        if not args.force and _stage_ok(variants_dir):
            log.info("[3/9] Variants: checkpoint OK, skip")
        else:
            log.info("-" * 40)
            log.info("[3/9] Variant Analysis")
            summary_in = get_summary()
            check_tsv_schema(summary_in, STAGE1_SCHEMA, "Stage3→upstream", log)
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
                _mark_stage_ok(variants_dir)

    # ═══════════════════════════════════════════
    # Stage 6: Post-hoc Visualization
    # ═══════════════════════════════════════════
    if args.stage in ("all", "post"):
        _sh4 = add_stage_log(post_dir, "4_post")
        if not args.force and _stage_ok(post_dir):
            log.info("[4/9] Post-hoc: checkpoint OK, skip")
        else:
            log.info("-" * 40)
            log.info("[4/9] Post-hoc Visualization")

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
                            merge_flags = (f"-d {vcf_in} -o {vout / 'vcf_merge'} --prefix {vname} "
                                          f"--visualize --qc --snp-matrix --tree "
                                          f"--dist-metrics both --pca-method genotype --ld")
                            if args.variant_caller == "ivar":
                                merge_flags += " --ivar"
                            if getattr(args, 'meta', None):
                                merge_flags += f" --meta {args.meta}"
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
                    # - --metadata_tsv given -> publication-grade v2 engine (repo root)
                    # - otherwise            -> utils v1 (auto-fetches SRA metadata via --sra_list)
                    if not getattr(args, 'skip_metadata_assoc', False):
                        meta_v2 = script_dir / 'virus_metadata_plot.py'
                        meta_v1 = script_dir / 'utils' / 'virus_metadata_plot.py'
                        if getattr(args, 'metadata_tsv', None) and meta_v2.is_file():
                            run(f"python {meta_v2} "
                                f"-v {summary_for_post} "
                                f"-m {args.metadata_tsv} "
                                f"-o {post_dir / 'metadata_association'}",
                                log, "meta_association_v2")
                        elif meta_v1.is_file():
                            # sra.list 取自 summary 第一列 (Sample/Run 编号)
                            sra_col = next((c for c in ("Sample", "sample", "Run", "run") if c in df.columns), df.columns[0])
                            sra_list = post_dir / 'sra.list'
                            try:
                                df[sra_col].dropna().astype(str).str.strip() \
                                    .loc[lambda s: s != ""].drop_duplicates() \
                                    .to_csv(sra_list, index=False, header=False)
                            except Exception as e:
                                log.warning("  sra.list 生成失败: %s", e)
                            run(f"python {meta_v1} "
                                f"-v {summary_for_post} "
                                f"--sra_list {sra_list} "
                                f"-o {post_dir / 'metadata_association'}",
                                log, "meta_association")

                log.info("  Post-hoc complete -> %s", post_dir)
                _mark_stage_ok(post_dir)

    # ═══════════════════════════════════════════
    # Stage 4: Full-length Assembly
    # ═══════════════════════════════════════════
    if args.stage in ("all", "full"):
        _sh5 = add_stage_log(assembly_dir, "5_assembly")
        if not args.force and _stage_ok(assembly_dir):
            log.info("[5/9] Assembly: checkpoint OK, skip")
        else:
            log.info("-" * 40)
            log.info("[5/9] Full-length Assembly")
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
                f"--outdir {assembly_dir}",
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
                log.info("  Assemblies complete -> %s", assembly_dir)
                _mark_stage_ok(assembly_dir)

    # ═══════════════════════════════════════════
    # Stage 5: Extract Clean Assemblies
    # ═══════════════════════════════════════════
    if args.stage in ("all", "extract"):
        _sh6 = add_stage_log(extract_dir, "6_extract")
        if not args.force and _stage_ok(extract_dir):
            log.info("[6/9] Extract: checkpoint OK, skip")
        else:
            log.info("-" * 40)
            log.info("[6/9] Extract Longest Contigs")
            if assembly_dir.exists():
                parts = [
                    f"python {script_dir / 'utils/extract_full_fasta.py'}",
                    f"--dir {assembly_dir}",
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
                    _mark_stage_ok(extract_dir)
            else:
                log.warning("  Assembly dir not found, skipping extract")

    # Stage 7: Full-length Similarity Panorama
    # ═══════════════════════════════════════════
    if args.stage in ("all", "similarity"):
        _sh7 = add_stage_log(similarity_dir, "7_similarity")
        if not args.force and _stage_ok(similarity_dir):
            log.info("[7/9] Similarity: checkpoint OK, skip")
        else:
            log.info("-" * 40)
            log.info("[7/9] Full-length Similarity Panorama")

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
                _mark_stage_ok(similarity_dir)

    # ═══════════════════════════════════════════
    # Stage 8: DVG & Recombination Analysis
    # ═══════════════════════════════════════════
    if args.stage in ("all", "dvg"):
        _sh8 = add_stage_log(dvg_dir, "8_dvg")
        if not args.force and _stage_ok(dvg_dir):
            log.info("[8/9] DVG: checkpoint OK, skip")
        else:
            log.info("-" * 40)
            log.info("[8/9] DVG & Recombination Analysis")
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
                    _mark_stage_ok(dvg_dir)

    # ═══════════════════════════════════════════
    # Stage 9: Generate Summary Report
    # ═══════════════════════════════════════════
    if args.stage in ("all", "report"):
        _sh9 = add_stage_log(report_dir, "9_report")
        log.info("-" * 40)
        log.info("[9/9] Generate Pipeline Summary Report")
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


# ── Pipeline topology registry (documentation + programmatic introspection) ──
# Each entry: (stage_num, name, script, input_from, output_dir_key, checkpoint)
STAGES = [
    (1,  "detect",     "batch_virus_depth.py",     None,            "detect_dir",     "best.summary.tsv"),
    (2,  "filter",     "filter_summary.py",         "detect_dir",    "filter_dir",     "high_conf.summary.tsv"),
    (3,  "variants",   "batch_virus_variants.py",   "get_summary()", "variants_dir",   "all_summary.tsv"),
    (4,  "post",       "virus_vcf_pipeline.py",     "variants_dir",  "post_dir",       "figs/"),
    (5,  "assembly",   "batch_virus_full.py",       "variants_dir",  "assembly_dir",   "*/"),
    (6,  "extract",    "extract_full_fasta.py",     "assembly_dir",  "extract_dir",    "*.fasta"),
    (7,  "similarity", "virus_auto_pipeline.py",    "variants_dir",  "similarity_dir", "pipeline_results/"),
    (8,  "dvg",        "batch_virema_dvg.py",       "get_summary()", "dvg_dir",       "Summary_Analysis_Report/"),
    (9,  "report",     "generate_pipeline_report.py","*",             "report_dir",     "Pipeline_Summary_Report.html"),
]


if __name__ == "__main__":
    main()
