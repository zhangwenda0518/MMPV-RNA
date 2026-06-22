#!/usr/bin/env python3
"""
virome_analysis.py — 病毒基因组下游分析 + GenBank 提交准备 v2.1
============================================================
输入: all_plant_viruses.fasta (来自 08_Rescue)
输出: 09_Virome_Analysis/

阶段:
  1. suvtk taxonomy     — ICTV 分类学注释
  2. suvtk features     — CDS/tRNA/结构蛋白注释 (生成 .tbl)
  3. hypothetical       — DIAMOND/HMM 假想蛋白深度注释
  4. submission (Sequin) — .fsa + .tbl + .cmt → .sqn (tbl2asn)

用法:
  python virome_analysis.py -i out/08_Rescue/all_plant_viruses.fasta -o out/09_Virome_Analysis/ -t 40
"""

import argparse, os, sys, subprocess, logging, time, shutil
from pathlib import Path
from datetime import datetime
from shutil import which

SCRIPT_DIR = Path(__file__).resolve().parent


def setup_logger(out_dir):
    logger = logging.getLogger("VirusAnalysis")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('[%H:%M:%S] %(levelname)s %(message)s'))
    logger.addHandler(ch)
    os.makedirs(out_dir, exist_ok=True)
    fh = logging.FileHandler(os.path.join(out_dir, 'analysis.log'), encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s %(message)s'))
    logger.addHandler(fh)
    return logger


def run(cmd, log, step_name):
    log.info("[%s] %s", step_name, cmd[:150])
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        log.info("[%s] ✓", step_name)
        return True
    except subprocess.CalledProcessError as e:
        log.error("[%s] ✗ (exit=%d): %s", step_name, e.returncode, e.stderr[:200])
        return False


# ══════════════════════════════════════════════════════════════
# Stage 1: suvtk taxonomy
# ══════════════════════════════════════════════════════════════

def run_taxonomy(fasta, out_dir, db, threads, log):
    tax_out = out_dir / "suvtk_taxonomy"
    tax_tsv = tax_out / "taxonomy.tsv"
    if tax_tsv.is_file() and tax_tsv.stat().st_size > 100:
        log.info("[1/5] taxonomy — 跳过 (已存在)")
        return tax_out
    tax_out.mkdir(parents=True, exist_ok=True)
    ok = run(f"suvtk taxonomy -i {fasta} -o {tax_out} -d {db} -s 0.7 -t {threads}",
             log, "suvtk taxonomy")
    return tax_out if ok else None


# ══════════════════════════════════════════════════════════════
# Stage 2: suvtk features (CDS annotation)
# ══════════════════════════════════════════════════════════════

def run_features(fasta, out_dir, tax_dir, db, threads, log):
    feat_out = out_dir / "suvtk_features"
    tbl = feat_out / "featuretable.tbl"
    if tbl.is_file() and tbl.stat().st_size > 100:
        log.info("[2/5] features — 跳过 (已存在)")
        return feat_out
    feat_out.mkdir(parents=True, exist_ok=True)
    cmd = f"suvtk features -i {fasta} -o {feat_out} -d {db} --coding-complete -t {threads}"
    if tax_dir:
        tax_tsv = tax_dir / "taxonomy.tsv"
        if tax_tsv.is_file():
            cmd += f" --taxonomy {tax_tsv}"
    ok = run(cmd, log, "suvtk features")
    return feat_out if ok else None


# ══════════════════════════════════════════════════════════════
# Stage 3: Hypothetical protein analysis
# ══════════════════════════════════════════════════════════════

def run_hypothetical(out_dir, feat_dir, rvdb_dir, threads, log):
    hypo_out = out_dir / "hypothetical"
    updated_tbl = hypo_out / "featuretable_updated.tbl"
    if updated_tbl.is_file() and updated_tbl.stat().st_size > 100:
        log.info("[3/5] hypothetical — 跳过 (已存在)")
        return hypo_out

    tbl = feat_dir / "featuretable.tbl"
    faa = feat_dir / "proteins.faa"
    if not tbl.is_file() or not faa.is_file():
        log.warning("[3/5] hypothetical — 缺少 .tbl/.faa, 跳过")
        return None

    hypo_out.mkdir(parents=True, exist_ok=True)
    rvdb = Path(rvdb_dir)
    rvdb_fasta = rvdb / "U-RVDBv31.0-prot.EX.acc.fasta" if rvdb else None
    rvdb_hmm = rvdb / "U-RVDBv31.0-prot.hmm" if rvdb else None
    rvdb_annot = rvdb / "U-RVDBv31.0-prot.info.tab" if rvdb else None

    script = SCRIPT_DIR / "analyze_hypothetical.py"
    if not script.is_file():
        script = SCRIPT_DIR.parent / "suvtk_submission" / "analyze_hypothetical.py"

    cmd = f"python {script} -t {tbl} -f {faa} -o {hypo_out} --threads {threads}"
    if rvdb_fasta and rvdb_fasta.is_file():
        cmd += f" --diamond -d {rvdb_fasta}"
    if rvdb_hmm and rvdb_hmm.is_file():
        cmd += f" --hmmer --hmmer-db {rvdb_hmm}"
    if rvdb_annot and rvdb_annot.is_file():
        cmd += f" --annot {rvdb_annot}"

    ok = run(cmd, log, "analyze_hypothetical")
    return hypo_out if ok else None


# ══════════════════════════════════════════════════════════════
# Stage 4: Sequin submission files
# ══════════════════════════════════════════════════════════════

def run_submission(fasta, out_dir, feat_dir, tax_dir, hypo_dir, log):
    sub_out = out_dir / "submission"
    sqn = sub_out / "submission.sqn"
    if sqn.is_file() and sqn.stat().st_size > 100:
        log.info("[4/5] submission — 跳过 (已存在)")
        return sub_out
    sub_out.mkdir(parents=True, exist_ok=True)

    fsa = sub_out / "sequences.fsa"
    tbl = sub_out / "featuretable.tbl"
    cmt = sub_out / "comments.cmt"

    # Build .fsa
    if not fsa.is_file():
        import shutil
        shutil.copy(fasta, fsa)
    # Build .tbl
    if not tbl.is_file():
        src_tbl = feat_dir / "featuretable.tbl"
        if hypo_dir:
            updated = hypo_dir / "featuretable_updated.tbl"
            if updated.is_file():
                src_tbl = updated
        if src_tbl.is_file():
            import shutil
            shutil.copy(src_tbl, tbl)
    # Build .cmt — GenBank submission comments
    if not cmt.is_file():
        with open(cmt, 'w') as f:
            f.write("Submitter Comments\n")
            f.write("Assembly Method: MMPV-RNA v2.3 pipeline\n")
            f.write(f"Sequencing Technology: Illumina\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d')}\n")

    # tbl2asn: generate .sqn
    cmd = f"tbl2asn -p {sub_out} -t template.sbt -i {fsa} -o {sub_out}"
    log.info("[4/5] tbl2asn: %s", cmd)

    return sub_out


# ══════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════

def main():
    p = argparse.ArgumentParser(description="病毒基因组下游分析 + NCBI 提交")
    p.add_argument("-i", "--input", required=True, help="输入: all_plant_viruses.fasta")
    p.add_argument("-o", "--output", default="./09_Virome_Analysis", help="输出目录")
    p.add_argument("-t", "--threads", type=int, default=40)
    p.add_argument("--suvtk-db", default=os.path.expanduser("~/database/virus-db/suvtk_db/"))
    p.add_argument("--rvdb-dir", default=os.path.expanduser("~/database/virus-db/RVDB-v31/"))
    p.add_argument("--skip-suvtk", action="store_true", help="跳过 suvtk 步骤")
    args = p.parse_args()

    inp = Path(args.input)
    if not inp.is_file():
        sys.exit(f"输入文件不存在: {inp}")
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    log = setup_logger(out)
    log.info("=" * 60)
    log.info("Virome Analysis Pipeline v1.0")
    log.info("  Input:  %s (%d 条序列)", inp, sum(1 for _ in open(inp) if _.startswith('>')))
    log.info("  Output: %s", out)
    log.info("  Threads: %d", args.threads)
    log.info("=" * 60)

    t0 = time.time()

    # 1-4. suvtk pipeline
    if not args.skip_suvtk:
        # 1. Taxonomy
        tax_dir = run_taxonomy(inp, out, args.suvtk_db, args.threads, log)
        # 2. Features
        feat_dir = run_features(inp, out, tax_dir, args.suvtk_db, args.threads, log)
        # 3. Hypothetical
        hypo_dir = None
        if feat_dir:
            hypo_dir = run_hypothetical(out, feat_dir, args.rvdb_dir, args.threads, log)
        # 4. Submission
        if feat_dir:
            sub_dir = run_submission(inp, out, feat_dir, tax_dir, hypo_dir, log)
    else:
        log.info("[1-4/5] suvtk pipeline — 跳过")

    # 5. 生成 ref_info.tsv (供 virome_analysis_pipeline/auto_known_virus.py 使用)
    ref_info = out / "ref_info.tsv"
    if not ref_info.is_file():
        ref_builder = SCRIPT_DIR / "build_ref_info.py"
        if ref_builder.is_file():
            pipeline_root = Path(args.output).parent if Path(args.output).name == "09_Virome_Analysis" else Path(args.output)
            run(f"python {ref_builder} -o {pipeline_root}", log, "build ref_info.tsv")

    elapsed = time.time() - t0
    log.info("=" * 60)
    log.info("完成! 耗时 %.0fs → %s", elapsed, out)
    log.info("  suvtk_taxonomy/       ICTV 分类注释")
    log.info("  suvtk_features/       基因注释 (CDS/tRNA)")
    log.info("  hypothetical/          假想蛋白深度注释")
    log.info("  submission/            GenBank 提交文件 (.sqn)")
    log.info("  ref_info.tsv           已知病毒分析 info 文件")
    log.info("=" * 60)


if __name__ == "__main__":
    main()
