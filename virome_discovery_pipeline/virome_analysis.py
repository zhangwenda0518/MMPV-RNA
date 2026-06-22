#!/usr/bin/env python3
"""
virome_analysis.py — 病毒基因组下游分析 + GenBank 提交准备 v2.0
============================================================
输入: all_plant_viruses.fasta (来自 08_Rescue)
输出: 09_Virome_Analysis/

阶段:
  1. Cenote-Taker3      — 病毒 hallmark 基因检测 + 功能注释 + 分类
  2. suvtk taxonomy     — ICTV 分类学注释
  3. suvtk features     — CDS/tRNA/结构蛋白注释 (生成 .tbl)
  4. hypothetical       — DIAMOND/HMM 假想蛋白深度注释
  5. submission (Sequin) — .fsa + .tbl + .cmt → .sqn (tbl2asn)
  6. summary report     — 分析总结

用法:
  python virome_analysis.py -i out/08_Rescue/all_plant_viruses.fasta -o out/09_Virome_Analysis/ -t 40
  python virome_analysis.py -i out/08_Rescue/all_plant_viruses.fasta -o out/09_Virome_Analysis/ -t 40 --skip-cenote
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
# Stage 0: Cenote-Taker3 (病毒 hallmark 基因 + 功能注释 + 分类)
# ══════════════════════════════════════════════════════════════

def run_cenote(fasta, out_dir, threads, log):
    """Cenote-Taker3: virus-specific annotation pipeline"""
    cenote_bin = which("cenotetaker3") or which("cenote-taker3")
    if not cenote_bin:
        log.warning("[0/5] Cenote-Taker3 未安装, 跳过")
        return None

    # 自动检测数据库
    cenote_dbs = os.environ.get("CENOTE_DBS", "")
    if not cenote_dbs:
        for p in [os.path.expanduser("~/database/virus-db/ct3_DBs"),
                  os.path.expanduser("~/database/virus-db/ct3_DB")]:
            if os.path.isdir(p):
                cenote_dbs = p; break
    if not cenote_dbs:
        log.warning("[0/5] Cenote-Taker3 DB 未找到, 跳过 (设置 CENOTE_DBS 或安装到 ~/database/virus-db/ct3_DBs)")
        return None

    ct3_out = out_dir / "cenote_taker3"
    summary_file = ct3_out / "run_summary.tsv"
    if summary_file.is_file() and summary_file.stat().st_size > 100:
        log.info("[0/5] Cenote-Taker3 — 跳过 (已存在)")
        return ct3_out

    ct3_out.mkdir(parents=True, exist_ok=True)
    wt = Path(ct3_out) / "cenote_workdir"
    wt.mkdir(exist_ok=True)

    n_seqs = sum(1 for _ in open(fasta) if _.startswith('>'))
    log.info("[0/5] Cenote-Taker3: %d 病毒序列, DB=%s", n_seqs, cenote_dbs)

    cmd = (f"{cenote_bin} -c {fasta} -r plant_virus_analysis "
           f"-p False -t {threads} -am True "
           f"-wd {wt} --cenote-dbs {cenote_dbs} "
           f"--minimum_length_circular 500 --minimum_length_linear 500 "
           f"--molecule_type {args.molecule_type} --seqtech {args.seqtech} "
           f"--caller prodigal-gv --taxdb hallmark "
           f"--circ_minimum_hallmark_genes 0 --lin_minimum_hallmark_genes 1")
    if args.isolation_source:
        cmd += f" --isolation_source {args.isolation_source}"
    if args.collection_date:
        cmd += f" --collection_date {args.collection_date}"
    if args.assembler_info:
        cmd += f" --assembler {args.assembler_info}"
    ok = run(cmd, log, "Cenote-Taker3")
    if ok:
        # Collect output
        for f in wt.glob("*.tsv"):
            import shutil
            shutil.copy(f, ct3_out / f.name)
        for f in wt.glob("*.gbf"):
            import shutil
            shutil.copy(f, ct3_out / f.name)
    return ct3_out if ok else None


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
    p.add_argument("--skip-cenote", action="store_true", help="跳过 Cenote-Taker3")
    p.add_argument("--skip-suvtk", action="store_true", help="跳过 suvtk 步骤")
    # Cenote-Taker3 / GenBank 元数据
    p.add_argument("--molecule-type", default="RNA", choices=["DNA","RNA"], help="分子类型 (默认: RNA)")
    p.add_argument("--seqtech", default="Illumina", help="测序平台 (默认: Illumina)")
    p.add_argument("--isolation-source", help="样本地理来源")
    p.add_argument("--collection-date", help="采集日期 (DD-Mmm-YYYY)")
    p.add_argument("--assembler-info", help="组装工具及版本")
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

    # 0. Cenote-Taker3 (病毒 hallmark 基因检测)
    if not args.skip_cenote:
        cenote_out = run_cenote(inp, out, args.threads, log)
    else:
        log.info("[0/5] Cenote-Taker3 — 跳过")
        cenote_out = None

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

    elapsed = time.time() - t0
    log.info("=" * 60)
    log.info("完成! 耗时 %.0fs → %s", elapsed, out)
    log.info("  cenote_taker3/       病毒 hallmark 基因 + 功能注释 + 分类")
    log.info("  suvtk_taxonomy/       ICTV 分类注释")
    log.info("  suvtk_features/       基因注释 (CDS/tRNA)")
    log.info("  hypothetical/          假想蛋白深度注释")
    log.info("  submission/            GenBank 提交文件 (.sqn)")
    log.info("=" * 60)


if __name__ == "__main__":
    main()
