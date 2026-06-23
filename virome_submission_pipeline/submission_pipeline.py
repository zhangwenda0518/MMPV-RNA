#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
submission_pipeline.py — MMPV-RNA → NCBI 提交统一入口 v1.0
==========================================================

接收上游 MMPV-RNA 管道 (08_Rescue) 的输出, 自动检测已有文件,
按需运行 suvtk taxonomy/features/hypothetical,
最终生成两种格式的提交文件。

上游输入:
  08_Rescue/
  ├── all_plant_viruses.fasta         (必需)
  ├── suvtk.taxonomy_output/          (可选, 已存在则跳过)
  ├── suvtk.features_output/          (可选, 已存在则跳过)
  ├── analyze_hypothetical/           (可选, 已存在则跳过)
  └── metadata_output/                (可选, 自动填充 source.src)

两种输出模式:
  A) suvtk 模式: source.src + suvtk comments → table2asn → .sqn
  B) Sequin 模式: .fsa + .tbl + .cmt → tbl2asn → .sqn (Cenote-Taker3 风格)

用法:
  # 从 08_Rescue 一键运行
  python submission_pipeline.py \\
      --work-dir $OUT/08_Rescue/ \\
      --run-title ningxiagouqi_plant_virome \\
      --mode both \\
      --suvtk-db ~/database/virus-db/suvtk_db/ \\
      --rvdb-db ~/database/virus-db/RVDB-v31/ \\
      -t 40
"""

import argparse
import os
import re
import sys
import subprocess
import logging
from pathlib import Path
from datetime import datetime

SCRIPT_DIR = Path(__file__).resolve().parent

# ══════════════════════════════════════════════════════════════
# 管道参数自动检测 (从 pipeline_config.yaml)
# ══════════════════════════════════════════════════════════════

def load_pipeline_params(work_dir):
    """从 pipeline_config.yaml 自动提取组装/测序/富集参数"""
    import yaml as _yaml
    params = {
        'assembler': 'MEGAHIT;1.2.9;default parameters',
        'sequencer': 'Illumina NovaSeq 6000',
        'enrichment': 'none',
        'source_uvig': 'metatranscriptome (not viral targeted)',  # 默认公共转录组
        'metagenome_source': 'plant virome',
        'pipeline': 'MMPV-RNA v2.3 + suvtk v0.1.1',
    }

    # 从 work_dir 往上找 pipeline_config.yaml
    for parent in [work_dir, work_dir.parent, work_dir.parent.parent,
                   work_dir.parent.parent.parent]:
        config_path = parent / 'pipeline_config.yaml'
        if not config_path.exists():
            continue
        try:
            with open(config_path) as f:
                config = _yaml.safe_load(f)
            profiles = config.get('profiles', {})
            default = profiles.get('default', {})
            # 组装
            assembly = default.get('assembly', {})
            if assembly.get('assembler'):
                params['assembler'] = f"{assembly['assembler']};auto;default parameters"
            # 工具
            tools = default.get('tools', {})
            if tools.get('diamond'):
                pass  # 已全局使用
            # 运行时
            runtime = default.get('runtime', {})
            # 默认用 public data 的 source_uvig
            break
        except Exception:
            continue

    return params


# ══════════════════════════════════════════════════════════════
# 日志
# ══════════════════════════════════════════════════════════════

def setup_logger(out_dir):
    logger = logging.getLogger("SubPipeline")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    os.makedirs(out_dir, exist_ok=True)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(ch)

    fh = logging.FileHandler(os.path.join(out_dir, 'submission_pipeline.log'), encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s %(message)s'))
    logger.addHandler(fh)

    return logger


def run(cmd, log, step_name, check=False):
    """执行命令, 容忍失败"""
    log.info("[%s] %s", step_name, cmd[:120])
    try:
        subprocess.run(cmd, shell=True, check=check, executable='/bin/bash')
        log.info("[%s] ✓", step_name)
        return True
    except subprocess.CalledProcessError as e:
        log.error("[%s] ✗ (exit=%d)", step_name, e.returncode)
        return False


# ══════════════════════════════════════════════════════════════
# 自动检测 + 补运行缺失步骤
# ══════════════════════════════════════════════════════════════

def ensure_taxonomy(work_dir, fasta, suvtk_db, threads, log):
    """确保 suvtk taxonomy 已完成"""
    tax_out = work_dir / "suvtk.taxonomy_output"
    tax_tsv = tax_out / "taxonomy.tsv"
    if tax_tsv.exists() and os.path.getsize(tax_tsv) > 100:
        log.info("[1/5] taxonomy — 已存在, 跳过")
        return tax_out

    log.info("[1/5] 运行 suvtk taxonomy ...")
    run(
        f"suvtk taxonomy -i {fasta} -o {tax_out} -d {suvtk_db} -s 0.7 -t {threads}",
        log, "suvtk taxonomy"
    )
    return tax_out if tax_tsv.exists() else None


def ensure_features(work_dir, fasta, tax_dir, suvtk_db, threads, log):
    """确保 suvtk features 已完成"""
    feat_out = work_dir / "suvtk.features_output"
    tbl = feat_out / "featuretable.tbl"
    if tbl.exists() and os.path.getsize(tbl) > 100:
        log.info("[2/5] features — 已存在, 跳过")
        return feat_out

    tax_tsv = tax_dir / "taxonomy.tsv" if tax_dir else None
    cmd = f"suvtk features -i {fasta} -o {feat_out} -d {suvtk_db} --coding-complete -t {threads}"
    if tax_tsv and tax_tsv.exists():
        cmd += f" --taxonomy {tax_tsv}"

    log.info("[2/5] 运行 suvtk features ...")
    run(cmd, log, "suvtk features")
    return feat_out if tbl.exists() else None


def ensure_hypothetical(work_dir, feat_dir, rvdb_dir, threads, log):
    """确保 analyze_hypothetical 已完成"""
    hypo_out = work_dir / "analyze_hypothetical"
    updated_tbl = hypo_out / "featuretable_updated.tbl"
    if updated_tbl.exists() and os.path.getsize(updated_tbl) > 100:
        log.info("[2.5/5] hypothetical — 已存在, 跳过")
        return hypo_out

    tbl = feat_dir / "featuretable.tbl"
    faa = feat_dir / "proteins.faa"
    if not tbl.exists() or not faa.exists():
        log.warning("[2.5/5] hypothetical — 缺少 input, 跳过")
        return None

    # 自动检测 RVDB 数据库
    rvdb_fasta = Path(rvdb_dir) / "U-RVDBv31.0-prot.EX.acc.fasta" if rvdb_dir else None
    rvdb_hmm = Path(rvdb_dir) / "U-RVDBv31.0-prot.hmm" if rvdb_dir else None
    rvdb_annot = Path(rvdb_dir) / "U-RVDBv31.0-prot.info.tab" if rvdb_dir else None

    script = SCRIPT_DIR / "analyze_hypothetical.py"

    cmd = (
        f"python {script} "
        f"-t {tbl} -f {faa} "
        f"-o {hypo_out} "
        f"--blast {hypo_out}/merged_blast.txt "
        f"--diamond -d {rvdb_fasta} "
        f"--hmmer --hmmer-db {rvdb_hmm} "
        f"--annot {rvdb_annot} "
        f"--threads {threads}"
    )

    log.info("[2.5/5] 运行 analyze_hypothetical (diamond + HMM) ...")
    run(cmd, log, "analyze_hypothetical")
    return hypo_out if updated_tbl.exists() else None


# ══════════════════════════════════════════════════════════════
# 已知病毒管道 (INFO → suvtk)
# ══════════════════════════════════════════════════════════════

def gather_metadata_from_summary(all_summary_path, output_dir, metadata_dir=None, log=None):
    """INFO 步骤: 从 all_summary.tsv 提取所有 SRA/CRR, 运行 gsa_sra.info.py 生成元数据。

    流程:
      1. 从 all_summary.tsv 提取唯一 SRA/CRR 编号
      2. 运行 gsa_sra.info.py -i <id_list> -o <output_dir>/0_info/
         → 生成 Global_Unified_Metadata_Full.tsv / Core13.tsv
      3. 合并额外的 Corrected_Metadata.tsv 和 runinfo.csv (如存在)
      4. 返回合并后的元数据 lookup

    Returns: dict {SRA_run: {collection_date, geo_loc_name, bioproject,
                              biosample, platform, host, tissue, source, center}}
    """
    import pandas as pd

    # 1. 提取唯一 SRA/CRR
    summary = pd.read_csv(all_summary_path, sep='\t')
    all_runs = set()
    for _, row in summary.iterrows():
        sra = str(row.get('Sample', '')).strip()
        if sra and sra.lower() not in ('nan', ''):
            all_runs.add(sra)

    if log:
        log.info("[INFO] 从 all_summary 提取 %d 个唯一 SRA/CRR", len(all_runs))

    # 2. 运行 gsa_sra.info.py 生成元数据
    info_out = Path(output_dir) / '0_info'
    info_out.mkdir(parents=True, exist_ok=True)
    full_tsv = info_out / 'Global_Unified_Metadata_Full.tsv'

    if not (full_tsv.exists() and os.path.getsize(full_tsv) > 1000):
        id_file = info_out / 'run_ids.txt'
        with open(id_file, 'w') as f:
            for rid in sorted(all_runs):
                f.write(rid + '\n')

        gsa_script = SCRIPT_DIR.parent / 'public_metadata_pipeline' / 'gsa_sra.info.py'
        if not gsa_script.exists():
            gsa_script = SCRIPT_DIR.parent.parent / 'public_metadata_pipeline' / 'gsa_sra.info.py'

        if gsa_script.exists():
            if log:
                log.info("  运行 gsa_sra.info.py — 获取 %d 个 Run 的元数据 ...", len(all_runs))
            cmd = f'python {gsa_script} -i {id_file} -o {info_out} -m local'
            run(cmd, log, 'gsa_sra.info')
        else:
            log.warning("  未找到 gsa_sra.info.py, 跳过在线元数据获取")

    # 3. 读取生成的元数据
    meta_lookup = {}

    def _ingest(df, source_name, log):
        added = 0
        for _, row in df.iterrows():
            run = str(row.get('Run', '')).strip()
            if not run or run not in all_runs:
                continue
            if run not in meta_lookup:
                meta_lookup[run] = {}

            for out_key, col_name in [
                ('collection_date', 'CollectionDate'),
                ('geo_loc_name', 'Location'),
                ('bioproject', 'BioProject'),
                ('biosample', 'BioSample'),
                ('platform', 'Platform'),
                ('host', 'ScientificName'),
                ('tissue', 'Tissue'),
                ('source', 'Source'),
                ('center', 'CenterName'),
            ]:
                if col_name in df.columns:
                    v = str(row.get(col_name, '')).strip()
                    if v and v.lower() not in ('nan', 'not_provided', '', ' '):
                        meta_lookup[run][out_key] = v
                        added += 1
        if log:
            log.info("  %s: 贡献 %d 字段", source_name, added)

    # Load Full.tsv (from gsa_sra.info.py output)
    if full_tsv.exists():
        df_full = pd.read_csv(full_tsv, sep=None, engine='python')
        df_full.columns = [c.lstrip('\ufeff') for c in df_full.columns]
        _ingest(df_full, 'Full', log)

    # Core13
    core13 = info_out / 'Global_Unified_Metadata_Core13.tsv'
    if core13.exists():
        df_core = pd.read_csv(core13, sep=None, engine='python')
        df_core.columns = [c.lstrip('\ufeff') for c in df_core.columns]
        _ingest(df_core, 'Core13', log)

    # 4. 额外合并 Corrected_Metadata.tsv 和 runinfo.csv (优先覆盖)
    if metadata_dir and os.path.isdir(metadata_dir):
        corrected = os.path.join(metadata_dir, 'Corrected_Metadata.tsv')
        if os.path.exists(corrected):
            df_cm = pd.read_csv(corrected, sep='\t')
            _ingest(df_cm, 'Corrected_Metadata', log)

        runinfo = os.path.join(metadata_dir, 'runinfo.csv')
        if os.path.exists(runinfo):
            df_ri = pd.read_csv(runinfo, sep=',')
            _ingest(df_ri, 'runinfo', log)

    if log:
        n_with_date = sum(1 for v in meta_lookup.values() if v.get('collection_date'))
        n_with_bs = sum(1 for v in meta_lookup.values() if v.get('biosample'))
        n_with_plat = sum(1 for v in meta_lookup.values() if v.get('platform'))
        log.info("[INFO] 元数据合并: %d/%d 个 Run, 含日期=%d, BioSample=%d, Platform=%d",
                 len(meta_lookup), len(all_runs), n_with_date, n_with_bs, n_with_plat)

    return meta_lookup


def build_taxonomy_from_summary(all_summary_path, fasta_dir=None, log=None):
    """从 all_summary.tsv 构建 contig-level taxonomy, 按 FASTA 过滤。

    Returns:
        tax_tsv_path, n_seqs, n_viruses
    """
    import pandas as pd

    df = pd.read_csv(all_summary_path, sep='\t')
    contig_tax = {}
    for _, row in df.iterrows():
        sample = str(row.get('Sample', '')).strip()
        acc = str(row.get('Accession', '')).strip()
        taxonomy = str(row.get('Taxonomy', '')).strip()
        if not sample or not acc:
            continue
        contig = f'{sample}_{acc}'
        if taxonomy and contig not in contig_tax:
            contig_tax[contig] = taxonomy

    # Filter by FASTA
    if fasta_dir and os.path.isdir(fasta_dir):
        fa_contigs = set()
        from pathlib import Path
        for fpath in Path(fasta_dir).rglob('*.fasta'):
            for line in open(fpath):
                if line.startswith('>'):
                    fa_contigs.add(line[1:].strip())
        before = len(contig_tax)
        contig_tax = {k: v for k, v in contig_tax.items() if k in fa_contigs}
        if log:
            log.info("[INFO] FASTA 过滤: %d → %d contigs (匹配 %d FASTA 头)",
                     before, len(contig_tax), len(fa_contigs))

    # Write taxonomy.tsv
    out_dir = Path(all_summary_path).parent.parent  # go up from summary/
    tax_path = out_dir / 'taxonomy.tsv'
    with open(tax_path, 'w') as f:
        f.write('contig\ttaxonomy\n')
        for contig, tax in sorted(contig_tax.items()):
            f.write(f'{contig}\t{tax}\n')

    return str(tax_path), len(contig_tax)


def prepare_known_viruses(all_summary_path, fasta_dir, metadata_dir, output_dir,
                          run_title, suvtk_db, threads, sequencer, assembler,
                          authors, pipeline, enrichment, source_uvig, metagenome_source,
                          log):
    """已知病毒全流程: INFO → suvtk taxonomy → suvtk features → unified_metadata → 导出。

    参数:
      all_summary_path : 2_Virus_variants_Results/summary/all_summary.tsv
      fasta_dir        : 4_assemblies_clean/
      metadata_dir     : 6.Virus_variants_Results/summary/ (含 Full.tsv, runinfo.csv 等)
      output_dir       : 提交输出目录
    """
    import pandas as pd
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    log.info("=" * 60)
    log.info("已知病毒 GenBank 提交准备 (Known Virus Pipeline)")
    log.info("  all_summary: %s", all_summary_path)
    log.info("  fasta_dir:   %s", fasta_dir)
    log.info("  output:      %s", out)
    log.info("=" * 60)

    # ── Step 0: INFO — gather metadata ──
    log.info("[0/4] INFO — 收集元数据 ...")
    meta_lookup = gather_metadata_from_summary(
        all_summary_path, str(out), metadata_dir, log)

    # ── Step 1: Build taxonomy + combine FASTA ──
    log.info("[1/4] 构建 taxonomy + 合并 FASTA ...")
    tax_path, n_contigs = build_taxonomy_from_summary(all_summary_path, fasta_dir, log)

    # Combine FASTA
    combined_fa = out / 'combined_known_viruses.fasta'
    with open(combined_fa, 'w') as cf:
        for fpath in Path(fasta_dir).rglob('*.fasta'):
            with open(fpath) as inf:
                cf.write(inf.read())
    n_fa = sum(1 for _ in open(combined_fa) if _.startswith('>'))
    log.info("  合并 FASTA: %d 条序列 → %s", n_fa, combined_fa)

    # ── Step 2: suvtk taxonomy ──
    if suvtk_db:
        log.info("[2/4] suvtk taxonomy ...")
        tax_dir = out / '1_taxonomy'
        tax_dir.mkdir(exist_ok=True)
        tax_out_tsv = tax_dir / 'taxonomy.tsv'
        if not tax_out_tsv.exists() or os.path.getsize(tax_out_tsv) < 100:
            cmd = (f'suvtk taxonomy -i {combined_fa} -o {tax_dir} '
                   f'-d {suvtk_db} -s 0.7 -t {threads}')
            run(cmd, log, 'suvtk taxonomy')
        else:
            log.info("  taxonomy 已存在, 跳过")
    else:
        log.info("[2/4] suvtk taxonomy — 跳过 (未指定 --suvtk-db)")

    # ── Step 3: suvtk features ──
    if suvtk_db:
        log.info("[3/4] suvtk features ...")
        feat_dir = out / '2_features'
        feat_dir.mkdir(exist_ok=True)
        tbl = feat_dir / 'featuretable.tbl'
        if not tbl.exists() or os.path.getsize(tbl) < 100:
            tax_tsv = out / '1_taxonomy' / 'taxonomy.tsv'
            cmd = (f'suvtk features -i {combined_fa} -o {feat_dir} '
                   f'-d {suvtk_db} --coding-complete -t {threads}')
            if tax_tsv.exists():
                cmd += f' --taxonomy {tax_tsv}'
            run(cmd, log, 'suvtk features')
        else:
            log.info("  features 已存在, 跳过")
    else:
        log.info("[3/4] suvtk features — 跳过 (未指定 --suvtk-db)")

    # ── Step 4: unified_metadata + 导出 ──
    log.info("[4/4] 生成统一元数据 + 导出 ...")
    meta_dir = out / 'submission'
    meta_dir.mkdir(exist_ok=True)

    meta_script = SCRIPT_DIR / 'unified_metadata.py'

    # Use the contig-taxonomy from summary (not suvtk's taxonomy)
    # suvtk taxonomy.tsv might have different format
    summary_tax = tax_path  # from build_taxonomy_from_summary

    cmd = (
        f'python {meta_script} '
        f'--taxonomy {summary_tax} '
        f'--run-title {run_title} '
        f'--assembler "{assembler}" '
        f'--sequencer "{sequencer}" '
        f'--enrichment "{enrichment}" '
        f'--source-uvig "{source_uvig}" '
        f'--metagenome-source "{metagenome_source}" '
        f'--pipeline "{pipeline}" '
        f'--authors "{authors}" '
        f'-o {meta_dir}/'
    )
    if metadata_dir:
        # Point to metadata dir so load_metadata can auto-detect Full.tsv
        pass
    run(cmd, log, 'unified_metadata')

    # ── Enrich CSV with metadata + abundance ──
    csv_path = meta_dir / 'unified_metadata.csv'
    if csv_path.exists():
        log.info("  充实 CSV 元数据 ...")
        df = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
        enriched = 0

        # Build abundance lookup from all_summary.tsv (contig → {MeanDepth, Covered%, ...})
        abu_lookup = {}
        from_summary = pd.read_csv(all_summary_path, sep='\t')
        for _, r2 in from_summary.iterrows():
            s = str(r2.get('Sample', '')).strip()
            a = str(r2.get('Accession', '')).strip()
            if s and a:
                contig_key = f'{s}_{a}'
                abu_lookup[contig_key] = {
                    'depth': str(r2.get('Rep_MeanDepth', '')).strip(),
                    'covered': str(r2.get('Covered%', '')).strip(),
                    'reads': str(r2.get('Reads', '')).strip(),
                    'rel_abund': str(r2.get('Asm_Rel_Abund(%)', '')).strip(),
                }

        for idx, row in df.iterrows():
            contig = str(row.get('sequence_name', '')).strip()
            sra = re.match(r'([SC]RR\d+)', contig)
            sra = sra.group(1) if sra else contig.split('_')[0] if '_' in contig else contig

            # ── Fill from gsa_sra.info.py metadata ──
            if sra in meta_lookup:
                m = meta_lookup[sra]
                for csv_col, meta_key in [
                    ('collection_date', 'collection_date'),
                    ('src-geo_loc_name', 'geo_loc_name'),
                    ('bioproject', 'bioproject'),
                    ('biosample', 'biosample'),
                    ('src-Host', 'host'),
                    ('src-Tissue_type', 'tissue'),
                    ('src-Cultivar', 'source'),
                    ('src-Collected_by', 'center'),
                ]:
                    if csv_col in df.columns and m.get(meta_key):
                        current = str(row.get(csv_col, '')).strip()
                        is_ph = (not current or current.lower() in ('nan', 'not_provided', '') or
                                 'XXXX' in current or 'YYYY' in current or
                                 'Country:Region' in current or 'PRJNAXXXX' in current or
                                 'SAMNXXXXXXXX' in current or 'XX.' in current)
                        if is_ph:
                            df.at[idx, csv_col] = m[meta_key]
                            enriched += 1
                # Platform → Sequencing Technology
                if m.get('platform') and 'cmt-Sequencing_Technology' in df.columns:
                    cur = str(row.get('cmt-Sequencing_Technology', '')).strip()
                    if cur == sequencer or not cur or 'Illumina NovaSeq' in cur:
                        df.at[idx, 'cmt-Sequencing_Technology'] = m['platform']
                        enriched += 1

            # ── Fill abundance from all_summary ──
            if contig in abu_lookup and 'cmt-Genome_Coverage' in df.columns:
                a = abu_lookup[contig]
                cur_cov = str(row.get('cmt-Genome_Coverage', '')).strip()
                if not cur_cov:
                    depth = a['depth']
                    covered = a['covered']
                    parts = []
                    if depth and depth not in ('nan', ''):
                        parts.append(f'{depth}x')
                    if covered and covered not in ('nan', ''):
                        parts.append(f'{covered}%')
                    if parts:
                        df.at[idx, 'cmt-Genome_Coverage'] = ' '.join(parts)
                        enriched += 1

        df.to_csv(csv_path, index=False, encoding='utf-8-sig')
        log.info("  充实 %d 个单元格", enriched)

    # ── Topology ──
    fna = feat_dir / 'reoriented_nucleotide_sequences.fna'
    if fna.exists():
        topo_script = SCRIPT_DIR / 'viral_topology.py'
        topo_out = meta_dir / 'topology.tsv'
        run(f'python {topo_script} --taxonomy {summary_tax} --fasta {fna} -o {topo_out}',
            log, 'viral_topology')

    # ── HTML report ──
    if csv_path.exists():
        report_script = SCRIPT_DIR / 'report_html.py'
        report_out = meta_dir / 'report.html'
        run(f'python {report_script} --csv {csv_path} --run-title {run_title} -o {report_out}',
            log, 'report_html')

    log.info("=" * 60)
    log.info("完成! 输出: %s", meta_dir)
    log.info("  序列: %d, 元数据覆盖: %d/%d", n_contigs, len(meta_lookup), n_contigs)
    log.info("=" * 60)

    return str(meta_dir)


# ══════════════════════════════════════════════════════════════
# 主流程
# ══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="submission_pipeline.py — MMPV-RNA → NCBI 提交统一入口",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest='command', help='管道模式')

    # ── known 模式: 已知病毒 (INFO → suvtk) ──
    p_known = sub.add_parser('known', help='已知病毒提交 (INFO → suvtk)',
        epilog="""
示例:
  python submission_pipeline.py known \\
      --all-summary known_virus_pipeline/2_Virus_variants_Results/summary/all_summary.tsv \\
      --fasta-dir known_virus_pipeline/4_assemblies_clean/ \\
      --metadata-dir 6.Virus_variants_Results/summary/ \\
      --suvtk-db ~/database/virus-db/suvtk_db/ \\
      -o genbank_submission/test_known/ \\
      -t 40
""")
    p_known.add_argument('--all-summary', required=True,
                         help='all_summary.tsv 路径 (2_Virus_variants_Results/summary/)')
    p_known.add_argument('--fasta-dir', required=True,
                         help='4_assemblies_clean/ 目录')
    p_known.add_argument('--metadata-dir',
                         help='元数据目录 (默认自动从 all_summary 推断: 同 data-* 父级下的 6.Virus_variants_Results/summary/)')
    p_known.add_argument('--suvtk-db', help='suvtk 数据库路径')
    p_known.add_argument('-o', '--output', default='./genbank_submission/known',
                         help='输出目录')
    p_known.add_argument('--run-title', default='known_virus_submission',
                         help='运行标题')
    p_known.add_argument('-t', '--threads', type=int, default=40)
    p_known.add_argument('--assembler', default='MEGAHIT;1.2.9;default parameters')
    p_known.add_argument('--sequencer', default='Illumina NovaSeq 6000')
    p_known.add_argument('--enrichment', default='none')
    p_known.add_argument('--pipeline', default='MMPV-RNA v2.3 + suvtk v0.1.1')
    p_known.add_argument('--source-uvig', default='metatranscriptome (not viral targeted)')
    p_known.add_argument('--metagenome-source', default='plant virome')
    p_known.add_argument('--authors', default='Zhang, Wenda')

    # ── 默认模式: discovery pipeline (08_Rescue) ──
    p_disc = sub.add_parser('discovery', help='新病毒提交 (discovery pipeline, 08_Rescue)',
        epilog="""
示例:
  python submission_pipeline.py discovery --work-dir $OUT/08_Rescue/ --suvtk-db ~/.../suvtk_db/ -t 40
""")
    p_disc.add_argument('--work-dir', required=True,
                        help='上游输出目录 (含 all_plant_viruses.fasta)')
    p_disc.add_argument('--fasta',
                        help='输入 FASTA (默认 work-dir/all_plant_viruses.fasta)')
    p_disc.add_argument('--suvtk-db', help='suvtk 数据库路径')
    p_disc.add_argument('--rvdb-db', help='RVDB 数据库目录')
    p_disc.add_argument('--run-title', default='viral_submission', help='运行标题')
    p_disc.add_argument('--mode', choices=['suvtk', 'sequin', 'both'], default='both')
    p_disc.add_argument('-t', '--threads', type=int, default=40)
    p_disc.add_argument('--fetch-metadata', action='store_true')
    p_disc.add_argument('--metadata-dir', help='Core13 元数据目录')
    p_disc.add_argument('--config', help='YAML 配置文件')
    p_disc.add_argument('--interactive', action='store_true')
    p_disc.add_argument('--pipeline-type', choices=['auto', 'discovery', 'analysis'],
                        default='auto')
    p_disc.add_argument('--ref-info', help='已知病毒参考信息')
    p_disc.add_argument('--force-taxonomy', action='store_true')
    p_disc.add_argument('--force-features', action='store_true')
    p_disc.add_argument('--force-hypo', action='store_true')
    p_disc.add_argument('--export-metadata', action='store_true')

    args = parser.parse_args()

    # ── known 模式 ──
    if args.command == 'known':
        out = Path(args.output)
        log = setup_logger(str(out))

        # Auto-detect metadata_dir from all_summary path
        metadata_dir = args.metadata_dir
        if not metadata_dir:
            summary_dir = Path(args.all_summary).resolve().parent  # .../2_Virus_variants_Results/summary
            # Go up 3 levels to data-{year}/sra_rna.dataN/
            for up_to in [summary_dir.parent.parent.parent,  # data-*/sra_rna.*
                          summary_dir.parent.parent.parent.parent]:
                candidate = up_to / '6.Virus_variants_Results' / 'summary'
                if candidate.is_dir():
                    metadata_dir = str(candidate)
                    log.info("  自动检测 metadata-dir: %s", metadata_dir)
                    break
            if not metadata_dir:
                log.warning("  未找到 6.Virus_variants_Results/summary/, 元数据将为空")

        prepare_known_viruses(
            all_summary_path=args.all_summary,
            fasta_dir=args.fasta_dir,
            metadata_dir=metadata_dir or '',
            output_dir=str(out),
            run_title=args.run_title,
            suvtk_db=args.suvtk_db or '',
            threads=args.threads,
            sequencer=args.sequencer,
            assembler=args.assembler,
            authors=args.authors,
            pipeline=args.pipeline,
            enrichment=args.enrichment,
            source_uvig=args.source_uvig,
            metagenome_source=args.metagenome_source,
            log=log,
        )
        return

    # ── discovery 模式 ──
    if args.command != 'discovery':
        parser.print_help()
        sys.exit(0)

    work_dir = Path(args.work_dir)
    fasta = args.fasta or str(work_dir / "all_plant_viruses.fasta")

    if not os.path.exists(fasta):
        # 尝试自动找
        candidates = list(work_dir.glob("*.fasta")) + list(work_dir.glob("*.fna"))
        if candidates:
            fasta = str(candidates[0])
            print(f"[*] 自动检测 FASTA: {fasta}")
        else:
            print(f"[!] 未找到 FASTA 文件: {fasta}")
            sys.exit(1)

    log = setup_logger(str(work_dir))

    # ── 自动检测管道类型 ──
    pipeline_type = args.pipeline_type
    if pipeline_type == 'auto':
        # 检测 discovery pipeline 特征 (Plant/, centroids/, all_plant_viruses.fasta)
        is_discovery = (work_dir / "Plant").exists() or (work_dir / "centroids").exists()
        # 检测 analysis pipeline 特征 (4_assemblies_clean/, known_viruses/)
        is_analysis = (work_dir.parent / "4_assemblies_clean").exists() or \
                      (work_dir / "4_assemblies_clean").exists() or \
                      (len(list(work_dir.glob("*_assemblies_clean"))) > 0)

        if is_discovery and not is_analysis:
            pipeline_type = 'discovery'
        elif is_analysis and not is_discovery:
            pipeline_type = 'analysis'
        elif is_discovery and is_analysis:
            pipeline_type = 'discovery'
        else:
            pipeline_type = 'discovery'

    # ── 收集病毒列表 ──
    virus_list = []  # [{name, accession, fastas: [path, ...], n_samples}]
    if pipeline_type == 'analysis':
        # 检测 analysis pipeline 输出: 4_assemblies_clean/
        analysis_dir = work_dir.parent / "4_assemblies_clean"
        if not analysis_dir.exists():
            analysis_dir = work_dir / "4_assemblies_clean"
        if not analysis_dir.exists():
            # 尝试从 work_dir 的兄弟目录找
            for parent in [work_dir.parent, work_dir.parent.parent]:
                candidate = parent / "4_assemblies_clean"
                if candidate.exists():
                    analysis_dir = candidate
                    break

        if analysis_dir.exists():
            for virus_dir in sorted(analysis_dir.iterdir()):
                if virus_dir.is_dir() and not virus_dir.name.startswith('.'):
                    fastas = sorted(virus_dir.glob("*.full.fasta")) + sorted(virus_dir.glob("*.fasta"))
                    # 从目录名解析病毒名和 accession
                    parts = virus_dir.name.rsplit('_', 1)
                    virus_name = parts[0] if len(parts) > 1 else virus_dir.name
                    accession = parts[1] if len(parts) > 1 else ''
                    virus_list.append({
                        'name': virus_name,
                        'accession': accession,
                        'dir': str(virus_dir),
                        'fastas': [str(f) for f in fastas],
                        'n_samples': len(fastas),
                    })
            log.info("  检测到 analysis pipeline: 1 种病毒 × %d 个样本 (全长序列)", len(virus_list) if virus_list else 0)
            for v in virus_list:
                log.info("    %s (%s): %d samples", v['name'], v['accession'], v['n_samples'])
        else:
            log.warning("  未找到 4_assemblies_clean/ 目录")

    elif pipeline_type == 'discovery':
        log.info("  检测到 discovery pipeline: %s", fasta)

    log.info("=" * 60)
    log.info("MMPV-RNA → NCBI 提交流水线 v1.0")
    log.info("  工作目录: %s", work_dir)
    log.info("  管道类型: %s — %s",
             pipeline_type,
             f"1种病毒 × {len(virus_list[0]['n_samples']) if virus_list else '?'} 个样本" if pipeline_type == 'analysis'
             else "多病毒, 每病毒1条代表序列")
    log.info("  输出模式: %s", args.mode)
    log.info("=" * 60)

    # ═══ 步骤 1-2.5: 自动补缺失 ═══

    # 1. taxonomy
    tax_dir = None
    if args.suvtk_db:
        tax_dir = work_dir / "suvtk.taxonomy_output"
        if args.force_taxonomy and tax_dir.exists():
            import shutil
            shutil.rmtree(tax_dir)
        tax_dir = ensure_taxonomy(work_dir, fasta, args.suvtk_db, args.threads, log)

    # 2. features
    feat_dir = None
    if args.suvtk_db:
        feat_dir = work_dir / "suvtk.features_output"
        if args.force_features and feat_dir.exists():
            import shutil
            shutil.rmtree(feat_dir)
        feat_dir = ensure_features(work_dir, fasta, tax_dir, args.suvtk_db, args.threads, log)

    # 2.5. hypothetical (仅 discovery pipeline 需要)
    hypo_dir = None
    if pipeline_type == 'discovery' and feat_dir and args.rvdb_db:
        hypo_dir = work_dir / "analyze_hypothetical"
        if args.force_hypo and hypo_dir.exists():
            import shutil
            shutil.rmtree(hypo_dir)
        hypo_dir = ensure_hypothetical(work_dir, feat_dir, args.rvdb_db, args.threads, log)
    elif pipeline_type == 'analysis':
        log.info("[2.5/5] hypothetical — 已知病毒跳过 (analysis pipeline 使用参考基因组注释)")

    # 2.6. 公共元数据获取 (public_metadata_pipeline)
    meta_dir = args.metadata_dir or str(work_dir / "metadata_output")
    core13 = Path(meta_dir) / "Global_Unified_Metadata_Full.tsv"
    if not core13.exists():
        core13 = Path(meta_dir) / "Global_Unified_Metadata_Core13.tsv"

    if args.fetch_metadata and tax_tsv and tax_tsv.exists():
        log.info("[2.6/5] 获取公共元数据 (NCBI SRA → 自动填 source.src)")

        # 提取 SRA/CRR 编号
        srr_list = work_dir / "srr_list.txt"
        if not srr_list.exists() or args.force_hypo:
            log.info("  提取 SRA/CRR 编号 ...")
            import subprocess as _sp
            result = _sp.run(
                f"grep -oP '[SC]RR\\d+' {tax_tsv} | sort -u > {srr_list}",
                shell=True, executable='/bin/bash'
            )
            n_sra = sum(1 for _ in open(srr_list)) if srr_list.exists() else 0
            log.info("  找到 %d 个唯一 SRA/CRR", n_sra)

        if not core13.exists():
            log.info("  运行 gsa_sra.info.py 获取元数据 ...")
            gsa_script = work_dir.parent.parent.parent / "MMPV-RNA" / "public_metadata_pipeline" / "gsa_sra.info.py"
            if not gsa_script.exists():
                gsa_script = Path.home() / "MMPV-RNA" / "public_metadata_pipeline" / "gsa_sra.info.py"
            cmd = (
                f"python {gsa_script} "
                f"-i {srr_list} "
                f"-o {meta_dir} "
                f"-m local "
                f"-t {args.threads}"
            )
            run(cmd, log, "gsa_sra.info", check=False)

    if core13.exists():
        log.info("  元数据已就绪: %s", core13)
    else:
        log.info("  无公共元数据文件")
        if pipeline_type == 'analysis':
            log.info("  ═══════════════════════════════════════════════")
            log.info("  📋 新测序数据提交指南:")
            log.info("  1. unified_metadata.csv 已生成 (占位符待填)")
            log.info("  2. 去 NCBI 注册 BioProject → 拿到 PRJNA-xxx")
            log.info("     https://submit.ncbi.nlm.nih.gov/")
            log.info("  3. 用 biosample_template.tsv 批量注册 BioSample")
            log.info("  4. 编辑 unified_metadata.csv 填入:")
            log.info("     collection_date, geo_loc_name, lat_lon,")
            log.info("     bioproject (你的 PRJNA), biosample (你的 SAMN)")
            log.info("  5. 重新运行 --export-metadata 不带 --fetch-metadata")
            log.info("  ═══════════════════════════════════════════════")
        else:
            log.info("  公共数据: source.src 从 Core13 自动填, 样本信息已有")
        meta_dir = None

    # ═══ 步骤 3-5: 生成提交文件 ═══

    # 确定使用的 featuretable
    updated_tbl = None
    updated_faa = None
    if hypo_dir and (hypo_dir / "featuretable_updated.tbl").exists():
        updated_tbl = str(hypo_dir / "featuretable_updated.tbl")
        updated_faa = str(hypo_dir / "proteins_updated.faa")
    elif feat_dir:
        updated_tbl = str(feat_dir / "featuretable.tbl")
        updated_faa = str(feat_dir / "proteins.faa")

    tax_tsv = tax_dir / "taxonomy.tsv" if tax_dir else None
    miuvig_tax = tax_dir / "miuvig_taxonomy.tsv" if tax_dir else None
    fna = feat_dir / "reoriented_nucleotide_sequences.fna" if feat_dir else None
    miuvig_feat = feat_dir / "miuvig_features.tsv" if feat_dir else None

    if args.mode in ('suvtk', 'both'):
        log.info("─" * 40 + "\n  suvtk 模式: 生成 source.src + .cmt → table2asn\n" + "─" * 40)

        suvtk_script = SCRIPT_DIR / "virome_submission.py"

        if tax_tsv and feat_dir:
            cmd = (
                f"python {suvtk_script} report "
                f"--work-dir {work_dir} "
                f"--updated-tbl {updated_tbl} "
                f"--updated-faa {updated_faa}"
            )
            run(cmd, log, "suvtk report")

        log.info("  下一步: 编辑 source.src 后运行 suvtk comments & table2asn")

    if args.mode in ('sequin', 'both'):
        log.info("─" * 40 + "\n  Sequin 模式: 生成 .fsa + .tbl + .cmt (Cenote-Taker3 风格)\n" + "─" * 40)

        sequin_script = SCRIPT_DIR / "sequin_builder.py"

        if tax_tsv and miuvig_tax and fna and updated_tbl and miuvig_feat:
            cmd = (
                f"python {sequin_script} "
                f"--taxonomy {tax_tsv} "
                f"--miuvig-tax {miuvig_tax} "
                f"--features {updated_tbl} "
                f"--proteins {updated_faa} "
                f"--fasta {fna} "
                f"--miuvig-feat {miuvig_feat} "
                f"--run-title {args.run_title} "
                f"-o {sub_out}/"
            )
            run(cmd, log, "sequin_builder")

    # ── 统一元数据 CSV 导出 (SeqSender 兼容) ──
    if args.export_metadata and tax_tsv:
        log.info("─" * 40 + "\n  导出统一元数据 CSV (SeqSender 兼容)\n" + "─" * 40)

        pipe_params = load_pipeline_params(work_dir)

        # 根据管道类型自动决定 source_uvig 和 enrichment
        if pipeline_type == 'discovery':
            params = pipe_params.copy()
        else:
            # analysis pipeline: 自己测序, 可能做了病毒富集
            params = pipe_params.copy()

        # 输出到 work_dir 平级，不污染分析目录
        sub_out = work_dir.parent / f"submission_{args.run_title}"
        sub_out.mkdir(parents=True, exist_ok=True)

        meta_script = SCRIPT_DIR / "unified_metadata.py"
        cmd = (
            f"python {meta_script} "
            f"--taxonomy {tax_tsv} "
            f"--run-title {args.run_title} "
            f"--assembler \"{params['assembler']}\" "
            f"--sequencer \"{params['sequencer']}\" "
            f"--enrichment \"{params['enrichment']}\" "
            f"--source-uvig \"{params['source_uvig']}\" "
            f"--metagenome-source \"{params['metagenome_source']}\" "
            f"--pipeline \"{params['pipeline']}\" "
            f"-o {sub_out}/"
        )

        if core13.exists():
            cmd += f" --reference-metadata {core13}"
        if args.config:
            cmd += f" --config {args.config}"
        if args.interactive:
            cmd += " --interactive"

        run(cmd, log, "unified_metadata")

        # ── 病毒拓扑判断 ──
        if fna and tax_tsv:
            topo_script = SCRIPT_DIR / "viral_topology.py"
            topo_out = sub_out / "topology.tsv"
            cmd_topo = (
                f"python {topo_script} --taxonomy {tax_tsv} --fasta {fna} -o {topo_out}"
            )
            run(cmd_topo, log, "viral_topology")

        # ── 交互式 HTML 报告 ──
        csv_path = sub_out / "unified_metadata.csv"
        miuvig_path = sub_out / "miuvig.tsv"
        asm_path = sub_out / "assembly.tsv"
        report_path = sub_out / "report.html"
        report_script = SCRIPT_DIR / "report_html.py"
        if csv_path.exists():
            cmd_report = (
                f"python {report_script} --csv {csv_path} "
                f"--run-title {args.run_title} -o {report_path}"
            )
            if miuvig_path.exists():
                cmd_report += f" --miuvig {miuvig_path}"
            if asm_path.exists():
                cmd_report += f" --assembly {asm_path}"
            run(cmd_report, log, "report_html")

    # ═══ 最终报告 ═══
    log.info("=" * 60)
    log.info("完成!")
    log.info("  分析目录: %s", work_dir)
    log.info("  提交输出: %s", sub_out)
    report_path = work_dir / "submission_metadata" / "report.html"
    if report_path.exists():
        log.info("  📊 交互报告:  %s", report_path)
    if tax_dir:
        log.info("  [1] taxonomy:    %s", tax_dir)
    if feat_dir:
        log.info("  [2] features:    %s", feat_dir)
    if hypo_dir:
        log.info("  [2.5] hypo:      %s", hypo_dir)
    if args.mode in ('suvtk', 'both'):
        log.info("  [suvtk]  report:  %s/report/", work_dir)
    if args.mode in ('sequin', 'both'):
        gb_dir = work_dir / "genbank_submission" / args.run_title
        log.info("  [sequin] output:  %s/", gb_dir)
    log.info("=" * 60)


if __name__ == '__main__':
    main()
