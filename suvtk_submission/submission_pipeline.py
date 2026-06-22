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
import sys
import subprocess
import logging
from pathlib import Path
from datetime import datetime

SCRIPT_DIR = Path(__file__).resolve().parent


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
# 主流程
# ══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="submission_pipeline.py — MMPV-RNA → NCBI 提交统一入口",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 从 08_Rescue 一键运行 (自动补缺失步骤)
  python submission_pipeline.py \\
      --work-dir $OUT/08_Rescue/ \\
      --run-title my_project \\
      --mode both \\
      --suvtk-db ~/database/virus-db/suvtk_db/ \\
      --rvdb-db ~/database/virus-db/RVDB-v31/ \\
      -t 40

  # 仅生成 suvtk 格式
  python submission_pipeline.py --work-dir ./ --mode suvtk

  # 仅生成 Sequin 格式 (Cenote-Taker3 风格)
  python submission_pipeline.py --work-dir ./ --mode sequin
"""
    )

    # 上游输入
    parser.add_argument('--work-dir', required=True,
                        help='上游输出目录 (含 all_plant_viruses.fasta)')
    parser.add_argument('--fasta',
                        help='输入 FASTA (默认 work-dir/all_plant_viruses.fasta)')

    # 数据库
    parser.add_argument('--suvtk-db', help='suvtk 数据库路径')
    parser.add_argument('--rvdb-db', help='RVDB 数据库目录 (U-RVDBv31.0-prot.*)')

    # 运行参数
    parser.add_argument('--run-title', default='viral_submission', help='运行标题')
    parser.add_argument('--mode', choices=['suvtk', 'sequin', 'both'], default='both',
                        help='输出模式 (默认 both)')
    parser.add_argument('-t', '--threads', type=int, default=40)

    # 公共元数据 (public_metadata_pipeline)
    parser.add_argument('--fetch-metadata', action='store_true',
                        help='自动从 NCBI/GSA 获取 SRA 元数据 (调用 gsa_sra.info.py)')
    parser.add_argument('--metadata-dir',
                        help='public_metadata_pipeline 输出目录 (含 Global_Unified_Metadata_Core13.tsv)')

    # 配置
    parser.add_argument('--config', help='提交配置文件 (YAML), 用于生成 authorset.sbt')

    # 交互
    parser.add_argument('--interactive', action='store_true',
                        help='交互式补全缺失元数据 (VAPiD 风格)')
    parser.add_argument('--pipeline-type', choices=['auto', 'discovery', 'analysis'], default='auto',
                        help='上游管道类型: auto=自动检测, discovery=新病毒, analysis=已知病毒')
    parser.add_argument('--ref-info', help='已知病毒参考信息 (auto_known_virus 的 ref_info.tsv)')

    # 重跑
    parser.add_argument('--force-taxonomy', action='store_true', help='强制重跑 taxonomy')
    parser.add_argument('--force-features', action='store_true', help='强制重跑 features')
    parser.add_argument('--force-hypo', action='store_true', help='强制重跑 hypothetical')
    parser.add_argument('--export-metadata', action='store_true',
                        help='额外导出统一元数据 CSV (SeqSender 兼容)')

    args = parser.parse_args()

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
        log.info("  无元数据文件 (source.src 留占位符待填)")
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

        suvtk_script = SCRIPT_DIR / "suvtk_submission.py"

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
                f"-o {work_dir}/genbank_submission/"
            )
            run(cmd, log, "sequin_builder")

    # ── 统一元数据 CSV 导出 (SeqSender 兼容) ──
    if args.export_metadata and tax_tsv:
        log.info("─" * 40 + "\n  导出统一元数据 CSV (SeqSender 兼容)\n" + "─" * 40)

        meta_script = SCRIPT_DIR / "unified_metadata.py"
        cmd = (
            f"python {meta_script} "
            f"--taxonomy {tax_tsv} "
            f"--run-title {args.run_title} "
            f"-o {work_dir}/submission_metadata/"
        )
        if core13.exists():
            cmd += f" --metadata {core13}"
        if args.config:
            cmd += f" --config {args.config}"
        if args.interactive:
            cmd += " --interactive"

        run(cmd, log, "unified_metadata")

    # ═══ 最终报告 ═══
    log.info("=" * 60)
    log.info("完成!")
    log.info("  输出目录: %s", work_dir)
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
