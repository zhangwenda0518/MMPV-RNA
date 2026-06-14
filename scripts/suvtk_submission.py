#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
suvtk_submission.py — MMPV-RNA → GenBank 提交准备管道 v1.0
==========================================================

基于 suvtk 工具链, 将 MMPV-RNA 管道产出的新病毒和已知病毒序列
准备为可直接上传 GenBank 的 .sqn 提交文件。

完整流程:
  1. suvtk taxonomy    → ICTV 分类学分配 + 基因组类型预测
  2. suvtk features    → ORF预测 (pyrodigal) + 功能注释 (BFVD) + .tbl特征表
  3. 准备 source.src   → 样本元数据 (需用户按模板填写)
  4. suvtk comments    → 整合 MIUVIG 结构化注释
  5. suvtk table2asn   → 生成最终 .sqn 提交文件

辅助:
  - hypothetical protein 分析
  - co-occurrence 分段病毒检测 (已知病毒)
  - CheckV 质量报告整合

依赖:
  pip install suvtk
  suvtk download-database  (已完成: ~/database/virus-db/suvtk_db/)

用法:
  # 新病毒 (plant novel viruses from rescue)
  python suvtk_submission.py novel \
      --fasta $OUT/08_Rescue/Plant/centroids/final_centroids.fasta \
      --taxonomy $OUT/05_Taxonomy/integrated/final_integrated_classification.tsv \
      --host $OUT/06_HostPrediction/ensemble_host_summary.tsv \
      --checkv $OUT/08_Rescue/checkv/ \
      --suvtk-db ~/database/virus-db/suvtk_db/ \
      --output ./genbank_submission/novel/ \
      -t 40

  # 已知病毒 (known viruses from auto_known_virus)
  python suvtk_submission.py known \
      --fasta $OUT/known_viruses/3_Virus_assemblies_final/ \
      --summary $OUT/known_viruses/1_FastViromeExplorer/summary/best.summary.tsv \
      --ref-info /db/ref_info.tsv \
      --suvtk-db ~/database/virus-db/suvtk_db/ \
      --output ./genbank_submission/known/ \
      -t 40
"""

import argparse
import os
import sys
import subprocess
import shutil
import logging
from pathlib import Path
from datetime import datetime

# ══════════════════════════════════════════════════════════════
# 日志
# ══════════════════════════════════════════════════════════════

def setup_logger(out_dir, level="INFO"):
    logger = logging.getLogger("SuvtkSubmit")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    os.makedirs(out_dir, exist_ok=True)

    ch = logging.StreamHandler()
    ch.setLevel(getattr(logging, level.upper(), logging.INFO))
    ch.setFormatter(logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(ch)

    fh = logging.FileHandler(os.path.join(out_dir, 'suvtk_submission.log'), encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s %(message)s'))
    logger.addHandler(fh)

    return logger


def run(cmd, log, step_name, check=True):
    """执行命令, 记录日志"""
    log.info("[%s] 执行...", step_name)
    log.debug("  CMD: %s", cmd)
    try:
        subprocess.run(cmd, shell=True, check=check, executable='/bin/bash')
        log.info("[%s] ✓ 完成", step_name)
        return True
    except subprocess.CalledProcessError as e:
        log.error("[%s] ✗ 失败 (exit=%d)", step_name, e.returncode)
        if not check:
            return False
        sys.exit(1)


# ══════════════════════════════════════════════════════════════
# 步骤 1: suvtk taxonomy — ICTV 分类分配
# ══════════════════════════════════════════════════════════════

def run_taxonomy(fasta_path, out_dir, suvtk_db, threads, log):
    """suvtk taxonomy: MMseqs2 LCA → ICTV taxonomy + genome type prediction"""
    tax_dir = os.path.join(out_dir, "1_taxonomy")
    os.makedirs(tax_dir, exist_ok=True)

    tax_tsv = os.path.join(tax_dir, "taxonomy.tsv")
    if os.path.exists(tax_tsv) and os.path.getsize(tax_tsv) > 100:
        log.info("[1/5] taxonomy — 已有结果, 跳过")
    else:
        cmd = (
            f"suvtk taxonomy "
            f"-i {fasta_path} "
            f"-o {tax_dir} "
            f"-d {suvtk_db} "
            f"-s 0.7 "
            f"-t {threads}"
        )
        run(cmd, log, "suvtk taxonomy")

    # 验证输出
    for f in ["taxonomy.tsv", "miuvig_taxonomy.tsv"]:
        fp = os.path.join(tax_dir, f)
        if os.path.exists(fp):
            n = sum(1 for _ in open(fp)) - 1
            log.info("  %s: %d 条记录", f, n)

    return tax_dir


# ══════════════════════════════════════════════════════════════
# 步骤 2: suvtk features — ORF预测 + 功能注释 + .tbl
# ══════════════════════════════════════════════════════════════

def run_features(fasta_path, tax_dir, out_dir, suvtk_db, threads, log):
    """suvtk features: pyrodigal ORF → MMseqs2 BFVD 注释 → .tbl"""
    feat_dir = os.path.join(out_dir, "2_features")
    os.makedirs(feat_dir, exist_ok=True)

    tax_tsv = os.path.join(tax_dir, "taxonomy.tsv")
    tbl_files = list(Path(feat_dir).glob("*.tbl"))

    if tbl_files and os.path.getsize(tax_tsv) > 100:
        log.info("[2/5] features — 已有结果, 跳过")
    else:
        cmd = (
            f"suvtk features "
            f"-i {fasta_path} "
            f"-o {feat_dir} "
            f"-d {suvtk_db} "
            f"--coding-complete "
            f"--taxonomy {tax_tsv} "
            f"-t {threads}"
        )
        run(cmd, log, "suvtk features")

    # 验证输出
    for pat in ["*.tbl", "*.fna", "*.faa"]:
        files = list(Path(feat_dir).glob(pat))
        if files:
            log.info("  %s: %d 个文件", pat, len(files))

    return feat_dir


# ══════════════════════════════════════════════════════════════
# 步骤 2.5: 假定蛋白分析 (可选)
# ══════════════════════════════════════════════════════════════

def run_hypothetical_analysis(feat_dir, out_dir, log):
    """分析假定蛋白, 可选 BLAST 进一步搜索"""
    hypo_dir = os.path.join(out_dir, "2.5_hypothetical")
    os.makedirs(hypo_dir, exist_ok=True)

    tbl_files = list(Path(feat_dir).glob("*.tbl"))
    faa_files = list(Path(feat_dir).glob("*.faa"))

    if not tbl_files or not faa_files:
        log.warning("[2.5] hypothetical — 缺少 .tbl 或 .faa, 跳过")
        return None

    for tbl in tbl_files:
        # 找对应的 faa
        base = tbl.stem
        faa = Path(feat_dir) / f"{base}.faa"
        if not faa.exists():
            faa = faa_files[0]  # fallback

        hypo_out = os.path.join(hypo_dir, f"{base}_hypothetical.faa")
        if os.path.exists(hypo_out) and os.path.getsize(hypo_out) > 0:
            log.info("[2.5] hypothetical — %s 已有结果, 跳过", base)
            continue

        cmd = (
            f"python ~/bin/analyze_hypothetical.py "
            f"-t {tbl} "
            f"-f {faa} "
            f"-o {hypo_out}"
        )
        run(cmd, log, f"hypothetical: {base}", check=False)

    return hypo_dir


# ══════════════════════════════════════════════════════════════
# 步骤 3: 生成 source.src 模板 (用户需手动补充样本信息)
# ══════════════════════════════════════════════════════════════

TEMPLATE_SOURCE_SRC = """# source.src — GenBank 提交源信息模板
# 由 suvtk_submission.py 自动生成 | {timestamp}
#
# ⚠️ 请根据实际实验记录修改以下占位符字段:
#   - isolate:        唯一分离株标识符 (同病毒的不同片段必须相同)
#   - collection_date: 样本采集日期 (格式: DD-Mmm-YYYY, 如 15-Jun-2024)
#   - geo_loc_name:    采集地点 (格式: Country:Region, 如 China:Jiangsu)
#   - lat_lon:         经纬度 (格式: 32.06 N 118.79 E)
#   - bioproject:      BioProject 登录号 (如 PRJNA123456)
#   - biosample:       BioSample 登录号 (如 SAMN12345678)
#   - sra:             SRA 登录号 (如有, 如 SRR12345678)
#   - metagenome_source: 宏基因组来源 (如 "soil metagenome")
#   - segment:         分段病毒的片段编号 (非分段留空)
#
# 字段说明: https://landerdc.github.io/suvtk/index.html

Sequence_ID\tOrganism\tIsolate\tCollection_date\tgeo_loc_name\tLat_Lon\tBioproject\tBiosample\tSRA\tMetagenomic\tMetagenome_source\tSegment
{source_lines}
"""


def generate_source_template(tax_tsv, out_dir, host_tsv=None, log=None):
    """从 taxonomy.tsv 生成 source.src 模板"""
    src_file = os.path.join(out_dir, "3_metadata", "source.src.template")
    meta_dir = os.path.join(out_dir, "3_metadata")
    os.makedirs(meta_dir, exist_ok=True)

    # 读取 taxonomy
    lines = []
    with open(tax_tsv) as f:
        header = f.readline().strip().split("\t")
        tax_idx = header.index("taxonomy") if "taxonomy" in header else 1
        contig_idx = header.index("contig") if "contig" in header else 0

        for line in f:
            cols = line.strip().split("\t")
            contig = cols[contig_idx]
            tax = cols[tax_idx] if len(cols) > tax_idx else "Viruses"
            # 默认值 (用户需修改)
            isolate = contig.split("_")[0] if "_" in contig else contig[:20]
            lines.append(
                f"{contig}\t{tax}\t{isolate}_isolate\tDD-Mmm-YYYY\t"
                f"Country:Region\tXX.XX_N_XXX.XX_E\t"
                f"PRJNAXXXXXX\tSAMNXXXXXXXX\tSRRXXXXXXXX\t"
                f"TRUE\tsoil_metagenome\t"
            )

    source_content = TEMPLATE_SOURCE_SRC.format(
        timestamp=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        source_lines="\n".join(lines)
    )

    with open(src_file, "w", encoding="utf-8") as f:
        f.write(source_content)

    if log:
        log.info("[3/5] source.src 模板 → %s (%d 条序列)", src_file, len(lines))
    return src_file


# ══════════════════════════════════════════════════════════════
# 步骤 3.5: 生成 MIUVIG 元数据文件
# ══════════════════════════════════════════════════════════════

MIUVIG_TEMPLATE = """# miuvig.tsv — MIUVIG 标准全局元数据
# 由 suvtk_submission.py 自动生成 | {timestamp}
# 参考: https://standardsingenomics.org/miuvig/
#
# ⚠️ 请根据实际实验修改以下字段

sample_id\t{miuvig_fields}
"""

ASSEMBLY_TEMPLATE = """# assembly.tsv — GenBank 组装注释信息
# 由 suvtk_submission.py 自动生成 | {timestamp}

Sequencing_Technology\tAssembly_Method\tAssembly_Name\tAssembly_Software\tCoverage
{assembly_lines}
"""


def generate_miuvig_metadata(out_dir, log, seq_type="rna-short", assembler="MEGAHIT"):
    """生成 miuvig.tsv 和 assembly.tsv 模板"""
    meta_dir = os.path.join(out_dir, "3_metadata")
    os.makedirs(meta_dir, exist_ok=True)

    # miuvig.tsv
    miuvig_file = os.path.join(meta_dir, "miuvig.tsv")
    if not os.path.exists(miuvig_file):
        with open(miuvig_file, "w") as f:
            f.write(f"""# miuvig.tsv — MIUVIG 标准全局元数据
# 参考: https://standardsingenomics.org/miuvig/
#
# ⚠️ 请根据实际实验修改以下占位字段
sample_id\tviral_enrichment\tsequencing_platform\tsequencing_method\tassembly_software\tassembly_method\tquality_check_software
ALL\trRNA_depletion\tIllumina_NovaSeq\t{seq_type.upper()}\t{assembler}\tmetaSPAdes_MEGAHIT\tCheckV
""")
        log.info("[3.5] miuvig.tsv → %s", miuvig_file)

    # assembly.tsv
    asm_file = os.path.join(meta_dir, "assembly.tsv")
    if not os.path.exists(asm_file):
        with open(asm_file, "w") as f:
            f.write(f"""# assembly.tsv — GenBank 组装注释信息
Sequencing_Technology\tAssembly_Method\tAssembly_Name\tAssembly_Software\tCoverage
Illumina_NovaSeq\tmetatranscriptomic_assembly\tMMPV-RNA_v2.3\t{assembler}\tNOT_PROVIDED
""")
        log.info("[3.5] assembly.tsv → %s", asm_file)

    return meta_dir


# ══════════════════════════════════════════════════════════════
# 步骤 4: suvtk comments — 整合 MIUVIG 注释
# ══════════════════════════════════════════════════════════════

def run_comments(tax_dir, feat_dir, meta_dir, out_dir, log, checkv_dir=None):
    """suvtk comments: 合并 taxonomy + features + MIUVIG 元数据 → .cmt"""
    cmt_dir = os.path.join(out_dir, "4_comments")
    os.makedirs(cmt_dir, exist_ok=True)

    tax_tsv = os.path.join(tax_dir, "miuvig_taxonomy.tsv")
    feat_tsv = os.path.join(feat_dir, "miuvig_features.tsv")
    miuvig_tsv = os.path.join(meta_dir, "miuvig.tsv")
    asm_tsv = os.path.join(meta_dir, "assembly.tsv")
    cmt_file = os.path.join(cmt_dir, "output.cmt")

    if os.path.exists(cmt_file) and os.path.getsize(cmt_file) > 100:
        log.info("[4/5] comments — 已有结果, 跳过")
        return cmt_dir

    # 检查必需输入
    missing = []
    for f, name in [(tax_tsv, "miuvig_taxonomy.tsv"),
                     (feat_tsv, "miuvig_features.tsv"),
                     (miuvig_tsv, "miuvig.tsv"),
                     (asm_tsv, "assembly.tsv")]:
        if not os.path.exists(f):
            missing.append(name)

    if missing:
        log.error("[4/5] comments — 缺少输入文件: %s", ", ".join(missing))
        log.error("  请先运行步骤 1/2, 并检查 3_metadata/ 中的模板文件")
        return None

    cmd = (
        f"suvtk comments "
        f"--taxonomy {tax_tsv} "
        f"--features {feat_tsv} "
        f"--miuvig {miuvig_tsv} "
        f"--assembly {asm_tsv} "
        f"-o {cmt_dir}"
    )
    if checkv_dir:
        qs = os.path.join(checkv_dir, "completeness.tsv")
        if os.path.exists(qs):
            cmd += f" --quality {qs}"

    run(cmd, log, "suvtk comments")
    return cmt_dir


# ══════════════════════════════════════════════════════════════
# 步骤 5: suvtk table2asn — 生成 .sqn
# ══════════════════════════════════════════════════════════════

def run_table2asn(feat_dir, meta_dir, cmt_dir, out_dir, log):
    """suvtk table2asn: 打包 → .sqn + 验证"""
    sqn_dir = os.path.join(out_dir, "5_submission")
    os.makedirs(sqn_dir, exist_ok=True)

    # 查找输入文件
    fna_files = list(Path(feat_dir).glob("*.fna"))
    tbl_files = list(Path(feat_dir).glob("*.tbl"))
    src_file = os.path.join(meta_dir, "source.src")  # 用户已填写的版本
    src_template = os.path.join(meta_dir, "source.src.template")

    src = src_file if os.path.exists(src_file) else src_template
    cmt_file = os.path.join(cmt_dir, "output.cmt") if cmt_dir else None

    if not fna_files:
        log.error("[5/5] table2asn — 缺少 .fna 文件 (来自步骤2)")
        return None

    if not tbl_files:
        log.error("[5/5] table2asn — 缺少 .tbl 文件 (来自步骤2)")
        return None

    sqn_out = os.path.join(sqn_dir, "submission.sqn")
    if os.path.exists(sqn_out) and os.path.getsize(sqn_out) > 1000:
        log.info("[5/5] table2asn — 已有结果, 跳过")
        return sqn_dir

    log.info("[5/5] table2asn — 生成 .sqn 提交文件")
    log.info("  ⚠️ 请确保已修改 source.src 中的占位符信息!")
    log.info("  ⚠️ 请确保已从 NCBI 下载 template.sbt 文件!")
    log.info("  参考: https://submit.ncbi.nlm.nih.gov/genbank/template/submission/")

    # 构建命令 (每个 .fna + 对应的 .tbl)
    for fna in fna_files:
        base = fna.stem
        tbl = Path(feat_dir) / f"{base}.tbl"
        if not tbl.exists():
            tbl = tbl_files[0]  # fallback to first .tbl

        cmd = (
            f"suvtk table2asn "
            f"--fasta {fna} "
            f"--features {tbl} "
            f"--source {src} "
            f"-o {sqn_dir}"
        )
        if cmt_file and os.path.exists(cmt_file):
            cmd += f" --comments {cmt_file}"

        run(cmd, log, f"table2asn: {base}", check=False)

    return sqn_dir


# ══════════════════════════════════════════════════════════════
# 主流程
# ══════════════════════════════════════════════════════════════

def prepare_novel_viruses(args, log):
    """新病毒提交准备 (来自 rescue pipeline)"""
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)

    log.info("=" * 60)
    log.info("新病毒 GenBank 提交准备 (Plant Novel Viruses)")
    log.info("  FASTA:  %s", args.fasta)
    log.info("  Output: %s", out)
    log.info("=" * 60)

    # Step 1: taxonomy
    tax_dir = run_taxonomy(args.fasta, str(out), args.suvtk_db, args.threads, log)

    # Step 2: features
    feat_dir = run_features(args.fasta, tax_dir, str(out), args.suvtk_db, args.threads, log)

    # Step 2.5: hypothetical protein analysis (optional)
    run_hypothetical_analysis(feat_dir, str(out), log)

    # Step 3: metadata templates
    tax_tsv = os.path.join(tax_dir, "taxonomy.tsv")
    generate_source_template(tax_tsv, str(out), log=log)
    meta_dir = generate_miuvig_metadata(str(out), log)

    # Step 4: comments
    cmt_dir = run_comments(tax_dir, feat_dir, meta_dir, str(out), log,
                           checkv_dir=args.checkv if hasattr(args, 'checkv') else None)

    # Step 5: table2asn — 需要用户先填写 source.src
    src_filled = os.path.join(meta_dir, "source.src")
    if os.path.exists(src_filled):
        run_table2asn(feat_dir, meta_dir, cmt_dir, str(out), log)
    else:
        log.info("")
        log.info("=" * 60)
        log.info("⚠️  请先完成以下步骤再运行 table2asn:")
        log.info("  1. 编辑 %s/source.src.template → source.src", meta_dir)
        log.info("  2. 填写: isolate, collection_date, geo_loc_name, lat_lon,")
        log.info("           bioproject, biosample, sra, metagenome_source")
        log.info("  3. 从 NCBI 下载 template.sbt:")
        log.info("     https://submit.ncbi.nlm.nih.gov/genbank/template/submission/")
        log.info("  4. 重新运行: python suvtk_submission.py novel --step table2asn ...")
        log.info("=" * 60)

    log.info("")
    log.info("完成! 输出目录: %s", out)
    _print_output_tree(out)


def prepare_known_viruses(args, log):
    """已知病毒提交准备 (来自 auto_known_virus pipeline)"""
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)

    log.info("=" * 60)
    log.info("已知病毒 GenBank 提交准备 (Known Viruses)")
    log.info("  FASTA:  %s", args.fasta)
    log.info("  Output: %s", out)
    log.info("=" * 60)

    # 收集所有全长组装 FASTA
    fasta_dir = Path(args.fasta)
    all_fastas = list(fasta_dir.rglob("final.fasta")) + list(fasta_dir.rglob("*.fasta"))

    if not all_fastas:
        log.error("未找到全长组装结果! 请检查 --fasta 路径")
        sys.exit(1)

    # 合并所有 FASTA
    combined_fa = out / "combined_known_viruses.fasta"
    with open(combined_fa, "w") as cf:
        for fa in all_fastas:
            sample_tag = fa.parent.name if fa.parent.name != fasta_dir.name else ""
            with open(fa) as inf:
                for line in inf:
                    if line.startswith(">") and sample_tag:
                        cf.write(f">{sample_tag}|{line[1:]}")
                    else:
                        cf.write(line)
    n = sum(1 for l in open(combined_fa) if l.startswith(">"))
    log.info("  合并 %d 个 FASTA → %d 条序列 → %s", len(all_fastas), n, combined_fa)

    # Step 1: taxonomy
    tax_dir = run_taxonomy(str(combined_fa), str(out), args.suvtk_db, args.threads, log)

    # 对已知病毒, 可选 co-occurrence 分析
    if hasattr(args, 'summary') and args.summary:
        log.info("[1.5] co-occurrence — 分段病毒关联分析")
        summary_file = Path(args.summary)
        if summary_file.exists():
            cooc_dir = out / "1.5_cooccurrence"
            os.makedirs(cooc_dir, exist_ok=True)
            cmd = (
                f"suvtk co-occurrence "
                f"--abundance {summary_file} "
                f"-o {cooc_dir}"
            )
            run(cmd, log, "suvtk co-occurrence", check=False)

    # Step 2: features
    feat_dir = run_features(str(combined_fa), tax_dir, str(out), args.suvtk_db, args.threads, log)

    # Step 3: metadata
    tax_tsv = os.path.join(tax_dir, "taxonomy.tsv")
    generate_source_template(tax_tsv, str(out), log=log)
    meta_dir = generate_miuvig_metadata(str(out), log)

    # Step 4: comments
    cmt_dir = run_comments(tax_dir, feat_dir, meta_dir, str(out), log)

    # Step 5: table2asn
    src_filled = os.path.join(meta_dir, "source.src")
    if os.path.exists(src_filled):
        run_table2asn(feat_dir, meta_dir, cmt_dir, str(out), log)
    else:
        log.info("")
        log.info("⚠️  请编辑 source.src.template 后重新运行")

    log.info("完成! 输出目录: %s", out)
    _print_output_tree(out)


def _print_output_tree(out_dir):
    """打印输出目录树"""
    out = Path(out_dir)
    print(f"\n{'='*60}")
    print(f"输出目录: {out}")
    for d in sorted(out.rglob("*")):
        if d.is_file():
            size = d.stat().st_size
            if size > 1024 * 1024:
                s = f"{size/1024/1024:.1f}MB"
            elif size > 1024:
                s = f"{size/1024:.1f}KB"
            else:
                s = f"{size}B"
            rel = d.relative_to(out)
            print(f"  {rel} ({s})")
    print(f"{'='*60}")


# ══════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════

def main():
    p = argparse.ArgumentParser(
        description="MMPV-RNA → GenBank 提交准备管道 (基于 suvtk)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    sub = p.add_subparsers(dest="mode", help="提交模式")

    # ── novel 模式 ──
    p_novel = sub.add_parser("novel", help="新病毒提交准备 (来自 rescue pipeline)")
    p_novel.add_argument("--fasta", required=True,
                         help="centroids FASTA (e.g. 08_Rescue/Plant/centroids/final_centroids.fasta)")
    p_novel.add_argument("--taxonomy", help="final_integrated_classification.tsv (可选, 用于交叉验证)")
    p_novel.add_argument("--host", help="ensemble_host_summary.tsv (可选, 用于宿主信息)")
    p_novel.add_argument("--checkv", help="CheckV 结果目录 (可选, 用于质量报告)")
    p_novel.add_argument("--suvtk-db", required=True, help="suvtk 数据库路径")
    p_novel.add_argument("--output", "-o", default="./genbank_submission/novel", help="输出目录")
    p_novel.add_argument("--threads", "-t", type=int, default=40, help="线程数")

    # ── known 模式 ──
    p_known = sub.add_parser("known", help="已知病毒提交准备 (来自 auto_known_virus pipeline)")
    p_known.add_argument("--fasta", required=True,
                         help="全长组装 FASTA 目录 (e.g. known_viruses/3_Virus_assemblies_final/)")
    p_known.add_argument("--summary", help="best.summary.tsv (可选, 用于 co-occurrence)")
    p_known.add_argument("--ref-info", help="ref_info.tsv (可选)")
    p_known.add_argument("--suvtk-db", required=True, help="suvtk 数据库路径")
    p_known.add_argument("--output", "-o", default="./genbank_submission/known", help="输出目录")
    p_known.add_argument("--threads", "-t", type=int, default=40, help="线程数")

    # ── 分步模式 (在已有输出目录上继续) ──
    p_step = sub.add_parser("step", help="从已有输出目录继续某一步")
    p_step.add_argument("--step", required=True,
                        choices=["taxonomy", "features", "hypothetical", "metadata", "comments", "table2asn"],
                        help="执行步骤")
    p_step.add_argument("--work-dir", required=True, help="已有输出目录")
    p_step.add_argument("--suvtk-db", required=True, help="suvtk 数据库路径")
    p_step.add_argument("--fasta", help="输入 FASTA (taxonomy/features 步骤需要)")
    p_step.add_argument("--threads", "-t", type=int, default=40)

    args = p.parse_args()

    if not args.mode:
        p.print_help()
        sys.exit(0)

    out = Path(args.output if hasattr(args, 'output') else args.work_dir).resolve()
    log = setup_logger(str(out))

    if args.mode == "novel":
        prepare_novel_viruses(args, log)
    elif args.mode == "known":
        prepare_known_viruses(args, log)
    elif args.mode == "step":
        # 分步执行 (在已有目录上继续)
        work = Path(args.work_dir)
        if args.step == "taxonomy":
            run_taxonomy(args.fasta, str(work), args.suvtk_db, args.threads, log)
        elif args.step == "features":
            tax_dir = work / "1_taxonomy"
            run_features(args.fasta, str(tax_dir), str(work), args.suvtk_db, args.threads, log)
        elif args.step == "metadata":
            tax_tsv = work / "1_taxonomy" / "taxonomy.tsv"
            generate_source_template(str(tax_tsv), str(work), log=log)
            generate_miuvig_metadata(str(work), log)
        elif args.step == "comments":
            run_comments(str(work / "1_taxonomy"), str(work / "2_features"),
                        str(work / "3_metadata"), str(work), log)
        elif args.step == "table2asn":
            cmt_dir = work / "4_comments" if (work / "4_comments").exists() else None
            run_table2asn(str(work / "2_features"), str(work / "3_metadata"),
                         str(cmt_dir) if cmt_dir else None, str(work), log)


if __name__ == "__main__":
    main()
