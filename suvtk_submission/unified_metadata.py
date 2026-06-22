#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unified_metadata.py — 统一元数据模板生成器 v1.0
================================================

借鉴 SeqSender 的 "一张表覆盖 GenBank + BioSample" 设计:
  一次填写 ~25 个核心字段, 自动同时生成:
    - source.src   (GenBank source modifiers)
    - miuvig.tsv    (MIUVIG 全局参数)
    - assembly.tsv  (Assembly 注释)
    - BioSample CSV (可选, 用于 NCBI BioSample 批量注册)
    - authorset.sbt (可选, 作者信息)

字段前缀约定 (同 SeqSender):
  src-*  → source.src 列 (如 src-Isolate, src-geo_loc_name)
  cmt-*  → .cmt 注释 (如 cmt-Assembly_Method)
  bs-*   → BioSample 字段 (如 bs-isolate, bs-geo_loc_name)
  无前缀 → 通用字段 (organism, collection_date, bioproject)

用法:
  # 从 suvtk 产出生成统一模板
  python unified_metadata.py \\
      --taxonomy 1_taxonomy/taxonomy.tsv \\
      --metadata Global_Unified_Metadata_Core13.tsv \\
      --run-title my_project \\
      -o ./submission/

输出:
  submission/
  ├── unified_metadata.csv      ← 一张表 (可 Excel 编辑)
  ├── source.src                ← 自动从 CSV 生成
  ├── source_individual/        ← 每个病毒独立 source.src
  ├── biosample_template.tsv    ← NCBI BioSample 批量模板
  ├── miuvig.tsv
  └── assembly.tsv
"""

import argparse
import os
import sys
import csv
import re
import logging
from pathlib import Path
from datetime import datetime
from collections import OrderedDict

import pandas as pd
from tqdm import tqdm


# ══════════════════════════════════════════════════════════════
# 统一元数据表定义 (精简自 SeqSender, 保留核心 25 字段)
# ══════════════════════════════════════════════════════════════

UNIFIED_COLUMNS = OrderedDict([
    # === 通用 (Required) ===
    ("organism",               {"required": True,  "group": "general", "desc": "NCBI Taxonomy 物种名"}),
    ("sequence_name",          {"required": True,  "group": "general", "desc": "FASTA 中的序列 ID"}),
    ("authors",                {"required": True,  "group": "general", "desc": "引用作者: Last, First; ..."}),
    ("collection_date",        {"required": True,  "group": "general", "desc": "采集日期 (YYYY-MM-DD 或 YYYY-MM 或 YYYY)"}),
    ("bioproject",             {"required": True,  "group": "general", "desc": "NCBI BioProject ID (PRJNA...)"}),

    # === source.src 字段 (src-*) ===
    ("src-Isolate",            {"required": True,  "group": "src", "desc": "唯一分离株标识 (分段病毒各片段必须相同)"}),
    ("src-geo_loc_name",       {"required": True,  "group": "src", "desc": "采集地点 (格式 Country:Region, 如 China:Ningxia)"}),
    ("src-Lat_Lon",            {"required": True,  "group": "src", "desc": "经纬度 (如 38.47 N 106.27 E)"}),
    ("src-Host",               {"required": False, "group": "src", "desc": "宿主物种名"}),
    ("src-Segment",            {"required": False, "group": "src", "desc": "分段病毒的片段编号 (非分段留空)"}),
    ("src-Isolation-source",   {"required": False, "group": "src", "desc": "分离来源描述"}),
    ("src-Note",               {"required": False, "group": "src", "desc": "额外备注"}),
    ("src-Tissue_type",        {"required": False, "group": "src", "desc": "组织类型"}),
    ("src-Collected_by",       {"required": False, "group": "src", "desc": "采集人"}),

    # === GenBank 提交字段 ===
    ("gb-sample_name",         {"required": True,  "group": "gb", "desc": "GenBank 记录名 (≤50字符)"}),
    ("gb-title",               {"required": False, "group": "gb", "desc": "提交标题 (NCBI 门户显示)"}),
    ("sra",                    {"required": False, "group": "gb", "desc": "SRA 登录号 (SRR...)"}),
    ("biosample",              {"required": True,  "group": "gb", "desc": "BioSample ID (SAMN...)"}),

    # === 结构化注释字段 (cmt-*) ===
    ("cmt-Assembly_Method",    {"required": True,  "group": "cmt", "desc": "组装方法 (如 MEGAHIT v1.2.9)"}),
    ("cmt-Sequencing_Technology", {"required": True, "group": "cmt", "desc": "测序平台 (如 Illumina NovaSeq 6000)"}),
    ("cmt-Genome_Coverage",    {"required": False, "group": "cmt", "desc": "基因组覆盖度 (如 42.5x)"}),
    ("cmt-Annotation_Pipeline", {"required": False, "group": "cmt", "desc": "注释流程 (如 MMPV-RNA v2.3 + suvtk v0.1.1)"}),

    # === BioSample 字段 (bs-*, 如需自动注册) ===
    ("bs-isolate",             {"required": False, "group": "bs", "desc": "BioSample: 分离株标识"}),
    ("bs-geo_loc_name",        {"required": False, "group": "bs", "desc": "BioSample: 采集地点"}),
    ("bs-host",                {"required": False, "group": "bs", "desc": "BioSample: 宿主"}),
    ("bs-isolation_source",    {"required": False, "group": "bs", "desc": "BioSample: 分离来源"}),
])

REQUIRED_COLS = [k for k, v in UNIFIED_COLUMNS.items() if v["required"]]


# ══════════════════════════════════════════════════════════════

def setup_logger(out_dir):
    logger = logging.getLogger("UnifiedMeta")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    ch = logging.StreamHandler()
    ch.setFormatter(logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(ch)
    os.makedirs(out_dir, exist_ok=True)
    return logger


def sanitize_filename(name, max_len=50):
    name = re.sub(r'[^\w\s-]', '', name)
    name = re.sub(r'\s+', '_', name)
    return name[:max_len]


def load_taxonomy(tax_tsv):
    """加载 taxonomy.tsv"""
    df = pd.read_csv(tax_tsv, sep='\t')
    col_map = {}
    for c in df.columns:
        if c.lower() in ('contig', 'seq_id', 'sequence_id'):
            col_map['contig'] = c
        elif c.lower() in ('taxonomy', 'tax'):
            col_map['taxonomy'] = c
    if 'contig' not in col_map:
        col_map['contig'] = df.columns[0]
    if 'taxonomy' not in col_map:
        col_map['taxonomy'] = df.columns[1] if len(df.columns) > 1 else df.columns[0]

    seqs = []
    for _, row in df.iterrows():
        seqs.append({
            'contig': str(row[col_map['contig']]),
            'taxonomy': str(row[col_map['taxonomy']]) if pd.notna(row[col_map['taxonomy']]) else 'unclassified viruses',
        })
    return seqs


def load_metadata(meta_file):
    """加载 Core13 元数据表"""
    if not meta_file or not os.path.exists(meta_file):
        return {}
    df = pd.read_csv(meta_file, sep=None, engine='python')
    lookup = {}
    for _, row in df.iterrows():
        run = str(row.get('Run', '')).strip()
        if not run or run.lower() in ('nan', 'not_provided', ''):
            continue
        lookup[run] = {
            'collection_date': str(row.get('CollectionDate', '')) if pd.notna(row.get('CollectionDate')) else '',
            'geo_loc_name': str(row.get('Location', '')) if pd.notna(row.get('Location')) else '',
            'bioproject': str(row.get('BioProject', '')) if pd.notna(row.get('BioProject')) else '',
            'biosample': str(row.get('BioSample', '')) if pd.notna(row.get('BioSample')) else '',
            'host': str(row.get('ScientificName', '')) if pd.notna(row.get('ScientificName')) else '',
            'tissue': str(row.get('Tissue', '')) if pd.notna(row.get('Tissue')) else '',
        }
    return lookup


def extract_sra(contig):
    m = re.match(r'([SC]RR\d+)', contig)
    return m.group(1) if m else 'UNKNOWN'


def generate_metadata_csv(seqs, meta_lookup, args, log):
    """生成统一元数据 CSV"""

    rows = []
    for s in seqs:
        contig = s['contig']
        taxonomy = s['taxonomy']
        sra = extract_sra(contig)
        meta = meta_lookup.get(sra, {})

        # 自动填充已知值
        collection_date = meta.get('collection_date', '')
        if collection_date and collection_date.lower() in ('nan', 'not_provided', ''):
            collection_date = ''

        geo_loc = meta.get('geo_loc_name', '')
        if geo_loc and geo_loc.lower() in ('nan', 'not_provided', ''):
            geo_loc = ''

        bioproject = meta.get('bioproject', args.bioproject or '')
        if bioproject and bioproject.lower() in ('nan', 'prjnaXXXXXX', ''):
            bioproject = args.bioproject or ''

        biosample = meta.get('biosample', '')
        if biosample and biosample.lower() in ('nan', 'samnXXXXXXXX', ''):
            biosample = ''

        isolate = f"{sanitize_filename(taxonomy)}_{sra}"

        rows.append({
            "organism": taxonomy,
            "sequence_name": contig,
            "authors": args.authors or "Last, First",
            "collection_date": collection_date or "YYYY-MM-DD",
            "bioproject": bioproject or "PRJNAXXXXXX",

            "src-Isolate": isolate,
            "src-geo_loc_name": geo_loc or "Country:Region",
            "src-Lat_Lon": args.lat_lon or "XX.XX N XXX.XX E",
            "src-Host": meta.get('host', args.host or ''),
            "src-Segment": "",
            "src-Isolation-source": args.isolation_source or "",
            "src-Note": "",
            "src-Tissue_type": meta.get('tissue', args.tissue or ''),
            "src-Collected_by": args.collected_by or "",

            "gb-sample_name": isolate[:50],
            "gb-title": args.title or f"{taxonomy} genome sequencing",
            "sra": sra,
            "biosample": biosample or "SAMNXXXXXXXX",

            "cmt-Assembly_Method": args.assembler or "MEGAHIT v1.2.9",
            "cmt-Sequencing_Technology": args.sequencer or "Illumina NovaSeq 6000",
            "cmt-Genome_Coverage": args.coverage or "",
            "cmt-Annotation_Pipeline": args.pipeline or "MMPV-RNA v2.3 + suvtk v0.1.1",

            "bs-isolate": isolate,
            "bs-geo_loc_name": geo_loc or "Country:Region",
            "bs-host": meta.get('host', ''),
            "bs-isolation_source": args.metagenome_source or "plant virome",
        })

    df = pd.DataFrame(rows, columns=list(UNIFIED_COLUMNS.keys()))
    return df


def export_source_src(df, out_dir, log):
    """从统一 CSV 生成 source.src"""
    src_dir = Path(out_dir) / "source_individual"
    src_dir.mkdir(parents=True, exist_ok=True)

    # 完整 source.src
    src_all = Path(out_dir) / "source.src"
    with open(src_all, 'w', encoding='utf-8') as f:
        f.write("Sequence_ID\tOrganism\tIsolate\tCollection_date\tgeo_loc_name\tLat_Lon\t"
                "Bioproject\tBiosample\tSRA\tMetagenomic\tMetagenome_source\tSegment\n")
        for _, row in df.iterrows():
            f.write(f"{row['sequence_name']}\t{row['organism']}\t{row['src-Isolate']}\t"
                    f"{row['collection_date']}\t{row['src-geo_loc_name']}\t{row['src-Lat_Lon']}\t"
                    f"{row['bioproject']}\t{row['biosample']}\t{row['sra']}\tTRUE\t"
                    f"{row['src-Isolation-source'] or 'plant virome'}\t{row['src-Segment']}\n")
    log.info("  → %s (%d 条)", src_all, len(df))

    # 按 virus 拆分
    by_org = df.groupby('organism')
    for org, group in by_org:
        fname = f"source_{sanitize_filename(org)}.src"
        src_virus = src_dir / fname
        with open(src_virus, 'w', encoding='utf-8') as f:
            f.write("Sequence_ID\tOrganism\tIsolate\tCollection_date\tgeo_loc_name\tLat_Lon\t"
                    "Bioproject\tBiosample\tSRA\tMetagenomic\tMetagenome_source\tSegment\n")
            for _, row in group.iterrows():
                f.write(f"{row['sequence_name']}\t{row['organism']}\t{row['src-Isolate']}\t"
                        f"{row['collection_date']}\t{row['src-geo_loc_name']}\t{row['src-Lat_Lon']}\t"
                        f"{row['bioproject']}\t{row['biosample']}\t{row['sra']}\tTRUE\t"
                        f"{row['src-Isolation-source'] or 'plant virome'}\t{row['src-Segment']}\n")
        log.info("  → %s (%d 条)", src_virus, len(group))


def export_biosample_csv(df, out_dir, log):
    """生成 NCBI BioSample 批量提交模板 (Pathogen.cl.1.0 包)"""
    # BioSample Pathogen.cl.1.0 包的必填字段
    bs_mapping = {
        'sample_name':      'gb-sample_name',   # BioSample 记录名
        'organism':         'organism',          # 物种名
        'collection_date':  'collection_date',   # 采集日期
        'geo_loc_name':     'bs-geo_loc_name',   # 采集地点
        'isolation_source': 'bs-isolation_source', # 分离来源
        'isolate':          'bs-isolate',        # 分离株
        'host':             'bs-host',           # 宿主
        'bioproject':       'bioproject',        # BioProject
    }
    bs_df = pd.DataFrame()
    for bs_field, src_col in bs_mapping.items():
        if src_col in df.columns:
            if src_col.startswith('bs-'):
                bs_df[bs_field] = df[src_col]
            else:
                bs_df[bs_field] = df[src_col]
    bs_path = Path(out_dir) / "biosample_template.tsv"
    bs_df.to_csv(bs_path, sep='\t', index=False)
    log.info("  → %s (%d 条, Pathogen.cl.1.0 格式)", bs_path, len(bs_df))


def create_authorset_sbt(config_dict, out_dir, metadata_df, log):
    """从配置生成 authorset.sbt (SeqSender 风格)"""
    sbt_path = Path(out_dir) / "authorset.sbt"

    submitter = config_dict.get('Submitter', {})
    org = config_dict.get('Organization', {})
    addr = config_dict.get('Address', {})
    authors_str = metadata_df['authors'].iloc[0] if 'authors' in metadata_df.columns else 'Author, First'

    # 解析作者列表
    author_list = [a.strip() for a in authors_str.replace(';', ',').split(',') if a.strip()]

    with open(sbt_path, 'w', encoding='utf-8') as f:
        f.write("Submit-block ::= {\n")
        f.write("  contact {\n")
        f.write("    contact {\n")
        f.write("      name name {\n")
        f.write(f"        last \"{submitter.get('Last', 'LastName')}\",\n")
        f.write(f"        first \"{submitter.get('First', 'FirstName')}\"\n")
        f.write("      },\n")
        f.write("      affil std {\n")
        f.write(f"        affil \"{addr.get('Affil', 'Institution')}\",\n")
        f.write(f"        div \"{addr.get('Div', 'Department')}\",\n")
        f.write(f"        city \"{addr.get('City', 'City')}\",\n")
        f.write(f"        sub \"{addr.get('Sub', 'State')}\",\n")
        f.write(f"        country \"{addr.get('Country', 'Country')}\",\n")
        f.write(f"        street \"{addr.get('Street', '')}\",\n")
        f.write(f"        email \"{submitter.get('Email', 'email@example.com')}\",\n")
        f.write(f"        postal-code \"{addr.get('Postal_Code', '00000')}\"\n")
        f.write("      }\n")
        f.write("    }\n")
        f.write("  },\n")
        f.write("  cit {\n")
        f.write("    authors {\n")
        f.write("      names std {\n")
        for i, author in enumerate(author_list, 1):
            parts = author.strip().split()
            last = parts[-1] if parts else "Author"
            first = parts[0] if len(parts) >= 2 else ""
            f.write("        {\n")
            f.write("          name name {\n")
            f.write(f"            last \"{last}\",\n")
            f.write(f"            first \"{first}\"\n")
            f.write("          }\n")
            if i == len(author_list):
                f.write("        }\n")
            else:
                f.write("        },\n")
        f.write("      },\n")
        f.write("      affil std {\n")
        f.write(f"        affil \"{addr.get('Affil', 'Institution')}\",\n")
        f.write(f"        div \"{addr.get('Div', 'Department')}\",\n")
        f.write(f"        city \"{addr.get('City', 'City')}\",\n")
        f.write(f"        sub \"{addr.get('Sub', 'State')}\",\n")
        f.write(f"        country \"{addr.get('Country', 'Country')}\",\n")
        f.write(f"        street \"{addr.get('Street', '')}\"\n")
        f.write("      }\n")
        f.write("    },\n")
        f.write(f"    title \"{config_dict.get('Publication_Title', 'Viral genome sequencing and assembly')}\",\n")
        f.write(f"    status \"{config_dict.get('Publication_Status', 'Unpublished')}\"\n")
        f.write("  }\n")
        f.write("}\n")
    log.info("  → %s", sbt_path)


def create_submission_log(out_dir, run_title, log):
    """创建提交追踪日志 (submission_log.csv)"""
    log_path = Path(out_dir) / "submission_log.csv"
    exists = log_path.exists()
    with open(log_path, 'a' if exists else 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        if not exists:
            writer.writerow(['timestamp', 'submission_name', 'status', 'n_sequences', 'biosamples', 'sqn_file', 'notes'])
        writer.writerow([
            datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            run_title,
            'FILES_GENERATED',
            '', '', '', 'Metadata files ready. Edit unified_metadata.csv, then re-run with --validate to check.'
        ])
    log.info("  → %s (更新)", log_path)


def load_config(config_path):
    """加载提交配置文件 (YAML 或 JSON)"""
    config = {
        'Submitter': {'First': 'FirstName', 'Last': 'LastName', 'Email': 'email@example.com'},
        'Address': {'Affil': 'Institution', 'Div': 'Department', 'City': 'City',
                     'Sub': 'State', 'Country': 'Country', 'Street': '', 'Postal_Code': '00000'},
        'Publication_Title': 'Viral genome sequencing and assembly',
        'Publication_Status': 'Unpublished',
        'Spuid_Namespace': '',
        'GenBank_Auto_Remove_Failed_Samples': True,
        'Specified_Release_Date': '',
    }
    if config_path and os.path.exists(config_path):
        if config_path.endswith('.yaml') or config_path.endswith('.yml'):
            import yaml
            with open(config_path) as f:
                user_config = yaml.safe_load(f)
        elif config_path.endswith('.json'):
            import json
            with open(config_path) as f:
                user_config = json.load(f)
        else:
            return config
        # 扁平化嵌套结构
        if 'Submission' in user_config and 'NCBI' in user_config['Submission']:
            ncbi = user_config['Submission']['NCBI']
        else:
            ncbi = user_config
        if 'Description' in ncbi:
            desc = ncbi['Description']
            if 'Organization' in desc:
                org = desc['Organization']
                config['Submitter'] = org.get('Submitter', config['Submitter'])
                config['Address'] = org.get('Address', config['Address'])
                config['Address']['Affil'] = org.get('Name', config['Address']['Affil'])
            config['Publication_Title'] = ncbi.get('Publication_Title', config['Publication_Title'])
            config['Publication_Status'] = ncbi.get('Publication_Status', config['Publication_Status'])
            config['Specified_Release_Date'] = ncbi.get('Specified_Release_Date', '')
            config['Spuid_Namespace'] = ncbi.get('Spuid_Namespace', '')
    return config


def inject_interactive_metadata(out_dir, df, log):
    """交互式补全缺失元数据 (VAPiD 风格)"""
    import shutil
    # 备份原始 CSV
    csv_path = Path(out_dir) / "unified_metadata.csv"
    bak_path = Path(out_dir) / "unified_metadata.csv.bak"
    if not bak_path.exists():
        shutil.copy(csv_path, bak_path)
        log.info("  备份: %s", bak_path)

    modified = False
    for col in REQUIRED_COLS:
        if col not in df.columns:
            continue
        mask = df[col].isna() | (df[col].astype(str).str.strip() == '') | \
               df[col].astype(str).str.contains('XXXX|YYYY|PRJNAXXXX|Country:Region|SAMNXXXXXXXX|Author')
        if not mask.any():
            continue

        n_missing = mask.sum()
        unique_vals = df.loc[mask, col].unique() if hasattr(df.loc[mask, col], 'unique') else []

        print(f"\n{'─'*50}")
        print(f"  [{col}] — {UNIFIED_COLUMNS.get(col, {}).get('desc', '')}")
        print(f"  缺失 {n_missing}/{len(df)} 条, 当前值: {unique_vals[:3]}")
        print(f"\n  输入新值 (回车保留原值, '{col}=ALL' 批填所有):")
        val = input(f"  > ").strip()

        if val:
            if f"{col}=ALL" in val.replace(' ', ''):
                # 批量填充
                all_val = val.split('=', 1)[1].strip()
                df.loc[mask, col] = all_val
                log.info("  批量填充 %s = %s (%d 条)", col, all_val, n_missing)
            else:
                # 只更新有具体序列的 — 这里交互式逐条填太慢, 简化为一键全填
                df.loc[mask, col] = val
                log.info("  填充 %s = %s (%d 条)", col, val, n_missing)
            modified = True

    if modified:
        df.to_csv(csv_path, index=False, encoding='utf-8-sig')
        log.info("  → 已保存更新: %s", csv_path)
    else:
        log.info("  无需修改")
    return df


def export_miuvig_assembly(out_dir, log, assembler, sequencer, enrichment):
    """生成 miuvig.tsv 和 assembly.tsv"""
    miuvig_path = Path(out_dir) / "miuvig.tsv"
    with open(miuvig_path, 'w') as f:
        f.write("MIUVIG_parameter\tvalue\n")
        f.write(f"source_uvig\tviral fraction metagenome (virome)\n")
        f.write(f"assembly_software\t{assembler}\n")
        f.write(f"sequencing_platform\t{sequencer}\n")
        if enrichment:
            f.write(f"virus_enrich_appr\t{enrichment}\n")
    log.info("  → %s", miuvig_path)

    asm_path = Path(out_dir) / "assembly.tsv"
    with open(asm_path, 'w') as f:
        f.write("Assembly_parameter\tvalue\n")
        f.write("StructuredCommentPrefix\tAssembly-Data\n")
        f.write(f"Assembly Method\t{assembler}\n")
        f.write(f"Sequencing Technology\t{sequencer}\n")
    log.info("  → %s", asm_path)


def validate_csv(df, out_dir, log):
    """检查必填字段, 输出验证报告"""
    report_path = Path(out_dir) / "validation_report.txt"
    issues = []
    for col in REQUIRED_COLS:
        if col not in df.columns:
            issues.append(f"  [MISSING] 缺少必填列: {col}")
            continue
        missing = df[col].isna() | (df[col].astype(str).str.strip() == '') | \
                  df[col].astype(str).str.contains('XXXX|YYYY-MM-DD|Country:Region|PRJNAXXXX|SAMNXXXXXXXX')
        if missing.any():
            seqs = df.loc[missing, 'sequence_name'].tolist()
            issues.append(f"  [EMPTY] {col}: {len(seqs)} 条序列未填 ({seqs[:3]}...)")

    with open(report_path, 'w') as f:
        f.write("元数据验证报告\n")
        f.write(f"生成: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"总序列: {len(df)}\n")
        f.write(f"必填列: {len(REQUIRED_COLS)}\n")
        f.write("─" * 50 + "\n")
        if issues:
            f.write(f"⚠️  发现 {len(issues)} 个问题:\n")
            for i in issues:
                f.write(i + "\n")
            f.write("\n提交前请修复以上问题!\n")
        else:
            f.write("✓ 所有必填字段已填写, 可以提交\n")

    if issues:
        log.warning("  验证: %d 个问题 → %s", len(issues), report_path)
    else:
        log.info("  验证: 全部通过 ✓ → %s", report_path)


def main():
    parser = argparse.ArgumentParser(
        description="unified_metadata.py — 统一元数据模板生成器 v1.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--taxonomy', required=True, help='suvtk taxonomy.tsv')
    parser.add_argument('--metadata', help='Global_Unified_Metadata_Core13.tsv (自动填充)')
    parser.add_argument('--run-title', default='viral_submission', help='运行标题')
    parser.add_argument('-o', '--output', default='./submission/', help='输出目录')
    parser.add_argument('--assembler', default='MEGAHIT v1.2.9')
    parser.add_argument('--sequencer', default='Illumina NovaSeq 6000')
    parser.add_argument('--enrichment', default='rRNA depletion')
    parser.add_argument('--pipeline', default='MMPV-RNA v2.3 + suvtk v0.1.1')
    parser.add_argument('--authors', help='引用作者')
    parser.add_argument('--bioproject', help='BioProject ID')
    parser.add_argument('--host', help='宿主物种')
    parser.add_argument('--lat-lon', help='经纬度')
    parser.add_argument('--isolation-source', help='分离来源')
    parser.add_argument('--collected-by', help='采集人')
    parser.add_argument('--tissue', help='组织类型')
    parser.add_argument('--coverage', help='覆盖度')
    parser.add_argument('--title', help='提交标题')
    parser.add_argument('--metagenome-source', help='宏基因组来源')
    parser.add_argument('--config', help='提交配置文件 (YAML), 如 SeqSender 的 seqsender_config.yaml')
    parser.add_argument('--interactive', action='store_true',
                        help='交互式补全缺失元数据 (VAPiD 风格)')
    parser.add_argument('--skip-validate', action='store_true', help='跳过验证')
    args = parser.parse_args()

    out_dir = Path(args.output)
    log = setup_logger(str(out_dir))

    log.info("=" * 60)
    log.info("统一元数据模板生成: %s", args.run_title)

    # 加载数据
    seqs = load_taxonomy(args.taxonomy)
    meta_lookup = load_metadata(args.metadata) if args.metadata else {}
    log.info("  序列: %d, 元数据映射: %d", len(seqs), len(meta_lookup))

    # 生成统一 CSV
    df = generate_metadata_csv(seqs, meta_lookup, args, log)
    csv_path = out_dir / "unified_metadata.csv"
    df.to_csv(csv_path, index=False, encoding='utf-8-sig')
    log.info("  → %s (%d 行 × %d 列)", csv_path, len(df), len(df.columns))

    # 导出 source.src
    export_source_src(df, out_dir, log)

    # 导出 BioSample 模板
    export_biosample_csv(df, out_dir, log)

    # 导出 miuvig/assembly
    export_miuvig_assembly(out_dir, log, args.assembler, args.sequencer, args.enrichment)

    # 加载配置 → 生成 authorset.sbt
    config = load_config(args.config)
    create_authorset_sbt(config, out_dir, df, log)

    # 提交追踪日志
    create_submission_log(out_dir, args.run_title, log)

    # 交互式补全 (VAPiD 风格)
    if args.interactive:
        df = inject_interactive_metadata(out_dir, df, log)

    # 验证
    if not args.skip_validate:
        validate_csv(df, out_dir, log)

    # 打印结果
    auto_filled = sum(1 for c in REQUIRED_COLS if c in df.columns and not df[c].astype(str).str.contains('XXXX|YYYY|PRJNA|Country:Region').any())
    print(f"\n{'='*60}")
    print(f"  统一元数据文件已生成")
    print(f"  {'='*60}")
    print(f"  输出: {out_dir}/")
    print(f"    unified_metadata.csv     ← Excel 编辑 (一张表驱动一切)")
    print(f"    source.src               ← GenBank 全部序列")
    print(f"    source_individual/       ← 每个病毒独立")
    print(f"    biosample_template.tsv   ← NCBI BioSample 批量注册")
    print(f"    authorset.sbt            ← 作者信息 (从 config 生成)")
    print(f"    submission_log.csv       ← 提交追踪日志")
    print(f"    miuvig.tsv / assembly.tsv")
    print(f"    validation_report.txt    ← 字段校验")
    print(f"  {'='*60}")
    print(f"  序列数: {len(df)}")
    print(f"  必填字段自动填充: {auto_filled}/{len(REQUIRED_COLS)}")
    print(f"  {'='*60}")
    print(f"")
    print(f"  提交三步走:")
    print(f"  1. 编辑 unified_metadata.csv 占位符 (或用 --interactive 交互补全)")
    print(f"  2. 在 NCBI 注册 BioProject → 填入 bioproject 列")
    print(f"  3. 用 biosample_template.tsv 批量注册 BioSample → 填入 biosample 列")
    print(f"  4. 运行 suvtk table2asn → .sqn → 邮件发 gb-sub@ncbi.nlm.nih.gov")


if __name__ == '__main__':
    main()
