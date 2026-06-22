#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sequin_builder.py — GenBank Sequin 提交文件构建器 v1.0
======================================================

借鉴 Cenote-Taker3 的每病毒独立文件 + FASTA 头部嵌入元数据模式,
从 suvtk 产出构建可直接提交 NCBI 的 Sequin 文件包。

数据流:
  suvtk taxonomy + features + analyze_hypothetical
    │
    ├── 按病毒分组 → 每个病毒独立 .fsa/.tbl/.cmt/.gbf
    ├── .fsa 头部嵌入元数据 (无需单独 source.src)
    ├── .cmt 自动生成 (MIUVIG + Assembly)
    └── 可选 table2asn → .sqn

输出结构:
  {run_title}/
  ├── {run_title}_virus_summary.tsv              # 每个病毒的汇总表
  ├── {run_title}_run_arguments.txt              # 参数记录
  ├── {run_title}_virus_sequences.fna            # 所有核酸序列
  ├── {run_title}_virus_AA.faa                   # 所有蛋白序列
  │
  └── sequin_files/                              # 提交文件
      ├── {organism}_{isolate}.fsa               # Sequin FASTA
      ├── {organism}_{isolate}.tbl               # 5列特征表
      ├── {organism}_{isolate}.cmt               # 结构化注释
      ├── {organism}_{isolate}.src               # source modifier (备用)
      └── {organism}_{isolate}.sqn               # (可选) 最终提交文件

用法:
  python sequin_builder.py \\
      --taxonomy 1_taxonomy/taxonomy.tsv \\
      --miuvig-tax 1_taxonomy/miuvig_taxonomy.tsv \\
      --features 2_features/featuretable_updated.tbl \\
      --proteins 2_features/proteins_updated.faa \\
      --fasta 2_features/reoriented_nucleotide_sequences.fna \\
      --miuvig-feat 2_features/miuvig_features.tsv \\
      --metadata metadata/Global_Unified_Metadata_Core13.tsv \\
      --run-title ningxiagouqi_plant_virome \\
      --sequencer "Illumina NovaSeq 6000" \\
      --assembler "MEGAHIT v1.2.9" \\
      --enrichment "rRNA depletion" \\
      --metagenome-source "plant virome" \\
      -o ./genbank_submission/
"""

import argparse
import os
import sys
import json
import re
import random
import string
import logging
from pathlib import Path
from datetime import datetime
from collections import OrderedDict

import pandas as pd
from tqdm import tqdm

# ══════════════════════════════════════════════════════════════
# 日志
# ══════════════════════════════════════════════════════════════

def setup_logger(out_dir):
    logger = logging.getLogger("SequinBuilder")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    os.makedirs(out_dir, exist_ok=True)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(ch)

    fh = logging.FileHandler(os.path.join(out_dir, 'sequin_builder.log'), encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s %(message)s'))
    logger.addHandler(fh)

    return logger


# ══════════════════════════════════════════════════════════════
# 工具函数
# ══════════════════════════════════════════════════════════════

def sanitize_filename(name, max_len=50):
    """病毒名 → 合法文件名"""
    name = re.sub(r'[^\w\s-]', '', name)
    name = re.sub(r'\s+', '_', name)
    return name[:max_len]


def random_id(length=5):
    """Cenote-Taker3 风格的随机 5 位 ID"""
    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=length))


def load_taxonomy(tax_tsv):
    """加载 taxonomy.tsv, 返回 [(contig, taxonomy), ...]"""
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
    """加载 Core13 元数据表, 返回 {SRA: {collection_date, geo_loc_name, ...}}"""
    if not meta_file or not os.path.exists(meta_file):
        return {}

    df = pd.read_csv(meta_file, sep=None, engine='python')
    lookup = {}
    for _, row in df.iterrows():
        run = str(row.get('Run', '')).strip()
        if not run or run.lower() in ('nan', 'not_provided', ''):
            continue
        lookup[run] = {
            'collection_date': str(row.get('CollectionDate', '')),
            'geo_loc_name': str(row.get('Location', '')),
            'bioproject': str(row.get('BioProject', '')),
            'biosample': str(row.get('BioSample', '')),
            'tissue': str(row.get('Tissue', '')),
            'source': str(row.get('Source', '')),
        }
    return lookup


def load_featuretable(tbl_path):
    """解析 featuretable.tbl, 返回 {seq_id: [cds_lines]}"""
    records = OrderedDict()
    current_seq = None
    current_lines = []

    with open(tbl_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>Feature '):
                if current_seq and current_lines:
                    records[current_seq] = current_lines
                current_seq = line.split('>Feature ', 1)[1].strip()
                current_lines = []
            else:
                current_lines.append(line.rstrip('\n'))
        if current_seq and current_lines:
            records[current_seq] = current_lines

    return records


def extract_sra(contig):
    """从 contig 名提取 SRA/CRR 前缀"""
    m = re.match(r'([SC]RR\d+)', contig)
    return m.group(1) if m else 'UNKNOWN'


def build_genome_type(taxonomy):
    """根据 taxonomy 推断 genome type (简化版)"""
    tax = taxonomy.lower()
    if any(k in tax for k in ['rhabdoviridae', 'phenuiviridae', 'tospoviridae',
                                'fimoviridae', 'ophioviridae', 'aspiviridae',
                                'negarnaviricota']):
        return ('ssRNA(-)', 'undetermined')
    elif any(k in tax for k in ['partitiviridae', 'chrysoviridae', 'totiviridae',
                                  'endoraviridae', 'amalgaviridae', 'reoviridae']):
        return ('dsRNA', 'undetermined')
    elif any(k in tax for k in ['geminiviridae', 'nanoviridae']):
        return ('ssDNA', 'undetermined')
    else:
        return ('ssRNA(+)', 'undetermined')


# ══════════════════════════════════════════════════════════════
# 构建器
# ══════════════════════════════════════════════════════════════

class SequinBuilder:
    def __init__(self, args, log):
        self.args = args
        self.log = log

        # 路径
        self.tax_tsv = Path(args.taxonomy)
        self.miuvig_tax_tsv = Path(args.miuvig_tax)
        self.tbl_path = Path(args.features)
        self.faa_path = Path(args.proteins)
        self.fna_path = Path(args.fasta)
        self.miuvig_feat_tsv = Path(args.miuvig_feat)

        self.run_title = args.run_title
        self.out_dir = Path(args.output)
        self.sequin_dir = self.out_dir / self.run_title / 'sequin_files'
        self.sequin_dir.mkdir(parents=True, exist_ok=True)

        # 元数据参数
        self.sequencer = args.sequencer
        self.assembler = args.assembler
        self.enrichment = args.enrichment
        self.mol_type = args.mol_type
        self.metagenome_source = args.metagenome_source
        self.annotation_pipeline = args.annotation_pipeline

        # 加载数据
        self.seqs = load_taxonomy(self.tax_tsv)
        self.metadata = load_metadata(args.metadata) if args.metadata else {}
        self.tbl_records = load_featuretable(self.tbl_path)

        # 加载核酸序列
        self.fna_seqs = {}
        if self.fna_path.exists():
            from Bio import SeqIO
            for rec in SeqIO.parse(str(self.fna_path), 'fasta'):
                self.fna_seqs[rec.id] = str(rec.seq)

        # 按病毒分组
        self.virus_groups = self._group_by_virus()

    def _group_by_virus(self):
        """两级分组: 病毒 taxonomy → 病毒内序列"""
        viruses = OrderedDict()
        for s in self.seqs:
            tax = s['taxonomy'].strip()
            if tax not in viruses:
                viruses[tax] = {'taxonomy': tax, 'seqs': [], 'sras': set()}
            viruses[tax]['seqs'].append(s)
            viruses[tax]['sras'].add(extract_sra(s['contig']))
        return viruses

    def build_all(self):
        """构建所有病毒的文件"""
        summary_rows = []
        total_seqs = 0

        self.log.info("=" * 60)
        self.log.info("Sequin 文件构建: %d 种病毒, %d 条序列",
                      len(self.virus_groups), len(self.seqs))

        for tax_key, vgroup in tqdm(sorted(self.virus_groups.items()),
                                      desc="Building sequin files"):
            total_seqs += len(vgroup['seqs'])
            row = self._build_one_virus(tax_key, vgroup)
            if row:
                summary_rows.append(row)

        # 保存汇总表
        self._save_summary(summary_rows)

        # 保存合并 FASTA
        self._save_combined_fasta()

        # 保存 run_arguments.txt
        self._save_run_args(total_seqs)

        self.log.info("=" * 60)
        self.log.info("完成! 输出目录: %s", self.out_dir / self.run_title)
        self.log.info("  病毒种类: %d, 序列总数: %d", len(self.virus_groups), total_seqs)
        self.log.info("  提交文件: %s/", self.sequin_dir)

    def _build_one_virus(self, tax_key, vgroup):
        """为一个病毒构建全套 Sequin 文件"""
        taxonomy = vgroup['taxonomy']
        fname_base = sanitize_filename(taxonomy)
        isolate_id = f"{fname_base}_{random_id()}"
        ssras = sorted(vgroup['sras'])
        genome_type, genome_struc = build_genome_type(taxonomy)

        # 收集该病毒所有序列的元数据
        meta = self._collect_meta(ssras)

        # ── 1. .fsa (Sequin FASTA with embedded metadata) ──
        fsa_path = self.sequin_dir / f"{fname_base}_{isolate_id}.fsa"
        self._write_fsa(fsa_path, vgroup, taxonomy, isolate_id, meta, genome_type)

        # ── 2. .tbl (feature table for this virus) ──
        tbl_path = self.sequin_dir / f"{fname_base}_{isolate_id}.tbl"
        self._write_tbl(tbl_path, vgroup)

        # ── 3. .cmt (structured comment) ──
        cmt_path = self.sequin_dir / f"{fname_base}_{isolate_id}.cmt"
        self._write_cmt(cmt_path, vgroup, meta)

        # ── 4. .src (source modifier, 备用) ──
        src_path = self.sequin_dir / f"{fname_base}_{isolate_id}.src"
        self._write_src(src_path, vgroup, taxonomy, isolate_id, meta)

        # 汇总行
        return {
            'virus_taxonomy': taxonomy,
            'isolate': isolate_id,
            'n_sequences': len(vgroup['seqs']),
            'n_sras': len(vgroup['sras']),
            'sras': ','.join(ssras),
            'genome_type': genome_type,
            'genome_struc': genome_struc,
            'fsa': fsa_path.name,
            'tbl': tbl_path.name,
            'cmt': cmt_path.name,
            'src': src_path.name,
            'collection_date': meta.get('collection_date', ''),
            'geo_loc_name': meta.get('geo_loc_name', ''),
        }

    def _collect_meta(self, sras):
        """从 metadata 收集多个 SRA 的元数据, 取最早的 collection_date"""
        meta = {
            'collection_date': '',
            'geo_loc_name': '',
            'bioproject': '',
            'biosample': '',
        }
        dates = []
        for sra in sras:
            if sra in self.metadata:
                m = self.metadata[sra]
                if m.get('collection_date') and m['collection_date'].lower() not in ('nan', 'not_provided', ''):
                    dates.append(m['collection_date'])
                if not meta['geo_loc_name'] and m.get('geo_loc_name') and m['geo_loc_name'].lower() not in ('nan', 'not_provided', ''):
                    meta['geo_loc_name'] = m['geo_loc_name']
                if not meta['bioproject'] and m.get('bioproject') and m['bioproject'].lower() not in ('nan', 'prjnaXXXXXX'.lower(), ''):
                    meta['bioproject'] = m['bioproject']
                if not meta['biosample'] and m.get('biosample') and m['biosample'].lower() not in ('nan', 'samnXXXXXXXX'.lower(), ''):
                    meta['biosample'] = m['biosample']

        if dates:
            meta['collection_date'] = min(dates)
        return meta

    def _write_fsa(self, path, vgroup, taxonomy, isolate_id, meta, genome_type):
        """写 Sequin FASTA (.fsa) — 参考 Cenote-Taker3 格式"""
        with open(path, 'w', encoding='utf-8') as f:
            for s in vgroup['seqs']:
                contig = s['contig']
                sra = extract_sra(contig)
                seq = self.fna_seqs.get(contig, '')

                if not seq:
                    self.log.warning("  未找到核酸序列: %s", contig)
                    continue

                # 构建 Cenote-Taker3 风格的头部
                # [organism=...] [gcode=1] [topology=linear]
                # [isolation_source=...] [collection_date=...]
                # [metagenome_source=...] [SRA=...] [Biosample=...] [Bioproject=...]
                # [moltype=cRNA]
                gcode = 11 if 'ssRNA' in genome_type else 1
                mol_type = 'cRNA' if 'RNA' in genome_type else 'DNA'
                topology = 'linear'

                header_parts = [contig]
                header_parts.append(f"[organism={taxonomy}]")
                header_parts.append(f"[isolate={isolate_id}]")
                header_parts.append(f"[gcode={gcode}]")
                header_parts.append(f"[topology={topology}]")
                if meta.get('geo_loc_name'):
                    header_parts.append(f"[geo_loc_name={meta['geo_loc_name']}]")
                if meta.get('collection_date'):
                    header_parts.append(f"[collection_date={meta['collection_date']}]")
                if self.metagenome_source:
                    header_parts.append(f"[metagenome_source={self.metagenome_source}]")
                if sra != 'UNKNOWN':
                    header_parts.append(f"[SRA={sra}]")
                if meta.get('biosample'):
                    header_parts.append(f"[Biosample={meta['biosample']}]")
                if meta.get('bioproject'):
                    header_parts.append(f"[Bioproject={meta['bioproject']}]")
                header_parts.append(f"[moltype={mol_type}]")

                f.write(f">{' '.join(header_parts)}\n")
                f.write(f"{seq}\n")

    def _write_tbl(self, path, vgroup):
        """写特征表 (.tbl) — 只包含该病毒的 CDS"""
        with open(path, 'w', encoding='utf-8') as f:
            for s in vgroup['seqs']:
                contig = s['contig']
                if contig in self.tbl_records:
                    f.write(f">Feature {contig}\n")
                    for line in self.tbl_records[contig]:
                        f.write(line + '\n')

    def _write_cmt(self, path, vgroup, meta):
        """写结构化注释 (.cmt)"""
        with open(path, 'w', encoding='utf-8') as f:
            f.write("StructuredCommentPrefix\t##Genome-Assembly-Data-START##\n")
            f.write(f"Assembly Method\t{self.assembler}\n")
            f.write(f"Sequencing Technology\t{self.sequencer}\n")
            f.write(f"Annotation Pipeline\t{self.annotation_pipeline}\n")
            if self.enrichment:
                f.write(f"Viral Enrichment\t{self.enrichment}\n")
            f.write("\n")

            # MIUVIG 部分
            f.write("StructuredCommentPrefix\t##MIUVIG-Data-START##\n")
            f.write(f"pred_genome_type\t{build_genome_type(vgroup['taxonomy'])[0]}\n")
            f.write(f"metagenomic\tTRUE\n")
            if self.metagenome_source:
                f.write(f"metagenome_source\t{self.metagenome_source}\n")
            if meta.get('collection_date'):
                f.write(f"collection_date\t{meta['collection_date']}\n")
            n_contigs = len(vgroup['seqs'])
            f.write(f"number_contig\t{n_contigs}\n")

    def _write_src(self, path, vgroup, taxonomy, isolate_id, meta):
        """写 source modifier (.src) — 备用格式"""
        with open(path, 'w', encoding='utf-8') as f:
            f.write("Sequence_ID\tOrganism\tIsolate\tCollection_date\tgeo_loc_name\t"
                     "Bioproject\tBiosample\tSRA\tMetagenomic\tMetagenome_source\tSegment\n")
            for s in vgroup['seqs']:
                sra = extract_sra(s['contig'])
                date = meta.get('collection_date', '')
                geo = meta.get('geo_loc_name', '')
                bp = meta.get('bioproject', '')
                bs = meta.get('biosample', '')
                f.write(f"{s['contig']}\t{taxonomy}\t{isolate_id}\t{date}\t{geo}\t\t"
                         f"{bp}\t{bs}\t{sra}\tTRUE\t{self.metagenome_source}\t\n")

    def _save_summary(self, rows):
        """保存病毒汇总表"""
        if not rows:
            return
        df = pd.DataFrame(rows)
        path = self.out_dir / self.run_title / f"{self.run_title}_virus_summary.tsv"
        df.to_csv(path, sep='\t', index=False)
        self.log.info("  → %s (%d 种病毒)", path, len(rows))

    def _save_combined_fasta(self):
        """保存合并核酸和蛋白 FASTA"""
        # 核酸
        fna_out = self.out_dir / self.run_title / f"{self.run_title}_virus_sequences.fna"
        if self.fna_path.exists():
            with open(fna_out, 'w') as fout, open(self.fna_path) as fin:
                fout.write(fin.read())
            self.log.info("  → %s", fna_out)

        # 蛋白
        faa_out = self.out_dir / self.run_title / f"{self.run_title}_virus_AA.faa"
        if self.faa_path.exists():
            with open(faa_out, 'w') as fout, open(self.faa_path) as fin:
                fout.write(fin.read())
            self.log.info("  → %s", faa_out)

    def _save_run_args(self, total_seqs):
        """保存运行参数"""
        path = self.out_dir / self.run_title / f"{self.run_title}_run_arguments.txt"
        with open(path, 'w') as f:
            f.write(f"Sequin 提交文件构建参数\n")
            f.write(f"{'='*50}\n")
            f.write(f"运行时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"运行标题: {self.run_title}\n")
            f.write(f"病毒种类: {len(self.virus_groups)}\n")
            f.write(f"序列总数: {total_seqs}\n")
            f.write(f"\n输入文件:\n")
            f.write(f"  taxonomy:     {self.tax_tsv}\n")
            f.write(f"  features:     {self.tbl_path}\n")
            f.write(f"  proteins:     {self.faa_path}\n")
            f.write(f"  nucleotides:  {self.fna_path}\n")
            f.write(f"\n参数:\n")
            f.write(f"  sequencer:    {self.sequencer}\n")
            f.write(f"  assembler:    {self.assembler}\n")
            f.write(f"  enrichment:   {self.enrichment}\n")
            f.write(f"  mol_type:     {self.mol_type}\n")
            f.write(f"  metagenome_source: {self.metagenome_source}\n")
            f.write(f"  pipeline:     {self.annotation_pipeline}\n")
        self.log.info("  → %s", path)


# ══════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="sequin_builder.py — GenBank Sequin 提交文件构建器 v1.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python sequin_builder.py \\
      --taxonomy 1_taxonomy/taxonomy.tsv \\
      --miuvig-tax 1_taxonomy/miuvig_taxonomy.tsv \\
      --features 2_features/featuretable_updated.tbl \\
      --proteins 2_features/proteins_updated.faa \\
      --fasta 2_features/reoriented_nucleotide_sequences.fna \\
      --miuvig-feat 2_features/miuvig_features.tsv \\
      --metadata metadata/Global_Unified_Metadata_Core13.tsv \\
      --run-title my_project \\
      -o ./genbank_submission/
"""
    )

    # 必需输入
    parser.add_argument('--taxonomy', required=True, help='suvtk taxonomy.tsv')
    parser.add_argument('--miuvig-tax', required=True, help='miuvig_taxonomy.tsv')
    parser.add_argument('--features', required=True, help='featuretable.tbl (updated)')
    parser.add_argument('--proteins', required=True, help='proteins.faa (updated)')
    parser.add_argument('--fasta', required=True, help='reoriented_nucleotide_sequences.fna')
    parser.add_argument('--miuvig-feat', required=True, help='miuvig_features.tsv')

    # 可选输入
    parser.add_argument('--metadata', help='Global_Unified_Metadata_Core13.tsv (自动填充元数据)')

    # 运行参数
    parser.add_argument('--run-title', default='viral_submission', help='运行标题 (目录名)')
    parser.add_argument('-o', '--output', default='./genbank_submission/', help='输出根目录')

    # 元数据参数
    parser.add_argument('--sequencer', default='Illumina NovaSeq 6000', help='测序平台')
    parser.add_argument('--assembler', default='MEGAHIT v1.2.9', help='组装软件')
    parser.add_argument('--enrichment', default='rRNA depletion', help='病毒富集方法')
    parser.add_argument('--mol-type', default='cRNA', help='分子类型 (cRNA/DNA)')
    parser.add_argument('--metagenome-source', default='plant virome', help='宏基因组来源')
    parser.add_argument('--annotation-pipeline', default='MMPV-RNA v2.3 + suvtk v0.1.1', help='注释管线')

    args = parser.parse_args()

    out_dir = Path(args.output)
    log = setup_logger(str(out_dir / args.run_title))

    builder = SequinBuilder(args, log)
    builder.build_all()


if __name__ == '__main__':
    main()
