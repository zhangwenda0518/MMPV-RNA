#!/usr/bin/env python3
"""
tbl2gb.py — 将 featuretable_updated.tbl + FASTA 转为 GenBank (.gb) 文件
=====================================================================
输入:
  analyze_hypothetical/featuretable_updated.tbl  (5列特征表)
  suvtk.features_output/reoriented_nucleotide_sequences.fna  (核酸序列)

输出:
  virus-annotations/{contig_id}.gb

用法:
  python tbl2gb.py \
    --tbl analyze_hypothetical/featuretable_updated.tbl \
    --fasta suvtk.features_output/reoriented_nucleotide_sequences.fna \
    -o virus-annotations/
"""

import argparse, os, sys, re
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def parse_feature_table(tbl_path):
    """
    解析 5 列特征表，返回 {contig_id: [features]}
    格式: >Feature contig_id  (或 contig_id)
           start  end  type
                qualifier  value
    """
    features_by_contig = defaultdict(list)
    current_contig = None
    current_feat = None

    with open(tbl_path, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line.strip():
                continue

            # 新特征行: start<TAB>end<TAB>type
            if line[0].isdigit() or (line.startswith('<') and len(line.split('\t')) >= 3):
                if current_feat and current_contig:
                    features_by_contig[current_contig].append(current_feat)

                parts = line.split('\t')
                if len(parts) >= 3:
                    try:
                        # 处理 <1 或 >N 格式
                        start_str = parts[0].lstrip('<>')
                        end_str = parts[1].lstrip('<>')
                        start = int(start_str) - 1  # 5列格式是1-based → 0-based
                        end = int(end_str)
                        feat_type = parts[2].strip()
                        current_feat = {
                            'location': (start, end),
                            'type': feat_type,
                            'qualifiers': defaultdict(list),
                        }
                    except ValueError:
                        current_feat = None
                continue

            # 新 contig 标记
            if line.startswith('>Feature '):
                contig_name = line.split('>Feature ', 1)[1].strip()
                if current_contig is None or contig_name != current_contig:
                    if current_feat and current_contig:
                        features_by_contig[current_contig].append(current_feat)
                        current_feat = None
                    current_contig = contig_name
                continue

            # qualifier 行: \t\tqualifier_name\tvalue
            if line.startswith('\t\t') and current_feat:
                parts = line.split('\t')
                # 去除开头的空字符串
                qual_parts = [p for p in parts if p]
                if len(qual_parts) >= 2:
                    qual_name = qual_parts[0].strip()
                    qual_value = qual_parts[1].strip() if len(qual_parts) > 1 else ''
                    current_feat['qualifiers'][qual_name].append(qual_value)
                continue

        # 最后一条
        if current_feat and current_contig:
            features_by_contig[current_contig].append(current_feat)

    return features_by_contig


def build_gb_record(contig_id, seq, features, taxonomy_info=None):
    """从序列和特征列表构建 GenBank SeqRecord"""
    record = SeqRecord(
        seq,
        id=contig_id,
        name=contig_id,
        description=taxonomy_info.get(contig_id, '') if taxonomy_info else '',
        annotations={
            'molecule_type': 'RNA',
            'topology': 'linear',
        }
    )

    # 添加 source feature
    source_feat = SeqFeature(
        FeatureLocation(0, len(seq)),
        type='source',
        qualifiers={
            'mol_type': ['genomic RNA'],
            'organism': [taxonomy_info.get(contig_id, 'Unknown virus').replace('_', ' ')] if taxonomy_info else ['Unknown virus'],
        }
    )
    record.features.append(source_feat)

    for feat_data in features:
        start, end = feat_data['location']
        location = FeatureLocation(start, end, strand=1)  # assume + strand
        qualifiers = {k: v for k, v in feat_data['qualifiers'].items()}
        sf = SeqFeature(location, type=feat_data['type'], qualifiers=qualifiers)
        record.features.append(sf)

    return record


def load_taxonomy(taxonomy_tsv):
    """从 taxonomy.tsv 加载 contig→taxonomy 映射"""
    tax_map = {}
    if not taxonomy_tsv or not Path(taxonomy_tsv).exists():
        return tax_map
    with open(taxonomy_tsv) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                tax_map[parts[0]] = parts[1]
    return tax_map


def main():
    parser = argparse.ArgumentParser(description="将 featuretable + FASTA 转为 GenBank (.gb)")
    parser.add_argument('--tbl', required=True, help='featuretable_updated.tbl 路径')
    parser.add_argument('--fasta', required=True, help='FASTA 序列文件 (reoriented_nucleotide_sequences.fna)')
    parser.add_argument('-o', '--outdir', required=True, help='输出目录 (.gb 文件)')
    parser.add_argument('--taxonomy', help='taxonomy.tsv 路径 (可选, 用于 organism 注释)')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 解析特征表
    print(f"解析特征表: {args.tbl}")
    features_by_contig = parse_feature_table(args.tbl)
    print(f"  找到 {len(features_by_contig)} 个 contig, 共 {sum(len(v) for v in features_by_contig.values())} 个 CDS")

    # 加载序列
    print(f"加载序列: {args.fasta}")
    seqs = {}
    for rec in SeqIO.parse(args.fasta, 'fasta'):
        seqs[rec.id] = rec.seq

    # 加载 taxonomy
    taxonomy = load_taxonomy(args.taxonomy) if args.taxonomy else {}
    if taxonomy:
        print(f"加载 taxonomy: {len(taxonomy)} 条")

    # 生成 .gb 文件
    count = 0
    for contig_id, features in features_by_contig.items():
        if contig_id not in seqs:
            print(f"  警告: {contig_id} 的序列未找到，跳过")
            continue

        seq = seqs[contig_id]
        record = build_gb_record(contig_id, seq, features, taxonomy)

        # 清理 contig_id 中的非法文件名字符
        safe_name = re.sub(r'[^\w\-\.]', '_', contig_id)
        out_path = outdir / f"{safe_name}.gb"
        SeqIO.write(record, out_path, 'genbank')
        count += 1

    print(f"\n完成: {count} 个 .gb 文件 → {outdir}/")


if __name__ == '__main__':
    main()
