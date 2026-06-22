#!/usr/bin/env python3
"""
将 virome_discovery_pipeline 输出的 centroids + taxonomy
转换为 virome_analysis_pipeline 可识别的 --reference + --ref_info 格式。

Usage:
    python discovery2analysis.py \
        --centroids 04_CLUSTER/centroids/final_centroids.fasta \
        --taxonomy 05_Taxonomy/integrated/final_integrated_classification.tsv \
        --output_prefix my_project

输出:
    my_project.reference.fasta    # 参考基因组 FASTA
    my_project.ref_info.tsv       # 参考信息表

然后可直接喂给 virome_analysis_pipeline:
    python auto_known_virus.py \
        --reference my_project.reference.fasta \
        --ref_info my_project.ref_info.tsv \
        --reads_dir 00b_HostDepletion/ \
        --output_dir virus_analysis/
"""

import argparse
import hashlib
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="discovery → analysis pipeline 格式转换"
    )
    parser.add_argument("--centroids", required=True,
                        help="发现管线输出的 centroids FASTA (如 04_CLUSTER/centroids/final_centroids.fasta)")
    parser.add_argument("--taxonomy", required=True,
                        help="发现管线输出的 taxonomy TSV (如 05_Taxonomy/integrated/final_integrated_classification.tsv)")
    parser.add_argument("--output_prefix", required=True,
                        help="输出文件前缀 (生成 <prefix>.reference.fasta 和 <prefix>.ref_info.tsv)")
    parser.add_argument("--min_length", type=int, default=500,
                        help="最短序列长度 (default: 500)")
    parser.add_argument("--min_completeness", type=float, default=0.0,
                        help="最低分类完整度 0-1 (default: 0, 即不过滤)")
    parser.add_argument("--molecule_type", default="RNA",
                        help="分子类型 (default: RNA)")
    parser.add_argument("--molecule_type2", default="ssRNA",
                        help="分子子类型, 用于 DNA/RNA 双轨过滤 (default: ssRNA)")
    parser.add_argument("--taxid_offset", type=int, default=9000000,
                        help="人工 TaxID 起始偏移量 (default: 9000000)")
    parser.add_argument("--segment_col", default="genome",
                        help="Segment 列默认值 (default: genome)")
    return parser.parse_args()


def load_taxonomy(taxonomy_path):
    """加载分类 TSV, 返回 {contig_id: {Realm..Species, completeness}}"""
    tax_map = {}
    with open(taxonomy_path, 'r', encoding='utf-8') as f:
        header = None
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if header is None:
                header = [h.strip() for h in parts]
                continue

            row = dict(zip(header, [p.strip() for p in parts]))
            contig_id = row.get('contig_id', '')
            if not contig_id:
                continue

            completeness = float(row.get('completeness', 0) or 0)
            tax_map[contig_id] = {
                'Realm': row.get('Realm', ''),
                'Kingdom': row.get('Kingdom', ''),
                'Phylum': row.get('Phylum', ''),
                'Class': row.get('Class', ''),
                'Order': row.get('Order', ''),
                'Family': row.get('Family', ''),
                'Genus': row.get('Genus', ''),
                'Species': row.get('Species', ''),
                'completeness': completeness,
            }
    return tax_map


def resolve_species(tax_entry):
    """从分类层级中取最精细的非空分类名作为 Species 标签"""
    rank_priority = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom', 'Realm']
    for rank in rank_priority:
        val = tax_entry.get(rank, '')
        if val and val.lower() not in ('na', 'nan', 'none', 'unclassified', 'unknown', ''):
            return val
    return 'Unannotated'


def generate_taxid(contig_id, offset):
    """基于 contig_id 生成稳定的人工 TaxID"""
    h = hashlib.md5(contig_id.encode()).hexdigest()[:6]
    return str(offset + int(h, 16) % 900000)


def load_centroids(fasta_path, min_length):
    """读取 centroids FASTA, 返回 [(header, seq, length)]"""
    seqs = []
    current_header, current_seq = None, []
    with open(fasta_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                if current_header and current_seq:
                    seq_str = ''.join(current_seq)
                    if len(seq_str) >= min_length:
                        seqs.append((current_header, seq_str, len(seq_str)))
                current_header = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_header and current_seq:
            seq_str = ''.join(current_seq)
            if len(seq_str) >= min_length:
                seqs.append((current_header, seq_str, len(seq_str)))
    return seqs


def main():
    args = parse_args()

    # 验证输入
    for fpath, label in [(args.centroids, "centroids"), (args.taxonomy, "taxonomy")]:
        if not os.path.exists(fpath):
            print(f"[ERROR] Cannot find {label} file: {fpath}", file=sys.stderr)
            sys.exit(1)

    print(f"[INFO] Loading taxonomy: {args.taxonomy}")
    tax_map = load_taxonomy(args.taxonomy)
    print(f"       Loaded {len(tax_map)} records")

    print(f"[INFO] Loading centroids: {args.centroids}")
    centroids = load_centroids(args.centroids, args.min_length)
    print(f"       Loaded {len(centroids)} sequences (min_length={args.min_length})")

    # 输出路径
    ref_fasta = f"{args.output_prefix}.reference.fasta"
    ref_info = f"{args.output_prefix}.ref_info.tsv"

    # 写入
    matched, unmatched, skipped = 0, 0, 0
    with open(ref_fasta, 'w', encoding='utf-8') as f_fa, \
         open(ref_info, 'w', encoding='utf-8') as f_info:

        # ref_info 表头
        f_info.write('\t'.join([
            'Accession', 'Taxid', 'Species_NCBI', 'Species_ICTV',
            'Segment', 'Molecule_type', 'Molecule_Type2',
            'Genome_Length', 'Taxonomic_Completeness',
            'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'
        ]) + '\n')

        for header, seq, length in centroids:
            tax = tax_map.get(header)

            # 分类完整度过滤
            if tax and args.min_completeness > 0:
                if tax['completeness'] < args.min_completeness:
                    skipped += 1
                    continue

            if tax:
                matched += 1
                species_label = resolve_species(tax)
                completeness = str(tax['completeness'])
                realm, kingdom = tax['Realm'], tax['Kingdom']
                phylum, klass = tax['Phylum'], tax['Class']
                order, family = tax['Order'], tax['Family']
                genus, species = tax['Genus'], tax['Species']
            else:
                unmatched += 1
                species_label = 'Unannotated'
                completeness = '0'
                realm = kingdom = phylum = klass = ''
                order = family = genus = species = ''

            taxid = generate_taxid(header, args.taxid_offset)

            f_info.write('\t'.join([
                header,                          # Accession
                taxid,                           # Taxid (artificial)
                species_label,                   # Species_NCBI
                species_label,                   # Species_ICTV
                args.segment_col,                # Segment
                args.molecule_type,              # Molecule_type
                args.molecule_type2,             # Molecule_Type2
                str(length),                     # Genome_Length
                completeness,                    # Taxonomic_Completeness
                realm, kingdom, phylum, klass, order, family, genus, species
            ]) + '\n')

            f_fa.write(f">{header}\n{seq}\n")

    print(f"\n[DONE] Conversion complete:")
    print(f"   Reference FASTA: {ref_fasta} ({matched + unmatched} sequences)")
    print(f"   Reference info:  {ref_info}")
    print(f"   Matched taxonomy: {matched}")
    print(f"   Unmatched:        {unmatched}")
    if skipped:
        print(f"   Filtered (completeness): {skipped}")
    print(f"\nUsage with virome_analysis_pipeline:")
    print(f"   python virome_analysis_pipeline/auto_known_virus.py \\")
    print(f"       --reference {ref_fasta} \\")
    print(f"       --ref_info {ref_info} \\")
    print(f"       --reads_dir <00b_HostDepletion/> \\")
    print(f"       --output_dir <output/>")


if __name__ == "__main__":
    main()
