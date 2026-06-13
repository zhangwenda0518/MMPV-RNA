#!/usr/bin/env python
"""
C9: Classify contigs by host category using probability lookup tables.

Classification logic (cascading, most specific first):
  1. Species level lookup → species_host_probability.tsv
  2. Genus level fallback → genus_host_probability.tsv
  3. Family level fallback → family_host_probability.tsv
  4. Unknown if no match at any level

Outputs:
  - Per-category FASTA: {output_dir}/fasta/{Category}/contigs.fasta
  - Per-category TSV:  {output_dir}/classified/{Category}.tsv
  - Overall summary:   {output_dir}/classification_summary.tsv
  - Confidence report: {output_dir}/confidence_report.tsv

Usage:
  python C9_classify_contigs.py -i final_integrated_classification.tsv -f sequences.fasta
  python C9_classify_contigs.py -i ACVirus_results/step4_classification_eval.acvirus/final_result.tsv -f contigs.fasta -o classified/
"""
import polars as pl
import os, argparse, csv
from collections import defaultdict

CORE_CATS = ["Plant", "Insecta", "Arachnida", "Fungi",
             "Bacteria", "Animal_other", "Aves", "Oomycetes",
             "Protist", "Archaea", "Human"]


def load_prob_table(path, key_col="Family"):
    """Load a probability lookup table, return {taxon: {predicted_host, confidence, ...}}."""
    if not os.path.exists(path):
        return {}
    lookup = {}
    with open(path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            taxon = row.get(key_col, '').strip()
            if not taxon:
                continue
            lookup[taxon] = {
                'Predicted_Host': row.get('Predicted_Host', 'Unknown'),
                'Confidence_Level': row.get('Confidence_Level', 'Unknown'),
                'Integrated_Confidence': float(row.get('Integrated_Confidence', 0)),
                'Total_Records': int(row.get('Total_Records', 0)),
                'P_Max': float(row.get('P_Max', 0)),
            }
    return lookup


def load_contigs(path):
    """Load contig classification TSV (tab-separated)."""
    df = pl.read_csv(path, separator='\t', truncate_ragged_lines=True, ignore_errors=True)
    return df


def classify_contig(row, fam_lookup, gen_lookup, sp_lookup, ord_lookup):
    """
    Robust classification with column-shift tolerance.
    Cascade: Species → Genus → Family → Order (most specific first).
    Tries exact match first, then shifted match for column misalignment.
    """
    order = str(row.get('Order', '')).strip() if row.get('Order') else ''
    family = str(row.get('Family', '')).strip() if row.get('Family') else ''
    genus = str(row.get('Genus', '')).strip() if row.get('Genus') else ''
    species = str(row.get('Species', '')).strip() if row.get('Species') else ''

    result = {
        'Predicted_Host': 'Unknown',
        'Confidence_Level': 'Unknown',
        'Integrated_Confidence': 0.0,
        'Determination_Level': 'None',
    }

    candidates = []

    # Exact matches
    if species and species in sp_lookup:
        candidates.append(('Species', 'Species', sp_lookup[species]))
    if genus and genus in gen_lookup:
        candidates.append(('Genus', 'Genus', gen_lookup[genus]))
    if family and family in fam_lookup:
        candidates.append(('Family', 'Family', fam_lookup[family]))
    if order and order in ord_lookup:
        candidates.append(('Order', 'Order', ord_lookup[order]))

    # Shifted matches (column name ≠ actual level)
    if species and species in gen_lookup:
        candidates.append(('Species', 'Genus*', gen_lookup[species]))
    if species and species in fam_lookup:
        candidates.append(('Species', 'Family*', fam_lookup[species]))
    if species and species in ord_lookup:
        candidates.append(('Species', 'Order*', ord_lookup[species]))
    if genus and genus in sp_lookup:
        candidates.append(('Genus', 'Species*', sp_lookup[genus]))
    if genus and genus in fam_lookup:
        candidates.append(('Genus', 'Family*', fam_lookup[genus]))
    if genus and genus in ord_lookup:
        candidates.append(('Genus', 'Order*', ord_lookup[genus]))
    if family and family in gen_lookup:
        candidates.append(('Family', 'Genus*', gen_lookup[family]))
    if family and family in sp_lookup:
        candidates.append(('Family', 'Species*', sp_lookup[family]))
    if family and family in ord_lookup:
        candidates.append(('Family', 'Order*', ord_lookup[family]))
    if order and order in fam_lookup:
        candidates.append(('Order', 'Family*', fam_lookup[order]))
    if order and order in gen_lookup:
        candidates.append(('Order', 'Genus*', gen_lookup[order]))
    if order and order in sp_lookup:
        candidates.append(('Order', 'Species*', sp_lookup[order]))

    if candidates:
        # Sort by: actual level priority (Species=4, Genus=3, Family=2, Order=1), then confidence
        level_rank = {'Species': 4, 'Genus': 3, 'Family': 2, 'Order': 1}
        candidates.sort(key=lambda x: (
            level_rank.get(x[1].rstrip('*'), 0),
            x[2].get('Integrated_Confidence', 0)
        ), reverse=True)

        col_name, real_level, info = candidates[0]
        result['Predicted_Host'] = info['Predicted_Host']
        result['Confidence_Level'] = info['Confidence_Level']
        result['Integrated_Confidence'] = info['Integrated_Confidence']
        result['Determination_Level'] = f"{real_level}(via {col_name})"

    return result


def load_fasta_sequences(fasta_path):
    """Load contig sequences from FASTA. Returns {contig_id: sequence}."""
    if not fasta_path or not os.path.exists(fasta_path):
        return {}
    seqs = {}
    current_id = None
    current_seq = []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    seqs[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]  # first word after >
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            seqs[current_id] = ''.join(current_seq)
    return seqs


def main():
    p = argparse.ArgumentParser(description="C9: Classify contigs by host category")
    p.add_argument("-i", "--input", required=True,
                   help="Contig classification TSV (final_integrated_classification.tsv or ACVirus final_result.tsv)")
    p.add_argument("-f", "--fasta", default=None,
                   help="FASTA file with contig sequences (optional)")
    p.add_argument("-o", "--output_dir", default="classified/",
                   help="Output directory for classified .fasta and .tsv files")
    p.add_argument("--prob_dir", default="cross_analysis/",
                   help="Directory containing C7 output (*_host_probability.tsv)")
    p.add_argument("--mode", default="medium", choices=["all", "medium", "high"],
                   help="Output mode: all, medium (High+Medium), high (High only)")
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load probability tables (C7 output)
    fam_lookup = load_prob_table(os.path.join(args.prob_dir, "family_host_probability.tsv"), "Family")
    gen_lookup = load_prob_table(os.path.join(args.prob_dir, "genus_host_probability.tsv"), "Genus")
    sp_lookup = load_prob_table(os.path.join(args.prob_dir, "species_host_probability.tsv"), "Species")
    ord_lookup = load_prob_table(os.path.join(args.prob_dir, "order_host_probability.tsv"), "Order")
    print(f"  Order: {len(ord_lookup)}, Family: {len(fam_lookup)}, Genus: {len(gen_lookup)}, Species: {len(sp_lookup)}")

    # Load contigs
    print(f"Loading contigs: {args.input}")
    df = load_contigs(args.input)
    print(f"  {df.height} contigs")

    # Load FASTA
    fasta_seqs = {}
    if args.fasta:
        print(f"Loading FASTA: {args.fasta}")
        fasta_seqs = load_fasta_sequences(args.fasta)
        print(f"  {len(fasta_seqs)} sequences")

    # Classify each contig
    print("Classifying...")
    results = []
    # Detect contig_id column
    cols = df.columns
    cid_col = 'Nucleotide' if 'Nucleotide' in cols else cols[0]

    for r in df.iter_rows(named=True):
        cid = str(r.get(cid_col, ''))
        cls = classify_contig(r, fam_lookup, gen_lookup, sp_lookup, ord_lookup)
        cls['contig_id'] = cid
        results.append(cls)

    result_df = pl.DataFrame(results)

    # Merge with original classification
    out_df = df.with_columns([
        pl.Series('Predicted_Host', [r['Predicted_Host'] for r in results]),
        pl.Series('Confidence_Level', [r['Confidence_Level'] for r in results]),
        pl.Series('Integrated_Confidence', [r['Integrated_Confidence'] for r in results]),
        pl.Series('Determination_Level', [r['Determination_Level'] for r in results]),
    ])

    # Filter by mode
    conf_rank = {"High": 3, "Medium": 2, "Low (Singleton/Rare)": 1, "Low (Shared Genus)": 1, "Unknown": 0}
    min_rank = {"all": 0, "medium": 2, "high": 3}[args.mode]
    mask = pl.Series([conf_rank.get(str(l), 0) >= min_rank
                     for l in out_df['Confidence_Level'].to_list()])
    out_df = out_df.with_columns(mask.alias('_pass_conf'))
    print(f"  Mode={args.mode}: {mask.sum()} / {df.height} pass")

    # Save per-category TSV and FASTA (flat output, one directory)
    cat_summary = {}
    for cat in sorted(set(out_df['Predicted_Host'].to_list())):
        cat_df = out_df.filter(pl.col('Predicted_Host') == cat)
        if cat_df.height == 0:
            continue

        high_n = cat_df.filter(pl.col('_pass_conf')).height
        low_n = cat_df.height - high_n
        cat_summary[cat] = {'total': cat_df.height, 'classified': high_n, 'low_confidence': low_n}

        # TSV: {Category}.classified.tsv
        tsv_path = os.path.join(args.output_dir, f"{cat}.classified.tsv")
        cat_df.write_csv(tsv_path, separator='\t')

        # FASTA: {Category}.classified.fasta (only high/medium confidence)
        if fasta_seqs:
            pass_df = cat_df.filter(pl.col('_pass_conf'))
            fasta_path = os.path.join(args.output_dir, f"{cat}.classified.fasta")
            written = 0
            with open(fasta_path, 'w') as f:
                for row in pass_df.iter_rows(named=True):
                    cid = str(row['contig_id'])
                    if cid in fasta_seqs:
                        f.write(f">{cid}\n{fasta_seqs[cid]}\n")
                        written += 1
            print(f"  {cat}: {cat_df.height} contigs ({high_n} pass, {written} fasta)")

    # Save overall classified TSV
    out_df.write_csv(os.path.join(args.output_dir, "classification_result.tsv"), separator='\t')

    # Save summary
    with open(os.path.join(args.output_dir, "classification_summary.tsv"), 'w', newline='') as f:
        f.write("Category\tTotal_Contigs\tClassified\tLow_Confidence\tPercent_Classified\n")
        for cat in sorted(cat_summary.keys()):
            s = cat_summary[cat]
            pct = s['classified'] / s['total'] * 100 if s['total'] > 0 else 0
            f.write(f"{cat}\t{s['total']}\t{s['classified']}\t{s['low_confidence']}\t{pct:.1f}%\n")
        total_all = sum(s['total'] for s in cat_summary.values())
        total_cls = sum(s['classified'] for s in cat_summary.values())
        f.write(f"TOTAL\t{total_all}\t{total_cls}\t{total_all-total_cls}\t"
                f"{total_cls/total_all*100:.1f}%\n")

    # Confidence report
    with open(os.path.join(args.output_dir, "confidence_report.tsv"), 'w', newline='') as f:
        f.write("Confidence_Level\tCount\tPercent\n")
        for lvl in ["High", "Medium", "Low (Singleton/Rare)", "Low (Shared Genus)", "Unknown"]:
            cnt = out_df.filter(pl.col('Confidence_Level') == lvl).height
            pct = cnt / out_df.height * 100 if out_df.height > 0 else 0
            f.write(f"{lvl}\t{cnt}\t{pct:.1f}%\n")

    print(f"\nDone. Output in: {args.output_dir}/")


if __name__ == "__main__":
    main()
