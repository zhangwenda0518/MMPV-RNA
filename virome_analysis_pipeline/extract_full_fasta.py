#!/usr/bin/env python3
"""
extract_full_fasta.py
=====================
Two modes:
  1. --fill : N-fill with reference (RNA: fill only if N<5%; DNA: always fill)
     Requires --ref_info and --ref_dir
  2. Default: filter-only (--max_n, --min_len), no filling
"""

import argparse
import re
from pathlib import Path

try:
    from Bio import Align
    from Bio import SeqIO
    HAS_BIO = True
except ImportError:
    HAS_BIO = False


def get_longest_sequence(fasta_path):
    """Return (seq_str, n_ratio) of the longest contig in a FASTA."""
    longest_seq = ""
    current_seq = []

    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_seq:
                        seq_str = "".join(current_seq)
                        if len(seq_str) > len(longest_seq):
                            longest_seq = seq_str
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_seq:
                seq_str = "".join(current_seq)
                if len(seq_str) > len(longest_seq):
                    longest_seq = seq_str
    except Exception as e:
        print(f"读取文件错误 {fasta_path}: {e}")

    if longest_seq and len(longest_seq) > 0:
        n_ratio = (longest_seq.upper().count('N') / len(longest_seq)) * 100
        return longest_seq, n_ratio
    return "", 100.0


def extract_accession(dir_name):
    """Extract NCBI accession from dir name like 'Taxonomy_sp_lycii_OR489165.1'."""
    m = re.search(r'([A-Z]{1,2}_\d+\.\d+|[A-Z]{2}\d+\.\d+)$', dir_name)
    return m.group(1) if m else dir_name


def n_fill_with_reference(cons_seq, ref_fasta, logger=None):
    """
    Fill N positions in consensus using reference via global pairwise alignment.
    Returns (filled_seq, n_filled) or (cons_seq, 0) if no N or fill fails.
    """
    if not HAS_BIO:
        if logger: logger("Biopython unavailable, skip N-fill")
        return cons_seq, 0

    n_before = cons_seq.upper().count('N')
    if n_before == 0:
        return cons_seq, 0

    try:
        ref_rec = list(SeqIO.parse(ref_fasta, "fasta"))
        ref_seq = str(ref_rec[0].seq).upper()
    except Exception as e:
        if logger: logger(f"Failed to read reference {ref_fasta}: {e}")
        return cons_seq, 0

    try:
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -3
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -2

        alignments = aligner.align(ref_seq, cons_seq)
        if not alignments:
            return cons_seq, 0

        best = alignments[0]
        ref_aln = str(best[0])
        cons_aln = str(best[1])

        filled_chars = list(cons_seq)
        cons_pos = 0
        filled_n = 0
        for i in range(len(cons_aln)):
            r_base = ref_aln[i]
            c_base = cons_aln[i]
            if c_base == '-':
                continue
            if c_base == 'N' and r_base != '-' and r_base in 'ACGT':
                filled_chars[cons_pos] = r_base
                filled_n += 1
            cons_pos += 1

        return ''.join(filled_chars), filled_n
    except Exception as e:
        if logger: logger(f"N-fill alignment failed: {e}")
        return cons_seq, 0


def load_molecule_types(ref_info_path):
    """Return dict {accession: is_rna (True/False)} from ref_info TSV.
    Reads Molecule_Type2 / Molecule_type columns. Default: RNA if not found."""
    mol_map = {}
    if not ref_info_path or not Path(ref_info_path).exists():
        return mol_map
    with open(ref_info_path) as f:
        header = f.readline().strip().split('\t')
        try:
            idx_acc = header.index('Accession')
        except ValueError:
            return mol_map
        idx_mol = None
        for col in ['Molecule_Type2', 'Molecule_type']:
            if col in header:
                idx_mol = header.index(col); break
        if idx_mol is None:
            return mol_map
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) > idx_mol and cols[idx_mol]:
                mol_map[cols[idx_acc]] = 'RNA' in cols[idx_mol].upper()
    return mol_map


def main():
    parser = argparse.ArgumentParser(description="提取最长 FASTA 序列并重命名，支持参考序列 N 填补")
    parser.add_argument("-d", "--dir", required=True, help="输入的组装总目录")
    parser.add_argument("-o", "--outdir", required=True, help="输出目录")
    parser.add_argument("--target_file", default="11.Ultimate_Circular_Result.fasta", help="目标文件名")

    # Fill mode
    parser.add_argument("--fill", action="store_true", help="启用参考序列 N 填补")
    parser.add_argument("--ref_info", help="ref_info TSV (获取 Molecule_Type2, 默认RNA)")
    parser.add_argument("--ref_dir", help="参考 FASTA 目录 (如 2_Virus_variants_Results/, fill 模式必需)")
    parser.add_argument("--max_n_genome", type=float, default=5.0, help="RNA 病毒填补前全长 N 上限%% (默认: 5)")
    parser.add_argument("--no_dna_fill", action="store_true", help="DNA 病毒不填补 (默认: DNA 始终填补)")

    # Filter mode (no --fill)
    parser.add_argument("--max_n", type=float, default=5.0, help="过滤模式 N 含量阈值%% (默认: 5)")
    parser.add_argument("--min_len", type=int, default=150, help="最短序列长度 bp (默认: 150)")

    # Plot
    parser.add_argument("--plot", action="store_true", help="生成组装结果统计图")

    args = parser.parse_args()

    in_base = Path(args.dir)
    out_base = Path(args.outdir)

    if not in_base.exists():
        print(f"❌ 找不到输入目录: {in_base}")
        return

    out_base.mkdir(parents=True, exist_ok=True)

    # Load molecule type for fill mode
    mol_map = {}
    if args.fill:
        if not args.ref_dir:
            print("❌ --fill 模式需要 --ref_dir")
            return
        if args.ref_info:
            mol_map = load_molecule_types(args.ref_info)
            print(f"📋 加载 {len(mol_map)} 条分子类型 (默认RNA)")

    if args.fill:
        print(f"🔧 N 填补模式: RNA 病毒 N<{args.max_n_genome}% 才填补, DNA 病毒始终填补")
        ref_base = Path(args.ref_dir)
    else:
        print(f"📏 过滤模式: N > {args.max_n}% | 长度 < {args.min_len}bp → 跳过")

    print(f"🔍 扫描 {in_base} → {out_base}\n")

    success_count = 0
    missing_count = 0
    empty_count = 0
    n_fail_count = 0
    short_fail_count = 0
    n_filled_count = 0

    for tax_dir in sorted(in_base.iterdir()):
        if not tax_dir.is_dir():
            continue

        out_tax_dir = out_base / tax_dir.name
        out_tax_dir.mkdir(parents=True, exist_ok=True)

        acc = extract_accession(tax_dir.name)

        for sample_dir in sorted(tax_dir.iterdir()):
            if not sample_dir.is_dir():
                continue

            target_files = list(sample_dir.rglob(args.target_file))
            if not target_files:
                missing_count += 1
                continue

            for target_file in target_files:
                sample_accession = sample_dir.name
                file_prefix = sample_accession.replace('_', '.')
                out_name = f"{file_prefix}.full.fasta"
                out_path = out_tax_dir / out_name

                longest_seq, n_ratio = get_longest_sequence(target_file)

                if not longest_seq:
                    empty_count += 1
                    continue

                seq_len = len(longest_seq)

                # ---------- FILL MODE ----------
                if args.fill:
                    is_RNA = mol_map.get(acc, mol_map.get(acc.split('.')[0]))
                    # Default: RNA if not found
                    if is_RNA is None:
                        is_RNA = True

                    if is_RNA:
                        # RNA: fill only if N < max_n_genome
                        if n_ratio >= args.max_n_genome:
                            n_fail_count += 1
                            print(f"⏭️ RNA N超标不填补 (N={n_ratio:.1f}% ≥ {args.max_n_genome}%): {tax_dir.name}/{out_name}")
                            continue
                    else:
                        # DNA: always fill (unless --no_dna_fill)
                        if args.no_dna_fill:
                            pass  # skip fill

                    # Do N-fill
                    ref_fasta = ref_base / "virus-fasta" / f"ref_{acc}" / f"ref_{acc}.ref.fasta"
                    if ref_fasta.exists() and n_ratio > 0:
                        filled_seq, n_filled = n_fill_with_reference(
                            longest_seq, str(ref_fasta), logger=print
                        )
                        if n_filled > 0:
                            longest_seq = filled_seq
                            n_after = (longest_seq.upper().count('N') / len(longest_seq)) * 100
                            n_filled_count += 1
                            print(f"🔧 填补 {n_filled}N → {n_ratio:.1f}%→{n_after:.2f}%: {tax_dir.name}/{out_name}")
                        else:
                            print(f"✅ (无N或填补失败) {tax_dir.name}/{out_name}  ({seq_len} bp, N={n_ratio:.1f}%)")
                    else:
                        print(f"✅ (无需填补) {tax_dir.name}/{out_name}  ({seq_len} bp, N={n_ratio:.1f}%)")

                    try:
                        with open(out_path, 'w') as fout:
                            fout.write(f">{sample_accession}\n")
                            fout.write(f"{longest_seq}\n")
                        success_count += 1
                    except Exception as e:
                        print(f"⚠️ 写入 {out_path} 时出错: {e}")
                    continue

                # ---------- FILTER MODE ----------
                if seq_len < args.min_len:
                    short_fail_count += 1
                    print(f"⏭️ 太短 ({seq_len}bp < {args.min_len}): {tax_dir.name}/{out_name}")
                    continue

                if n_ratio > args.max_n:
                    n_fail_count += 1
                    print(f"⏭️ N超标 (N={n_ratio:.1f}% > {args.max_n}%): {tax_dir.name}/{out_name}")
                    continue

                try:
                    with open(out_path, 'w') as fout:
                        fout.write(f">{sample_accession}\n")
                        fout.write(f"{longest_seq}\n")
                    success_count += 1
                    print(f"✅ {tax_dir.name}/{out_name}  ({seq_len} bp, N={n_ratio:.1f}%)")
                except Exception as e:
                    print(f"⚠️ 写入 {out_path} 时出错: {e}")

    print("-" * 60)
    if args.fill:
        print(f"提取+填补完成！成功: {success_count} | 填补: {n_filled_count} | N超标跳过: {n_fail_count} | 缺失: {missing_count} | 空: {empty_count}")
    else:
        print(f"提取完成！成功: {success_count} | N超标: {n_fail_count} | 太短: {short_fail_count} | 缺失: {missing_count} | 空: {empty_count}")
    print(f"{out_base.absolute()}")

    # ── Assembly stats plot ──
    if args.plot and success_count > 0:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np

        # Collect all output contig lengths
        lengths = []
        labels = []
        for f in sorted(out_base.rglob("*.full.fasta")):
            seq, _ = get_longest_sequence(f)
            if seq:
                lengths.append(len(seq))
                labels.append(f.parent.name.split('_')[-1][:12] if '_' in f.parent.name else f.stem[:12])

        if lengths:
            lengths = np.array(sorted(lengths, reverse=True))
            n = len(lengths)
            total_bp = sum(lengths)
            cumsum = np.cumsum(lengths)
            n50_idx = np.searchsorted(cumsum, total_bp / 2)
            n90_idx = np.searchsorted(cumsum, total_bp * 0.9)
            n50 = lengths[min(n50_idx, n - 1)]
            n90 = lengths[min(n90_idx, n - 1)]

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

            # Histogram
            ax1.hist(lengths, bins=min(20, n), color='#1f77b4', edgecolor='white', alpha=0.8)
            for val, name, color in [(n50, 'N50', '#d62728'), (n90, 'N90', '#ff7f0e')]:
                ax1.axvline(val, color=color, linestyle='--', linewidth=2, label=f'{name}={val:,}bp')
            ax1.set_xlabel('Contig Length (bp)', fontweight='bold')
            ax1.set_ylabel('Count', fontweight='bold')
            ax1.set_title(f'Assembly Contig Length Distribution\n(n={n}, total={total_bp/1e6:.2f}Mb)',
                         fontweight='bold')
            ax1.legend()

            # Ranked lengths
            ax2.barh(range(n), lengths, color=plt.cm.viridis(np.linspace(0.2, 0.9, n)), edgecolor='#333')
            ax2.axvline(n50, color='#d62728', linestyle='--', linewidth=2, label=f'N50={n50:,}bp')
            ax2.axvline(n90, color='#ff7f0e', linestyle='--', linewidth=2, label=f'N90={n90:,}bp')
            ax2.set_xlabel('Length (bp)', fontweight='bold')
            ax2.set_title('Contigs Ranked by Length', fontweight='bold')
            ax2.legend(loc='lower right')

            plt.tight_layout()
            plot_path = out_base / "assembly_stats.pdf"
            plot_path_png = out_base / "assembly_stats.png"
            fig.savefig(plot_path, dpi=300, bbox_inches='tight')
            fig.savefig(plot_path_png, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Assembly plot -> {plot_path}")

            # Stats table
            with open(out_base / "assembly_stats.csv", 'w') as sf:
                sf.write(f"n_contigs,total_bp,longest,N50,N90\n")
                sf.write(f"{n},{total_bp},{lengths[0]},{n50},{n90}\n")


if __name__ == "__main__":
    main()
