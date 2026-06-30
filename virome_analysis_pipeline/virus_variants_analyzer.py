#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Unified Viral/Viroid Variant Analyzer (Peer-Reviewed Edition)
Highlights: Pool-Seq π, FORMAT-aware VCF Parser, Bulletproof Scalers
"""

import os
import glob
import re
import math
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import seaborn as sns
from pathlib import Path
import logging

try:
    from Bio import Entrez, SeqIO
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
sns.set_theme(style="whitegrid", rc={"axes.edgecolor": "black"})

class UnifiedVariantAnalyzer:
    def __init__(self, input_dir: str, output_dir: str, min_depth: int = 50, min_af: float = 0.05, virus_name: str = "", acc_id: str = ""):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.min_depth = min_depth; self.min_af = min_af
        self.master_df = None; self.popgen_df = None; self.total_samples = 0
        self.virus_name = virus_name if virus_name else self.input_dir.name.replace('_', ' ')
        self.acc_id = acc_id; self.gene_features = []; self.genome_length = 0

    def _fetch_ncbi_annotation(self):
        if not HAS_BIOPYTHON: return
        Entrez.email = os.environ.get("NCBI_EMAIL", "researcher@example.com")
        if not self.acc_id:
            match = re.search(r'\b[A-Z]{1,2}_?[0-9]{5,8}(?:\.[0-9]+)?\b', self.virus_name)
            if match: self.acc_id = match.group(0)
        if not self.acc_id: return
        logging.info(f"🌐 [NCBI] Accession: {self.acc_id}")
        try:
            handle = Entrez.efetch(db="nucleotide", id=self.acc_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank"); self.genome_length = len(record.seq)
            for f in record.features:
                if f.type in ['CDS', 'gene', 'mat_peptide']:
                    gname = f.qualifiers.get('gene', [''])[0] or f.qualifiers.get('product', ['Unnamed'])[0]
                    if gname != 'Unnamed':
                        self.gene_features.append({'name': gname, 'start': int(f.location.start),
                            'end': int(f.location.end), 'strand': f.location.strand})
            logging.info(f"✅ [NCBI] genome={self.genome_length}bp, genes={len(self.gene_features)}")
        except Exception as e:
            logging.warning(f"⚠️ [NCBI] {e}")

    @staticmethod
    def classify_molecular(ref: str, alt: str) -> str:
        ref, alt = str(ref).upper(), str(alt).upper()
        if len(ref) != len(alt) or '+' in alt or '-' in alt: return 'Indel'
        transitions = [{'A', 'G'}, {'C', 'T'}, {'C', 'U'}]
        if {ref, alt} in transitions: return 'Transition (Ts)'
        elif ref in list('ACGT') and alt in list('ACGT'): return 'Transversion (Tv)'
        return 'Complex/Unknown'

    @staticmethod
    def classify_functional(ref_aa: str, alt_aa: str) -> str:
        ref_aa, alt_aa = str(ref_aa), str(alt_aa)
        if ref_aa in ['nan', 'NA', 'None'] or alt_aa in ['nan', 'NA', 'None']: return 'Unannotated'
        return 'Synonymous' if ref_aa == alt_aa else 'Non-Synonymous'

    def _parse_tsv(self, filepath: Path, sample_name: str) -> pd.DataFrame:
        try:
            df = pd.read_csv(filepath, sep='\t')
            if 'POS' not in df.columns or 'ALT_FREQ' not in df.columns: return pd.DataFrame()
            df = df[(df['TOTAL_DP'] >= self.min_depth) & (df['ALT_FREQ'] >= self.min_af) & (df['PASS'] == True)].copy()
            if df.empty: return pd.DataFrame()
            df['Sample_ID'] = sample_name
            df['Molecular_Type'] = df.apply(lambda x: self.classify_molecular(x['REF'], x['ALT']), axis=1)
            if 'REF_AA' in df.columns and 'ALT_AA' in df.columns:
                df['Functional_Type'] = df.apply(lambda x: self.classify_functional(x['REF_AA'], x['ALT_AA']), axis=1)
            else: df['Functional_Type'] = 'Unannotated'
            return df[['Sample_ID', 'POS', 'REF', 'ALT', 'TOTAL_DP', 'ALT_FREQ', 'Molecular_Type', 'Functional_Type']]
        except Exception: return pd.DataFrame()

    def _parse_vcf(self, filepath: Path, sample_name: str) -> pd.DataFrame:
        """FORMAT-aware VCF parser with INFO fallback."""
        records = []
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    cols = line.strip().split('\t')
                    if len(cols) < 8: continue
                    ref = cols[3]; alt = cols[4].split(',')[0]
                    if cols[6] not in ['PASS', '.']: continue

                    dp_val, af_val = None, None

                    # 1) FORMAT / Genotype column (sample-specific DP/AF/AD)
                    if len(cols) >= 10:
                        fmt_keys = cols[8].split(':')
                        sample_vals = cols[9].split(':')
                        fmt_dict = dict(zip(fmt_keys, sample_vals))
                        if 'DP' in fmt_dict and fmt_dict['DP'] != '.':
                            try: dp_val = int(fmt_dict['DP'])
                            except ValueError: pass
                        if 'AF' in fmt_dict and fmt_dict['AF'] != '.':
                            try: af_val = float(fmt_dict['AF'].split(',')[0])
                            except ValueError: pass
                        elif 'AD' in fmt_dict and fmt_dict['AD'] != '.' and dp_val:
                            try:
                                ad_parts = fmt_dict['AD'].split(',')
                                if len(ad_parts) >= 2 and ad_parts[1].isdigit():
                                    af_val = float(ad_parts[1]) / dp_val
                            except Exception: pass

                    # 2) INFO field fallback
                    info = dict(item.split('=', 1) if '=' in item else (item, True) for item in cols[7].split(';'))
                    if dp_val is None and 'DP' in info and str(info['DP']).isdigit():
                        dp_val = int(info['DP'])
                    if af_val is None:
                        if 'AF' in info:
                            try: af_val = float(str(info['AF']).split(',')[0])
                            except ValueError: pass
                        elif 'AO' in info and 'DP' in info:
                            try:
                                ao = int(str(info['AO']).split(',')[0]); dp = int(info['DP'])
                                if dp > 0: af_val = ao / dp
                            except Exception: pass

                    dp_val = dp_val or 0; af_val = af_val or 0.0
                    if dp_val >= self.min_depth and af_val >= self.min_af:
                        records.append({'Sample_ID': sample_name, 'POS': int(cols[1]),
                            'REF': ref, 'ALT': alt, 'TOTAL_DP': dp_val, 'ALT_FREQ': af_val,
                            'Molecular_Type': self.classify_molecular(ref, alt), 'Functional_Type': 'Unannotated'})
            return pd.DataFrame(records)
        except Exception as e:
            logging.warning(f"⚠️ VCF parse error {filepath}: {e}")
            return pd.DataFrame()

    def load_and_merge(self):
        self._fetch_ncbi_annotation()
        logging.info(f"🚀 Scanning: {self.input_dir}")
        files = glob.glob(os.path.join(self.input_dir, "**", "*.variants.tsv"), recursive=True) + \
                glob.glob(os.path.join(self.input_dir, "**", "*.filtered.vcf"), recursive=True)
        all_data, parsed_samples = [], set()
        for f in set(files):
            sample = Path(f).name.split('.')[0].split('_')[0]
            if sample not in parsed_samples:
                df = self._parse_tsv(Path(f), sample) if f.endswith('.tsv') else self._parse_vcf(Path(f), sample)
                if not df.empty: all_data.append(df); parsed_samples.add(sample)
        if not all_data: logging.error("❌ No valid variants extracted!"); return False
        self.master_df = pd.concat(all_data, ignore_index=True)
        self.master_df['Variant_Signature'] = self.master_df['POS'].astype(str) + "_" + self.master_df['REF'] + ">" + self.master_df['ALT']
        self.total_samples = self.master_df['Sample_ID'].nunique()
        self.master_df.to_csv(self.output_dir / "Table1_Unified_HighQuality_Variants.csv", index=False)
        logging.info(f"✅ {self.total_samples} samples, {len(self.master_df)} variants")
        if self.genome_length == 0: self.genome_length = int(self.master_df['POS'].max() * 1.02)
        return True

    # ====== Population Genetics (Pool-Seq π) ======
    def _compute_tajimas_d(self, n, S, pi_window):
        if S < 3 or n < 2: return np.nan
        a1 = sum(1.0/i for i in range(1, n)); a2 = sum(1.0/(i**2) for i in range(1, n))
        b1 = (n + 1.0) / (3.0 * (n - 1.0)); b2 = (2.0 * (n**2 + n + 3.0)) / (9.0 * n * (n - 1.0))
        c1 = b1 - 1.0 / a1; c2 = b2 - (n + 2.0) / (a1 * n) + (a2 / (a1**2))
        e1 = c1 / a1; e2 = c2 / (a1**2 + a2)
        V = e1 * S + e2 * S * (S - 1.0)
        if V <= 1e-10: return np.nan
        return (pi_window - (S / a1)) / math.sqrt(V)

    def compute_popgen_sliding_window(self, win_size=0, step=0):
        """Pool-Seq expected heterozygosity: π per site from weighted VAF"""
        n = self.total_samples
        if n < 2: return
        max_pos = int(self.master_df['POS'].max())
        if win_size <= 0 or step <= 0:
            step = max(50, max_pos // 50) if max_pos > 500 else 25
            win_size = step * 2

        site_stats = []
        for pos, group in self.master_df.groupby('POS'):
            alt_alleles = group['ALT'].unique()
            p_alts = {alt: group[group['ALT'] == alt]['ALT_FREQ'].sum() / n for alt in alt_alleles}
            p_ref = max(0.0, 1.0 - sum(p_alts.values()))
            sum_sq = p_ref**2 + sum(p**2 for p in p_alts.values())
            pi_site = (1.0 - sum_sq) * (n / (n - 1.0)) if n > 1 else 0
            site_stats.append({'POS': pos, 'Pi': pi_site})

        site_df = pd.DataFrame(site_stats).set_index('POS')
        results = []
        for start in range(1, max_pos, step):
            end = start + win_size - 1
            in_window = site_df[(site_df.index >= start) & (site_df.index <= end)]
            S = len(in_window); pi_window = in_window['Pi'].sum()
            results.append({'BIN_START': start, 'BIN_END': end, 'BIN_MID': (start + end) / 2,
                'N_VARIANTS_S': S, 'PI': pi_window, 'TAJIMA_D': self._compute_tajimas_d(n, S, pi_window)})
        self.popgen_df = pd.DataFrame(results)
        self.popgen_df.to_csv(self.output_dir / "Table2_Sliding_Window_PopGen_Stats.csv", index=False)

    # ====== Plotting ======
    def _draw_gene_blocks(self, ax):
        ax.set_xlim(0, self.genome_length); ax.set_ylim(0, 1)
        ax.get_yaxis().set_visible(False)
        for sp in ['top', 'left', 'right']: ax.spines[sp].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.tick_params(axis='x', bottom=True, labelbottom=True)
        ax.hlines(0.5, 0, self.genome_length, color='black', linewidth=1.5, zorder=1)
        colors = ['#FF9F1C', '#2ECC71', '#3498DB', '#9B59B6', '#E74C3C', '#F1C40F']
        for i, feat in enumerate(self.gene_features):
            s, e, lbl = feat['start'], feat['end'], feat['name']
            ax.add_patch(patches.Rectangle((s, 0.2), e-s, 0.6, linewidth=1, edgecolor='black',
                facecolor=colors[i % len(colors)], zorder=2))
            ax.text((s+e)/2, 0.5, lbl, ha='center', va='center', fontsize=10, fontweight='bold', color='white', zorder=3)

    def plot_molecular_landscape(self):
        logging.info("🎨[1/8] Figure 1")
        colors = {'Transition (Ts)': '#d62728', 'Transversion (Tv)': '#1f77b4', 'Indel': '#ff7f0e', 'Complex/Unknown': 'gray'}
        plt.figure(figsize=(16, 6))
        pos_stats = self.master_df.groupby(['POS', 'Molecular_Type']).agg(Mean_AF=('ALT_FREQ', 'mean'), Sample_Count=('Sample_ID', 'nunique')).reset_index()
        pos_stats.to_csv(self.output_dir / "Figure1A_All_Variants_Landscape_Data.csv", index=False)
        sns.scatterplot(data=pos_stats, x='POS', y='Mean_AF', hue='Molecular_Type', size='Sample_Count', sizes=(20, 400), palette=colors, alpha=0.75, edgecolor="k", legend=False)
        plt.title(f'Figure 1A: Genomic Mutation Landscape\n({self.virus_name})', fontsize=18, fontweight='bold', pad=20)
        plt.xlabel('Genomic Position (bp)', fontsize=14); plt.ylabel('Mean Allele Frequency', fontsize=14)
        plt.tight_layout(); plt.savefig(self.output_dir / "Figure1A_All_Variants_Landscape.png", dpi=300, bbox_inches='tight')
        plt.savefig(self.output_dir / "Figure1A_All_Variants_Landscape.pdf", dpi=300, bbox_inches='tight'); plt.close()

        top_50 = self.master_df['Variant_Signature'].value_counts().head(50).index
        top_df = self.master_df[self.master_df['Variant_Signature'].isin(top_50)]
        plt.figure(figsize=(16, 6))
        ps50 = top_df.groupby(['POS', 'Molecular_Type']).agg(Mean_AF=('ALT_FREQ', 'mean'), Sample_Count=('Sample_ID', 'nunique')).reset_index()
        ps50.to_csv(self.output_dir / "Figure1B_Top50_Variants_Landscape_Data.csv", index=False)
        sns.scatterplot(data=ps50, x='POS', y='Mean_AF', hue='Molecular_Type', size='Sample_Count', sizes=(20, 400), palette=colors, alpha=0.75, edgecolor="k")
        plt.title(f'Figure 1B: Top 50 Shared Variants\n({self.virus_name})', fontsize=18, fontweight='bold', pad=20)
        plt.xlabel('Genomic Position (bp)', fontsize=14); plt.ylabel('Mean Allele Frequency', fontsize=14)
        plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left'); plt.tight_layout()
        plt.savefig(self.output_dir / "Figure1B_Top50_Variants_Landscape.png", dpi=300, bbox_inches='tight')
        plt.savefig(self.output_dir / "Figure1B_Top50_Variants_Landscape.pdf", dpi=300, bbox_inches='tight'); plt.close()

    def plot_tstv_pie(self):
        logging.info("🎨[2/8] Figure 2")
        plt.figure(figsize=(8, 8))
        counts = self.master_df['Molecular_Type'].value_counts()
        pd.DataFrame({'Molecular_Type': counts.index, 'Count': counts.values}).to_csv(self.output_dir / "Figure2_TsTv_Pie_Data.csv", index=False)
        ts, tv = counts.get('Transition (Ts)', 0), counts.get('Transversion (Tv)', 0)
        ratio = round(ts/tv, 2) if tv > 0 else "High"
        plt.pie(counts.values, labels=counts.index, autopct='%1.1f%%', startangle=140, colors=sns.color_palette("Set2"), wedgeprops={'edgecolor': 'black'})
        plt.title(f'Figure 2: Molecular Mutation Signature (Ts/Tv={ratio})', fontsize=16, fontweight='bold', pad=20)
        plt.tight_layout(); plt.savefig(self.output_dir / "Figure2_TsTv_Pie.png", dpi=300, bbox_inches='tight')
        plt.savefig(self.output_dir / "Figure2_TsTv_Pie.pdf", dpi=300, bbox_inches='tight'); plt.close()

    def plot_functional_pie(self):
        valid = self.master_df[self.master_df['Functional_Type'] != 'Unannotated']
        if valid.empty: return
        logging.info("🎨[3/8] Figure 3")
        plt.figure(figsize=(8, 8))
        counts = valid['Functional_Type'].value_counts()
        pd.DataFrame({'Functional_Type': counts.index, 'Count': counts.values}).to_csv(self.output_dir / "Figure3_Functional_Pie_Data.csv", index=False)
        plt.pie(counts.values, labels=counts.index, autopct='%1.1f%%', startangle=140, colors=['#ff9999', '#66b3ff'], wedgeprops={'edgecolor': 'black'})
        plt.title(f'Figure 3: Protein Functional Impact\n({self.virus_name})', fontsize=16, fontweight='bold', pad=20)
        plt.tight_layout(); plt.savefig(self.output_dir / "Figure3_Functional_Pie.png", dpi=300, bbox_inches='tight')
        plt.savefig(self.output_dir / "Figure3_Functional_Pie.pdf", dpi=300, bbox_inches='tight'); plt.close()

    def plot_heatmap(self):
        logging.info("🎨[4/8] Figure 4")
        top_vars = self.master_df['Variant_Signature'].value_counts().head(50).index
        pivot_df = self.master_df[self.master_df['Variant_Signature'].isin(top_vars)].pivot_table(
            index='Sample_ID', columns='Variant_Signature', values='ALT_FREQ', fill_value=0)
        if pivot_df.shape[0] < 2 or pivot_df.shape[1] < 2: return
        pivot_df = pivot_df.loc[:, pivot_df.var() > 0]
        if pivot_df.shape[1] < 2: return
        pivot_df.to_csv(self.output_dir / "Figure4_Clustermap_Data.csv")
        plt.figure(figsize=(14, max(8, len(pivot_df)*0.3)))
        cg = sns.clustermap(pivot_df, cmap="YlOrRd", metric="euclidean", method="ward", linewidths=0.5)
        cg.ax_col_dendrogram.set_title(f'Figure 4: Hierarchical Clustering of Top 50 Shared Variants\n({self.virus_name})', fontsize=16, fontweight='bold', pad=20)
        plt.savefig(self.output_dir / "Figure4_Clustermap.png", dpi=300, bbox_inches='tight')
        plt.savefig(self.output_dir / "Figure4_Clustermap.pdf", dpi=300, bbox_inches='tight'); plt.close()

    def plot_allele_frequency_spectrum(self):
        logging.info("🎨[5/8] Figure 5")
        self.master_df[['Sample_ID', 'Variant_Signature', 'ALT_FREQ']].to_csv(self.output_dir / "Figure5_AFS_Data.csv", index=False)
        plt.figure(figsize=(10, 6))
        sns.histplot(data=self.master_df, x='ALT_FREQ', bins=50, kde=True, color='purple')
        plt.title(f'Figure 5: Allele Frequency Spectrum\n({self.virus_name})', fontsize=16, fontweight='bold', pad=20)
        plt.xlabel('Intra-host AF', fontsize=14); plt.ylabel('Count', fontsize=14)
        plt.tight_layout(); plt.savefig(self.output_dir / "Figure5_AFS.png", dpi=300, bbox_inches='tight')
        plt.savefig(self.output_dir / "Figure5_AFS.pdf", dpi=300, bbox_inches='tight'); plt.close()

    def plot_variant_density(self):
        logging.info("🎨[6/8] Figure 6")
        self.master_df[['Sample_ID', 'Variant_Signature', 'POS']].to_csv(self.output_dir / "Figure6_Variant_Density_Data.csv", index=False)
        if self.gene_features:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 6.5), gridspec_kw={'height_ratios': [5, 1]}, sharex=True)
            sns.kdeplot(data=self.master_df, x='POS', bw_adjust=0.2, fill=True, color='teal', alpha=0.5, ax=ax1)
            ax1.set_title(f'Figure 6: Variant Density + Annotation\n({self.virus_name})', fontsize=16, fontweight='bold', pad=20)
            ax1.set_ylabel('Density', fontsize=14); sns.rugplot(data=self.master_df, x='POS', color='black', alpha=0.1, height=0.05, ax=ax1)
            self._draw_gene_blocks(ax2); ax2.set_xlabel('Genomic Position (bp)', fontsize=14, fontweight='bold')
        else:
            fig, ax1 = plt.subplots(figsize=(14, 5))
            sns.kdeplot(data=self.master_df, x='POS', bw_adjust=0.2, fill=True, color='teal', alpha=0.5, ax=ax1)
            ax1.set_title(f'Figure 6: Variant Density\n({self.virus_name})', fontsize=16, fontweight='bold', pad=20)
            ax1.set_xlabel('Genomic Position (bp)', fontsize=14); ax1.set_ylabel('Density', fontsize=14)
            sns.rugplot(data=self.master_df, x='POS', color='black', alpha=0.1, height=0.05, ax=ax1)
        plt.tight_layout(); plt.savefig(self.output_dir / "Figure6_Variant_Density.png", dpi=300, bbox_inches='tight')
        plt.savefig(self.output_dir / "Figure6_Variant_Density.pdf", dpi=300, bbox_inches='tight'); plt.close()

    def plot_af_violin(self):
        logging.info("🎨[7/8] Figure 7")
        self.master_df[['Sample_ID', 'Variant_Signature', 'Molecular_Type', 'ALT_FREQ']].to_csv(self.output_dir / "Figure7_AF_Violin_Data.csv", index=False)
        plt.figure(figsize=(10, 6))
        sns.violinplot(data=self.master_df, x='Molecular_Type', y='ALT_FREQ', hue='Molecular_Type', palette='Set2', inner='quartile', legend=False)
        sns.stripplot(data=self.master_df, x='Molecular_Type', y='ALT_FREQ', color='black', alpha=0.3, size=3, jitter=True)
        plt.title(f'Figure 7: AF Distribution by Mutation Type\n({self.virus_name})', fontsize=16, fontweight='bold', pad=20)
        plt.xlabel('Mutation Type', fontsize=14); plt.ylabel('Intra-host AF', fontsize=14)
        plt.tight_layout(); plt.savefig(self.output_dir / "Figure7_AF_Violin.png", dpi=300, bbox_inches='tight')
        plt.savefig(self.output_dir / "Figure7_AF_Violin.pdf", dpi=300, bbox_inches='tight'); plt.close()

    def plot_popgen_dynamics(self):
        logging.info("🎨[8/8] Figure 8")
        if self.popgen_df is None or self.popgen_df.empty: return
        plot_df = self.popgen_df.dropna(subset=['TAJIMA_D']).copy()
        if plot_df.empty: return
        plot_df.to_csv(self.output_dir / "Figure8_PopGen_Dynamics_Data.csv", index=False)
        if self.gene_features:
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(16, 10.5), sharex=True, gridspec_kw={'height_ratios': [4, 4, 1]})
        else:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 9), sharex=True)
        fig.suptitle(f"Figure 8: Evolutionary Dynamics\n({self.virus_name})", fontsize=18, fontweight='bold', y=0.96)
        sns.lineplot(data=plot_df, x='BIN_MID', y='PI', ax=ax1, color='#e74c3c', linewidth=2.5)
        ax1.fill_between(plot_df['BIN_MID'], plot_df['PI'], color='#e74c3c', alpha=0.2)
        ax1.set_ylabel(r'$\pi$', fontsize=14, fontweight='bold'); ax1.set_title("A: Diversity Profile", loc='left', fontsize=14)
        ax1.grid(True, linestyle='--', alpha=0.6)
        sns.lineplot(data=plot_df, x='BIN_MID', y='TAJIMA_D', ax=ax2, color='#2980b9', linewidth=2.5)
        ax2.fill_between(plot_df['BIN_MID'], plot_df['TAJIMA_D'], color='#2980b9', alpha=0.2)
        ax2.axhline(0, color='black', linestyle='--', linewidth=1.5, label='Neutrality (D=0)')
        tmin = plot_df['TAJIMA_D'].min()
        if not pd.isna(tmin) and tmin < -1.5:
            ax2.axhspan(tmin, -1.5, color='salmon', alpha=0.1, label='Directional Selection')
        ax2.set_ylabel("Tajima's D", fontsize=14, fontweight='bold'); ax2.set_title("B: Selection Inference", loc='left', fontsize=14)
        ax2.legend(loc='upper right'); ax2.grid(True, linestyle='--', alpha=0.6)
        if self.gene_features: self._draw_gene_blocks(ax3); ax3.set_xlabel('Position (bp)', fontsize=14, fontweight='bold')
        else: ax2.set_xlabel('Position (bp)', fontsize=14, fontweight='bold')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.savefig(self.output_dir / "Figure8_PopGen_Dynamics.png", dpi=300, bbox_inches='tight')
        plt.savefig(self.output_dir / "Figure8_PopGen_Dynamics.pdf", dpi=300, bbox_inches='tight'); plt.close()

    def plot_pca_lineages(self):
        if not HAS_SKLEARN: return
        logging.info("🌟[BONUS] Figure 9: PCA")
        var_counts = self.master_df['Variant_Signature'].value_counts()
        shared_vars = var_counts[var_counts > 1].index
        if len(shared_vars) < 5 or self.master_df['Sample_ID'].nunique() < 6:
            logging.info("  -> [跳过] 样本<6 或突变太少"); return
        pivot_df = self.master_df[self.master_df['Variant_Signature'].isin(shared_vars)].pivot_table(
            index='Sample_ID', columns='Variant_Signature', values='ALT_FREQ', fill_value=0)
        pivot_df = pivot_df.loc[:, pivot_df.var() > 0]
        if pivot_df.shape[1] < 2: logging.info("  -> [跳过] 变异特征方差为0"); return
        X_scaled = StandardScaler().fit_transform(pivot_df)
        n_pcs_full = min(X_scaled.shape[0]-1, X_scaled.shape[1], 10)
        pca_full = PCA(n_components=n_pcs_full); X_pca_full = pca_full.fit_transform(X_scaled)
        cumvar = np.cumsum(pca_full.explained_variance_ratio_)
        n_pcs_use = max(3, min(int(np.searchsorted(cumvar, 0.80)+1), n_pcs_full))
        X_pca_use = X_pca_full[:, :n_pcs_use]
        max_k = max(2, min(7, len(pivot_df)//3))
        best_k, best_score = 2, -1.0
        for k in range(2, max_k+1):
            labels_tmp = KMeans(n_clusters=k, random_state=42, n_init=10).fit_predict(X_pca_use)
            score = silhouette_score(X_pca_use, labels_tmp)
            if score > best_score: best_score, best_k = score, k
        logging.info(f"  -> K={best_k} (Silhouette={best_score:.3f})")
        kmeans = KMeans(n_clusters=best_k, random_state=42, n_init=10); clusters = kmeans.fit_predict(X_pca_use)
        X_pca_2d = X_pca_full[:, :2]; var_r = pca_full.explained_variance_ratio_
        pca_df = pd.DataFrame({'PC1': X_pca_2d[:, 0], 'PC2': X_pca_2d[:, 1],
            'Cluster': [f"Lineage {c+1}" for c in clusters], 'Sample_ID': pivot_df.index})
        pca_df.to_csv(self.output_dir / "Figure9_PCA_Lineages_Data.csv", index=False)
        plt.figure(figsize=(10, 8))
        sns.kdeplot(data=pca_df, x="PC1", y="PC2", hue="Cluster", fill=True, alpha=0.2, palette="Dark2", legend=False)
        sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Cluster", palette="Dark2", s=100, edgecolor='black', alpha=0.9)
        plt.xlabel(f'PC1 ({var_r[0]*100:.1f}%)', fontsize=14, fontweight='bold')
        plt.ylabel(f'PC2 ({var_r[1]*100:.1f}%)', fontsize=14, fontweight='bold')
        plt.title(f'Figure 9: PCA of Viral Populations (Auto-K={best_k})\n({self.virus_name})', fontsize=16, fontweight='bold', pad=20)
        plt.legend(title='Sub-lineages', bbox_to_anchor=(1.05, 1), loc='upper left'); plt.tight_layout()
        plt.savefig(self.output_dir / "Figure9_PCA_Lineages.png", dpi=300, bbox_inches='tight')
        plt.savefig(self.output_dir / "Figure9_PCA_Lineages.pdf", dpi=300, bbox_inches='tight'); plt.close()
        if n_pcs_use >= 3:
            X_pca_3d = X_pca_full[:, :3]
            pca_3d = pd.DataFrame({'PC1': X_pca_3d[:, 0], 'PC2': X_pca_3d[:, 1], 'PC3': X_pca_3d[:, 2],
                'Cluster': [f"Lineage {c+1}" for c in clusters], 'Sample_ID': pivot_df.index})
            palette = sns.color_palette("Dark2", n_colors=best_k)
            color_dict = {f"Lineage {i+1}": c for i, c in enumerate(palette)}
            fig = plt.figure(figsize=(10, 8)); ax3d = fig.add_subplot(111, projection='3d')
            for label in sorted(pca_3d['Cluster'].unique()):
                s = pca_3d[pca_3d['Cluster'] == label]
                ax3d.scatter(s['PC1'], s['PC2'], s['PC3'], label=label, c=[color_dict[label]], s=80, edgecolors='w', alpha=0.9)
            ax3d.set_title(f'Figure 9 (3D): PCA (Auto-K={best_k})\n({self.virus_name})', fontsize=15, fontweight='bold')
            ax3d.set_xlabel(f"PC1 ({var_r[0]*100:.1f}%)", fontweight='bold')
            ax3d.set_ylabel(f"PC2 ({var_r[1]*100:.1f}%)", fontweight='bold')
            ax3d.set_zlabel(f"PC3 ({var_r[2]*100:.1f}%)", fontweight='bold')
            ax3d.view_init(elev=20, azim=45); ax3d.grid(True, linestyle=':', alpha=0.6)
            ax3d.legend(title='Sub-lineages', bbox_to_anchor=(1.15, 0.8), loc='upper left')
            plt.tight_layout(); plt.savefig(self.output_dir / "Figure9_PCA_Lineages_3D.png", dpi=300, bbox_inches='tight')
            plt.savefig(self.output_dir / "Figure9_PCA_Lineages_3D.pdf", dpi=300, bbox_inches='tight'); plt.close()


def main():
    parser = argparse.ArgumentParser(description="Unified Viral Variant Analyzer")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-d", "--min-depth", type=int, default=50)
    parser.add_argument("-f", "--min-af", type=float, default=0.05)
    parser.add_argument("-v", "--virus-name", type=str, default="")
    parser.add_argument("-a", "--acc", type=str, default="")
    parser.add_argument("--win-size", type=int, default=0)
    parser.add_argument("--win-step", type=int, default=0)
    args = parser.parse_args()
    try:
        a = UnifiedVariantAnalyzer(input_dir=args.input, output_dir=args.output,
            min_depth=args.min_depth, min_af=args.min_af, virus_name=args.virus_name, acc_id=args.acc)
        if a.load_and_merge():
            a.plot_molecular_landscape(); a.plot_tstv_pie(); a.plot_functional_pie()
            a.plot_heatmap(); a.plot_allele_frequency_spectrum(); a.plot_variant_density()
            a.plot_af_violin(); a.compute_popgen_sliding_window(win_size=args.win_size, step=args.win_step)
            a.plot_popgen_dynamics(); a.plot_pca_lineages()
            logging.info(f"🏆 Done -> {args.output}")
    except Exception as e:
        import traceback; traceback.print_exc()
        logging.error(f"Critical error: {e}")

if __name__ == "__main__":
    main()
