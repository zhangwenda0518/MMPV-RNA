# Known Virus Analysis Pipeline

10-stage automated pipeline for known virus detection, variant analysis, full-length assembly, and post-hoc characterization.

## Quick Start

```bash
python auto_known_virus.py --stage all --filter \
  --reads_dir clean_reads/ \
  --output_dir ./virus_analysis/ \
  --ref_info ref_info.tsv \
  --reference ref.fasta \
  --genes_cov virus_genes_cov.tsv \
  --tool salmon --threads 60 --align_threads 8 \
  --variant_caller freebayes --snpeff --snpgenie \
  --resume
```

## Pipeline Flow

```
Input: clean_reads/ + ref_info.tsv + ref.fasta

[1/10] batch_virus_depth40.py       → 1_FastViromeExplorer/
[2/10] filter_summary.py             → high_conf.summary.tsv
[3/10] batch_virus_variants.py       → 2_Virus_variants_Results/
[4/10] batch_virus_full.py           → 3_Virus_assemblies_final/
       └── virus-full.py (engine)
[5/10] extract_full_fasta.py         → 4_assemblies_clean/
[6/10] batch_plot_virus_depth.py     → 5_post_analysis/
       ├── virus_variants_analyzer   VCF viz
       ├── virus_vcf_pipeline        merge + PCA
       ├── snpeff_analysis           SnpEff macro
       ├── snpeff2maf + viral_maftools  MAF waterfall
       └── snpgenie_master           dN/dS + 3D PCA
[7/10] capheine_pipeline.py          → 6_capheine/
       ├── gbk_extractor             CDS extraction
       └── visual_codon_miner        codon viz
[8/10] virus_auto_pipeline3.py       → 7_similarity/
[9/10] batch_virema_dvg.py           → 8_virema_dvg/
       └── virema_summary_report.R   Circos 4-track
[10/10] generate_pipeline_report.py  → Pipeline_Summary_Report.html
```

## Stage Control

```bash
# Single stage
--stage detect|filter|variants|full|extract|post|capheine|similarity|dvg|report

# Full pipeline with filtering
--stage all --filter

# Resume from checkpoint
--resume

# Preview only
--dry_run
```

## Script Index

| Script | Stage | Lines | Function |
|--------|:----:|:----:|----------|
| `auto_known_virus.py` | — | 755 | 10-stage orchestrator |
| `batch_virus_depth40.py` | 1 | 832 | Salmon pseudo-alignment + Poisson filtering |
| `filter_summary.py` | 2 | 236 | Coverage/depth/reads threshold filtering |
| `batch_virus_variants.py` | 3 | 1229 | Bowtie2 + FreeBayes + SnpEff + SnpGenie |
| `batch_virus_full.py` | 4 | 240 | Assembly task scheduler |
| `virus-full.py` | 4 | 1120 | 12-step de novo assembly engine |
| `extract_full_fasta.py` | 5 | 123 | Extract longest contigs |
| `batch_plot_virus_depth.py` | 6 | 784 | 3-in-1 viz (depth/freq/meta) |
| `virus_variants_analyzer.py` | 6 | 477 | VCF landscape + PCA |
| `virus_vcf_pipeline.py` | 6 | 186 | bcftools merge + VCF2Dis + VCF2PCA |
| `snpeff_analysis.py` | 6 | 310 | SnpEff macro stats + OncoPrint |
| `snpeff2maf.py` | 6 | 333 | VCF to MAF conversion |
| `viral_maftools.R` | 6 | 233 | maftools waterfall + lollipop |
| `snpgenie_master.py` | 6 | 671 | dN/dS + Auto-K 3D PCA + Kruskal-Wallis |
| `capheine_pipeline.py` | 7 | 471 | HyPhy positive selection (FEL/MEME/BUSTED/PRIME) |
| `gbk_extractor.py` | 7 | 122 | GenBank CDS extraction |
| `visual_codon_miner.py` | 7 | 234 | Positive selection site visualization |
| `virus_auto_pipeline3.py` | 8 | 1265 | Pairwise SDT similarity panorama |
| `batch_virema_dvg.py` | 9 | 544 | ViReMa DVG detection |
| `virema_summary_report.R` | 9 | 335 | Circos 4-track recombination plot |
| `generate_pipeline_report.py` | 10 | 374 | HTML summary + AI prompts |
| `consensus_extract.py` | util | 231 | Consensus sequence QA + N-filling |

## Dependencies

**Python packages**: polars, pysam, biopython, pandas, matplotlib, seaborn, scipy, scikit-learn, tqdm, colorlog

**External tools**: bowtie2, samtools, freebayes, bcftools, snpEff, salmon/kallisto, megahit/spades, mafft, iqtree, hyphy, cawlign, drhip, multiqc, pandepth, ViReMa, BBMerge

**R packages**: circlize, ComplexHeatmap, maftools, ggplot2, dplyr, tidyr, viridis

## Archives

- `archive/` — 52 historical versions of batch_virus_depth*.py
- `utils/` — 6 standalone tools (snpeff_build, consensus_to_proteins, etc.)
