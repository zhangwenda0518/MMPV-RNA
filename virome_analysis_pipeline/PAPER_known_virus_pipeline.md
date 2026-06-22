# Known Virus Analysis Pipeline (KVAP): A 10-Stage Integrated Framework for Comprehensive Plant Virus Detection, Variant Analysis, Assembly, and Evolutionary Characterization

## Abstract

Comprehensive characterization of known plant viruses from high-throughput sequencing data requires an integrated analytical framework that spans detection, variant calling, full-length genome assembly, functional annotation, and evolutionary analysis. We present KVAP (Known Virus Analysis Pipeline), a 10-stage automated bioinformatics pipeline that unifies these tasks within a single orchestrated workflow. KVAP integrates seven established bioinformatics methods—Salmon pseudo-alignment, FreeBayes variant calling, SnpEff functional annotation, SNPGenie population genetics, multi-tool de novo assembly with Divine Fusion refinement, HyPhy positive selection analysis, and ViReMa defective genome detection—under a modular, profile-driven architecture. The pipeline was validated on 7 Lycium transcriptome samples, identifying 9 virus species, calling 25 high-confidence variant records, assembling 22 complete viral genomes (20 achieving perfect quality), and performing comprehensive post-hoc evolutionary characterization. KVAP is implemented in Python with a YAML-configured orchestrator, supports parallel execution across stages, and generates an interactive HTML summary report with embedded charts. The pipeline is available at [repository URL] under an open-source license.

**Keywords**: virus detection, variant calling, genome assembly, positive selection, bioinformatics pipeline, plant virology, metatranscriptomics

---

## 1. Introduction

The increasing availability of public transcriptome datasets has created unprecedented opportunities for virus discovery and characterization in plant systems. However, extracting meaningful biological insights from these data requires navigating a fragmented landscape of specialized bioinformatics tools. A typical known-virus analysis workflow involves at least seven distinct computational tasks: (i) read quantification against reference genomes, (ii) quality filtering to remove false-positive detections, (iii) single-nucleotide variant (SNV) calling and annotation, (iv) de novo genome assembly for complete sequence recovery, (v) population genetic analysis to detect selection signals, (vi) codon-level positive selection testing, and (vii) defective genome or recombination detection.

Each of these tasks typically requires a separate tool with its own input/output conventions, parameter spaces, and dependency chains. Researchers must manually chain these tools together, manage intermediate files, and ensure consistent parameter settings across steps—a process that is error-prone, time-consuming, and difficult to reproduce. While workflow managers such as Nextflow and Snakemake can partially address these challenges, they do not provide opinionated, domain-specific parameter defaults or integrated quality-control checks tuned for the characteristics of viral genomes (e.g., small genome size, high mutation rates, haploid genotypes, and the presence of defective interfering particles).

Several existing pipelines address portions of this workflow. FastViromeExplorer provides rapid virus detection via pseudo-alignment (Tithi et al., 2018). DRHIP and the HyPhy package offer comprehensive positive selection analysis (Kosakovsky Pond et al., 2020). VIGA assembles viral genomes from metagenomic data (González-Tortuero et al., 2019). However, no single framework integrates detection, variant analysis, assembly, population genetics, positive selection, and defective genome analysis into a cohesive pipeline with automated data flow, consistent quality filtering, and publication-ready visualization output.

Here we present KVAP (Known Virus Analysis Pipeline), a 10-stage automated framework that addresses this gap. KVAP is designed around five design principles: (1) **modularity** — each stage can run independently; (2) **profile-driven configuration** — a single YAML file captures all parameters for reproducibility; (3) **defensive quality control** — Poisson ratio filtering, dual-track (genome + gene) validation, and dynamic depth thresholds reduce false-positive variant calls; (4) **parallelism** — stages and intra-stage tasks exploit multi-core architectures; and (5) **comprehensive output** — an interactive HTML report with embedded charts and structured AI interpretation prompts.

---

## 2. Methods

### 2.1 Pipeline Architecture

KVAP is implemented as a collection of 20 Python scripts and 2 R scripts orchestrated by `auto_known_virus.py`. The orchestrator defines a 10-stage workflow with explicit data dependencies between stages. Each stage invokes one or more specialized analysis scripts through subprocess calls with validated parameter passing. A YAML profile system (`default_profile.yaml`) captures all configurable parameters, enabling reproducible execution with a single command.

The pipeline architecture is organized into four logical layers:

1. **Detection & Filtering** (Stages 1–2): Virus identification and quality filtering
2. **Variant & Assembly** (Stages 3–5): Variant calling, functional annotation, and genome reconstruction
3. **Evolutionary Analysis** (Stages 6–9): Population genetics, positive selection, similarity analysis, and recombination detection
4. **Reporting** (Stage 10): Interactive HTML report generation

```
                    ┌──────────────────────────────────────────┐
                    │           Input: Clean Reads              │
                    │     + Reference FASTA + Ref Info TSV      │
                    └──────────────────┬───────────────────────┘
                                       │
        ┌──────────────────────────────┼──────────────────────────────┐
        │ Stage 1                      │                      Stage 9 │
        │ batch_virus_depth.py       │              batch_virema_dvg│
        │ Salmon pseudo-alignment      │              ViReMa + Circos │
        │ → best_summary.tsv           │              → 8_virema_dvg/ │
        └──────────────┬───────────────┘                              │
                       │                                              │
        ┌──────────────▼───────────────┐                              │
        │ Stage 2: filter_summary.py   │                              │
        │ Cov > 50%, Depth > 5, Reads  │                              │
        │ > 100 → high_conf.tsv        │                              │
        └──────────────┬───────────────┘                              │
                       │                                              │
        ┌──────────────▼───────────────────────────────────────┐      │
        │ Stage 3: batch_virus_variants.py                     │      │
        │ bowtie2 → FreeBayes → SnpEff → SNPGenie → Consensus  │      │
        │ → 2_Virus_variants_Results/                          │      │
        └───┬───────────┬──────────┬──────────┬────────────────┘      │
            │           │          │          │                       │
    ┌───────▼──┐ ┌──────▼───┐ ┌────▼────┐ ┌──▼──────────┐            │
    │ Stage 4  │ │ Stage 6  │ │Stage 7 │ │ Stage 8      │            │
    │ Assembly │ │Post-hoc  │ │Capheine│ │ Similarity   │            │
    │ → 3_asm/ │ │→ 5_post/ │ │→ 6_cap/│ │ → 7_sim/     │            │
    └──┬───────┘ └──────────┘ └────────┘ └──────────────┘            │
       │                                                              │
  ┌────▼──────┐                                                       │
  │ Stage 5   │                                                       │
  │ Extract   │                                                       │
  │→ 4_clean/ │                                                       │
  └───────────┘                                                       │
       │                                                              │
       └──────────────────────┬───────────────────────────────────────┘
                              │
                    ┌─────────▼──────────┐
                    │ Stage 10: Report   │
                    │ Interactive HTML   │
                    └────────────────────┘
```

### 2.2 Stage 1: Virus Detection and Quantification

**Script**: `batch_virus_depth.py`

Virus detection employs Salmon (v1.10) pseudo-alignment (Patro et al., 2017) against a curated plant virus reference database. For each sample, Salmon quantifies read abundance per reference sequence using an expectation-maximization algorithm that probabilistically assigns multi-mapping reads. The following normalized abundance metrics are computed:

- **CPM** (Counts Per Million): Reads per reference normalized by total mapped reads × 10⁶
- **RPM / FPKM / TPM**: Standard transcript-level normalization metrics
- **Rel_Abund(%)**: Relative abundance as percentage of total viral reads

To distinguish genuine viral signal from spurious mapping artifacts, we implement a Poisson Ratio filter. For each virus, the observed coverage distribution across the reference genome is compared against the expected coverage under a Poisson model given the total read count and average read length. A Poisson Ratio ≥ 0.3 indicates that reads are distributed across the genome rather than concentrated at a few positions—a hallmark of true viral presence rather than non-specific mapping.

Additionally, a dual-track (genome + gene) validation strategy is employed when gene coverage data (`genes_cov.tsv`) is available. Viruses must pass either the genome-wide coverage threshold (Track A: coverage ≥ 10%, Poisson Ratio ≥ 0.3) or the gene-level coverage threshold (Track B: gene total coverage ≥ 80%, gene average coverage ≥ 5).

### 2.3 Stage 2: High-Confidence Filtering

**Script**: `filter_summary.py`

Raw detection results are filtered to retain only high-confidence viral records. The default criteria require: (i) Rep_Coverage(%) ≥ 50%, (ii) Rep_MeanDepth ≥ 5×, and (iii) Asm_EM_Reads ≥ 100. Additional optional filters include minimum TPM, Poisson Ratio, and keyword-based species selection. The filter outputs per-sample and per-virus summary statistics as tab-separated tables for downstream epidemiological analysis.

### 2.4 Stage 3: Variant Calling and Functional Annotation

**Script**: `batch_virus_variants.py`

For each high-confidence virus-sample pair, reads are re-aligned against the virus reference genome using Bowtie2 (v2.4.4) in local alignment mode (Langmead & Salzberg, 2012). Single-nucleotide variants (SNVs) are called using FreeBayes (v1.3.6) in haploid mode (Garrison & Marth, 2012), which is appropriate for viral genomes. Dynamic depth thresholds are applied: DP ≥ 10 for samples with mean depth < 50×, DP ≥ 20 for depth 50–1000×, and DP ≥ 100 for depth > 1000×. The minimum allele frequency threshold is 5%.

Functional annotation of variants is performed using SnpEff (v5.1) (Cingolani et al., 2012) with custom databases built from viral GenBank records. Variants are classified as HIGH (e.g., frameshift, stop-gained), MODERATE (missense), LOW (synonymous), or MODIFIER (non-coding) impact. Consensus sequences are generated using iVar (Grubaugh et al., 2019) with minimum base quality Q20 and minimum depth 5×.

Population genetic parameters—nucleotide diversity (π), πN/πS ratios, and dN/dS—are computed using SNPGenie (Nelson et al., 2015) with 50-bp sliding windows. These metrics are calculated at both the orthologous block level (comparing each sample consensus to the reference) and the intra-host level (quantifying within-sample quasispecies diversity).

### 2.5 Stage 4: De Novo Full-Length Assembly

**Scripts**: `batch_virus_full.py` → `virus-full.py`

Complete viral genomes are reconstructed through a 12-step assembly pipeline (OmniVirusAssembler, virus-full.py) that integrates three de novo assemblers (MEGAHIT, SPAdes/RNAviralSPAdes, and PenguIN) followed by iterative refinement:

1. **Steps 1–2**: Multi-tool de novo assembly with refineC split-merge for contig deduplication
2. **Step 3**: Shiver-like BLASTn-based orientation correction and non-homologous region trimming
3. **Steps 4–8**: Divine Fusion—progressive skeleton merging using MAFFT multiple sequence alignment with reference-guided gap resolution; PVGA read-based extension with gap-size evaluation (>100 bp gaps trigger split); rmDup deduplication
4. **Step 9**: Iterative read-level polishing via minimap2 alignment followed by viral_consensus error correction (3 rounds)
5. **Step 10**: Dual-engine gap filling using gmcloser structural backfill and abyss-sealer Bloom-filter-based reads filling
6. **Step 11**: Circularity detection via BLASTn self-alignment of terminal regions
7. **Step 12**: Coverage visualization showing assembly evolution across all 11 intermediate steps

A "fast-track" mechanism bypasses intermediate fusion steps when the Shiver-cleanup skeleton already achieves ≥98% reference length and N50 ≥95% reference length, directly entering iterative polishing.

### 2.6 Stage 5: Assembly Extraction

**Script**: `extract_full_fasta.py`

The longest contig from each Ultimate Circular Result is extracted into a clean FASTA file, removing redundant contigs that may arise from the assembly of small circular genomes (e.g., viroids).

### 2.7 Stage 6: Post-hoc Evolutionary Characterization

Five parallel analysis modules process the variant and annotation outputs:

**6.1 VCF Visualization** (`virus_variants_analyzer.py`): Produces 11 publication-ready figures including a whole-genome variant landscape, Ts/Tv ratio analysis, allele frequency spectrum, variant density heatmap, and PCA of inter-sample genetic distances.

**6.2 Population Structure** (`virus_vcf_pipeline.py`): Merges per-sample VCFs using bcftools merge, then computes genetic distance matrices via VCF2Dis and performs PCA via VCF2PCACluster. Outputs include 2D/3D PCA plots, a neighbor-joining tree, and a merged VCF file for external analysis.

**6.3 Functional Impact Profiling** (`snpeff_analysis.py`): Generates a mutation Manhattan plot across the genome, per-gene mutational burden with stacked HIGH/MODERATE/LOW categorization, intra-host quasispecies diversity metrics (iSNV vs consensus), and a ComplexHeatmap OncoPrint of shared mutations across samples. Six statistical matrices (CSV) provide detailed variant-level data.

**6.4 Protein-Structure-Level Analysis** (`snpeff2maf.py` + `viral_maftools.R`): Converts SnpEff-annotated VCFs to Mutation Annotation Format (MAF) and leverages the maftools R package (Mayakonda et al., 2018) for oncoplot waterfall visualization, protein-domain lollipop plots, somatic interaction networks, and cohort-level Ti/Tv summaries.

**6.5 Selection Pressure Quantification** (`snpgenie_master.py`): Produces 11 figures spanning multiple analytical levels. Key analyses include: (i) dN vs dS joint distribution plots for inter-host and intra-host comparisons, (ii) per-gene dN/dS violin plots with Wilcoxon signed-rank tests for deviation from neutrality, (iii) 10,000-iteration bootstrap confidence intervals for dN/dS, (iv) Kruskal-Wallis non-parametric tests for cross-gene dN/dS distribution differences, (v) mutational spectrum analysis to detect substitution biases (e.g., APOBEC/ADAR editing signatures), and (vi) unsupervised machine learning clustering via PCA on the evolutionary feature matrix [dN, dS, πN, πS] with automatic K-selection using the Silhouette Score to eliminate human bias in sub-lineage identification. A 3D PCA (PC1–PC3) reveals hidden substructure potentially masked in 2D projections. Optional external metadata (geographic origin, host species) can replace K-Means coloring for epidemiological interpretation.

### 2.8 Stage 7: Codon-Level Positive Selection Analysis

**Scripts**: `gbk_extractor.py` → `capheine_pipeline.py`

For viruses with annotated coding sequences in the virus-annotations/ directory, GenBank files are parsed to extract codon-aware CDS alignments. The capheine_pipeline.py implements a pure-Python reimplementation of the CAPHEINE Nextflow pipeline (Verdonk & Callan, 2024), executing:

1. **cawlign**: Codon-aware alignment of sample sequences to reference CDS
2. **IQ-TREE**: Maximum-likelihood phylogenetic tree inference with GTR+I+G model
3. **HyPhy analyses** (Kosakovsky Pond et al., 2020):
   - FEL (Fixed Effects Likelihood): Pervasive site-level selection
   - MEME (Mixed Effects Model of Evolution): Episodic diversifying selection
   - BUSTED: Gene-wide episodic selection
   - PRIME: Property-informed amino acid substitution analysis
   - CONTRASTFEL/RELAX: Branch-set specific selection intensity comparison
4. **DRHIP**: Aggregation of HyPhy JSON outputs into combined summary and site tables

For identified positive selection sites, `visual_codon_miner.py` extracts population-level codon usage from the codon-aware alignments and generates per-site amino acid frequency bar charts, with correct frequency denominators (total population size, not filtered subset) to avoid inflating substitution frequencies.

### 2.9 Stage 8: Sequence Similarity Panorama

**Script**: `virus_auto_pipeline.py`

Pairwise sequence similarity is computed using SDT (Species Demarcation Tool) across all consensus sequences for each virus. Results are visualized as a clustered heatmap with hierarchical dendrograms, enabling rapid identification of genetically identical isolates (indicative of clonal propagation) versus divergent lineages.

### 2.10 Stage 9: Defective Genome and Recombination Detection

**Script**: `batch_virema_dvg.py`

Defective viral genomes (DVGs) and recombination events are detected using ViReMa (v0.29) (Routh & Johnson, 2014) with seed length 25, micro-indel threshold 15 bp, and defuzz=0 for strict breakpoint resolution. Paired-end reads are first merged using BBMerge to increase anchor length. The analysis pipeline includes: (i) strict sanitization to remove ambiguous bases, (ii) independent sandbox execution for each sample-virus pair, (iii) aggregation of BED/BEDPE/TXT outputs, and (iv) automated R-based Circos 4-track visualization (`virema_summary_report.R`) showing recombination breakpoints, duplication/deletion coverage, gene annotations, and donor-acceptor junction networks.

### 2.11 Stage 10: Interactive Report Generation

**Script**: `generate_pipeline_report.py`

An interactive HTML report is generated with embedded base64-encoded charts (PNG/PDF), per-virus metrics cards (sample count, coverage, CPM, depth, Poisson ratio), and a left sidebar navigation panel. The report aggregates outputs from all completed stages and optionally includes structured AI interpretation prompts for automated scientific narrative generation via large language models (DeepSeek/OpenAI API).

### 2.12 Implementation

The pipeline is implemented in Python 3.10+ with a YAML-configured orchestrator. Key dependencies include: Polars, PySAM, Biopython, pandas, matplotlib, seaborn, scipy, and scikit-learn. External bioinformatics tools required: bowtie2, samtools, freebayes, bcftools, SnpEff, Salmon, MEGAHIT/SPAdes, MAFFT, IQ-TREE, HyPhy, cawlign, DRHIP, MultiQC, pandepth, ViReMa, and BBMerge. R packages: circlize, ComplexHeatmap, maftools, ggplot2, dplyr, tidyr, viridis.

The pipeline supports both full execution (`--stage all --filter --profile default_profile.yaml`) and individual stage execution (`--stage detect|variants|full|post|capheine|dvg|report`). Checkpoint-based resume (`--resume`) enables recovery from partial failures. A dry-run mode (`--dry_run`) previews the execution plan without running computations.

---

## 3. Results

### 3.1 Validation Dataset

The pipeline was validated using 7 Lycium (goji berry) leaf transcriptome samples (NCBI BioProject accession pending). Raw paired-end reads were pre-processed with fastp (quality filtering) and host reads were depleted by mapping to the Lycium reference genome with Bowtie2. The resulting clean reads (FASTA format, 7 samples, paired-end) were used as input.

Reference genomes comprised 6,035 plant virus sequences from a curated database including complete RefSeq genomes. Reference metadata included NCBI taxonomy, segment information, and gene coverage profiles for dual-track validation.

### 3.2 Detection and Filtering

Stage 1 (Salmon pseudo-alignment) completed in 15 seconds (60 threads), detecting 39 virus records across 9 unique species from the 7 samples. Stage 2 filtering (Cov ≥ 50%, Depth ≥ 5×, Reads ≥ 100) retained 25 high-confidence records representing 5 virus species: Cytorhabdovirus sp. 'lycii' (4 samples), Potato spindle tuber viroid (7 samples), Rubber viroid India/2009 (7 samples), Tomato chlorotic dwarf viroid (6 samples), and Citrus exocortis Yucatan viroid (1 sample).

### 3.3 Variant Analysis

Stage 3 (variant calling + annotation) completed in 3.4 minutes (60 threads, 4 parallel jobs), processing all 25 virus-sample pairs. For Cytorhabdovirus sp. 'lycii' (GenBank: OR489165.1, 14,812 bp), the pipeline called 106 variant positions in the GS-3 sample (mean depth 22.6×), including 51 missense, 48 synonymous, and 7 regulatory variants. The Ti/Tv ratio of 3.2 was within the expected range for authentic RNA virus variants. SnpEff annotation identified mutations across all 6 viral genes (N, P, P4, M, G, L), with the L (polymerase) gene carrying the highest absolute mutational burden (45 variants).

### 3.4 Genome Assembly

Stage 4 (de novo assembly) successfully reconstructed 22 complete viral genomes from 25 tasks (88% success rate), with 20 achieving perfect quality (Fully-assembled.ok, 91% perfect rate). Cytorhabdovirus assemblies ranged from 14,784–14,812 bp against a 14,812 bp reference, confirming near-complete genome recovery. Viroid assemblies (359–368 bp) matched reference lengths within ±5 bp. The fast-track mechanism triggered for 15 of 22 assemblies, reducing computation time by bypassing intermediate fusion steps when early-stage skeletons met quality thresholds.

### 3.5 Evolutionary Analysis

Stage 6 post-hoc analysis was completed for all 5 high-confidence viruses. For Cytorhabdovirus sp. 'lycii', Kruskal-Wallis testing revealed significant differences in dN/dS across the 6 viral genes (H = 11.86, p = 0.037), with the G (glycoprotein) gene showing the highest median dN/dS (0.231), consistent with its role as the primary host-interacting protein under diversifying selection. The remaining genes showed strong purifying selection (dN/dS range: 0.039–0.122).

Stage 7 (Capheine) extracted 6 CDS sequences (N: 1,515 bp, P: 1,548 bp, P4: 720 bp, M: 861 bp, G: 627 bp, L: 6,783 bp) from the OR489165.1 GenBank record. HyPhy FEL and MEME analyses identified 10 statistically significant positive selection sites (ω > 1, p < 0.05) distributed across G (positions 49, 194), N (367, 449), P (144, 416), and L (1, 2, 3, 1047) proteins.

Stage 8 (similarity analysis) generated pairwise SDT heatmaps for 3 viruses with ≥3 consensus sequences, revealing near-identical Potato spindle tuber viroid sequences across multiple samples (NT identity > 99.5%), consistent with clonal propagation through vegetative cuttings.

### 3.6 Performance Characteristics

The complete 10-stage pipeline executed in approximately 40 minutes (excluding Stage 7 HyPhy analysis which required ~13 minutes for 6 CDS genes) on a server with 60 CPU threads and 256 GB RAM. Memory usage peaked at approximately 32 GB during Stage 4 assembly. The most computationally intensive stages were Stage 4 (22 assembly tasks, ~22 minutes) and Stage 7 (HyPhy, ~13 minutes). Stages 1–3 and 5–6 completed within 5 minutes collectively, demonstrating the efficiency of the Salmon pseudo-alignment and parallel variant analysis architecture.

---

## 4. Discussion

KVAP addresses a critical gap in plant virus bioinformatics by providing an integrated, reproducible framework that spans the complete known-virus analysis lifecycle. Several design decisions distinguish this pipeline from existing approaches.

**Profile-driven reproducibility**: The YAML profile system ensures that all parameters are captured in a single configuration file, enabling exact reproduction of analyses and facilitating collaboration. This design choice was motivated by the recognition that parameter inconsistency across analysis steps is a common source of irreproducibility in multi-tool bioinformatics workflows.

**Virus-specific quality control**: The Poisson Ratio filter and dual-track validation strategy are specifically designed for viral genomes, where small genome size and high mutation rates can produce misleading alignment signals. Traditional quality metrics designed for eukaryotic genomes (e.g., mapping quality, insert size distribution) are insufficient for viral analyses.

**Automatic K-selection in PCA clustering**: The Silhouette Score-based automatic K-selection eliminates subjective decisions about the number of evolutionary sub-lineages. This is particularly important for viral populations where the true number of sub-populations is unknown and may vary across genomic regions. The 3D PCA extension reveals hidden substructure that may be masked in standard 2D projections.

**Limitations**: Several limitations should be acknowledged. First, the pipeline is designed for known virus analysis and does not perform de novo virus discovery; integration with tools such as VirSorter or Cenote-Taker 2 would extend its applicability. Second, Stage 7 (Capheine) requires CDS annotations, which are not available for viroids and other non-coding viral elements. Third, Stage 9 (DVG detection) requires raw FASTQ files with quality scores, limiting its applicability to datasets where only processed FASTA files are available. Fourth, the pipeline has been validated primarily on plant RNA virus datasets; performance on DNA viruses or animal virus samples may require parameter adjustments. Fifth, the current implementation does not support cloud-based execution or containerization (Docker/Singularity), which would improve portability.

---

## 5. Conclusion

KVAP provides a comprehensive, automated solution for known virus characterization from high-throughput sequencing data. By integrating 20 analysis modules under a unified orchestrator with profile-driven configuration, the pipeline reduces the barrier to reproducible virus analysis while maintaining the flexibility to accommodate diverse research questions. The modular architecture allows users to execute individual stages independently or run the complete workflow with a single command. Future development will focus on containerization, cloud deployment, and integration with de novo virus discovery tools.

---

## Data Availability

The pipeline source code is available at [repository URL]. The validation dataset (7 Lycium transcriptome samples) is available from NCBI SRA under accession [pending]. Reference virus genomes were obtained from NCBI RefSeq (accessed 2024-2025).

## Author Contributions

[To be completed]

## Funding

[To be completed]

## Conflict of Interest

The authors declare no conflict of interest.

## References

1. Cingolani P, Platts A, Wang LL, et al. A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff. *Fly*. 2012;6(2):80-92.
2. Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. *arXiv*. 2012;1207.3907.
3. Grubaugh ND, Gangavarapu K, Quick J, et al. An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. *Genome Biology*. 2019;20(1):8.
4. Kosakovsky Pond SL, Poon AFY, Velazquez R, et al. HyPhy 2.5—A Customizable Platform for Evolutionary Hypothesis Testing Using Phylogenies. *Molecular Biology and Evolution*. 2020;37(1):295-299.
5. Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. *Nature Methods*. 2012;9(4):357-359.
6. Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP. Maftools: efficient and comprehensive analysis of somatic variants in cancer. *Genome Research*. 2018;28(11):1747-1756.
7. Nelson CW, Moncla LH, Hughes AL. SNPGenie: estimating evolutionary parameters to detect natural selection using pooled next-generation sequencing data. *Bioinformatics*. 2015;31(22):3709-3711.
8. Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods*. 2017;14(4):417-419.
9. Routh A, Johnson JE. Discovery of functional genomic motifs in viruses with ViReMa—a Virus Recombination Mapper—for analysis of next-generation sequencing data. *Nucleic Acids Research*. 2014;42(2):e11.
10. Tithi SS, Aylward FO, Jensen RV, Zhang L. FastViromeExplorer: a pipeline for virus and phage identification and abundance profiling in metagenomics data. *PeerJ*. 2018;6:e4227.
