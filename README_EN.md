# MMPV v3.1 — Massive Meta-mining of Plant Viruses

**Large-scale meta-mining and analysis platform for plant viruses**

中文 | **English**

> A complete closed loop for high-throughput identification of plant and other viruses from public sequencing data — covering public data mining, de novo virus discovery, and in-depth quantification / variant / evolutionary analysis of known viruses.

[![Python](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.2-blue.svg)](https://www.r-project.org/)
[![pixi](https://img.shields.io/badge/pixi-enabled-green.svg)](https://pixi.sh/)
[![BioConda](https://img.shields.io/badge/bioconda-supported-brightgreen.svg)](https://bioconda.github.io/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## 📖 Platform Wiki / Documentation

**20-page full documentation** (architecture overview · stage-by-stage interpretation of all seven pipelines · `--help` parameter library for 198 core scripts · audit system · deployment & operations):
👉 [Platform Wiki](https://github.com/zhangwenda0518/MMPV-RNA/wiki) — start from [Home](https://github.com/zhangwenda0518/MMPV-RNA/wiki/Home)

---

## Seven Pipelines

```
public_metadata_pipeline/     → Data acquisition (SRA/GSA search & download & conversion + host reference genome acquisition/indexing)
       │
       ├── Sequencing data: SRA/GSA search → metadata cleaning → bulk download → visualization → SRA→FASTQ.GZ
       └── Host reference: multi-channel download or user FASTA → Kraken2/Bowtie2/HISAT2/Minimap2 four-index build
              │
              ├──▶ endogenous_virus_pipeline/ → Endogenous viral element (EVE) forward screening [ARCHIVED 2026-10-11, immature methodology]
              │      sliding-window blastx vs viral reference proteins → locus merging →
              │      plant/viral dual-side bitscore verdict → RVDB deep assignment
              │      (input is host genome FASTA; independent from the reads pipeline;
              │       moved to archive/20261011_endogenous_virus_pipeline/)
              │
              ▼
data_preprocessing_pipeline/  → Data cleaning (QC Fastp→Seqkit→Clumpify → host depletion Kraken2→alignment→Ribodetector + optional BBNorm normalization)
              │
              ▼
virome_discovery_pipeline/    → Virus discovery (de novo)
       │
       ├── assembly → 10-tool parallel identification → layered filtering → COBRA extension
       ├── merged co-assembly → CD-HIT reference-guided clustering → vclust dereplication → 8-tool taxonomy → host prediction
       └── CheckV assessment → two-branch cascaded rescue (A: direct pass → B: VSI extension; C/D archived) → HQ vOTU catalog
              │
              ▼
virome_analysis_pipeline/     → In-depth analysis of known viruses
       │
       ├── rapid quantification (Salmon/Bowtie2) → variant calling (FreeBayes/iVar)
       ├── SnpEff annotation → SNPGenie evolution → 12-step full-length assembly
       ├── full-length similarity panorama (SDT) → DVG/recombination detection (ViReMa)
       └── interactive HTML comprehensive report
              │
              ▼
virome_submission_pipeline/   → Data submission
       │
       ├── topology inference → metadata template → hypothetical protein annotation
       └── suvtk tbl2asn / Sequin tbl2asn dual mode → .sqn submission file
```

**Optional supplementary module (ARCHIVED 2026-10-11 — immature methodology; moved to
`archive/20261011_host_expression_pipeline/`; main chain and other pipelines unaffected):**

```
host_expression_pipeline/     → ⑧ Host expression analysis (host side: virus-positive/negative grouping) [ARCHIVED]
       ├── host CDS extraction (genome+GFF3 → gffread) → Salmon quantification → quantification QC
       ├── virus label derivation (④ detection table → +/- status + continuous load; Core14 tissue merge)
       ├── TEA-GCN differential co-expression (main analysis: dual networks → ΔzScore → network z-QC gate → DEG-restricted Louvain modules; sample-size/imbalance gates)
       ├── DESeq2 differential expression (cross-reference, optional tissue strata) → load × genome-wide scan (dose validation)
       └── HTML summary report    ← core chain extracted from the goji-virome differential-expression project
```

> See `archive/20261011_host_expression_pipeline/README.md` and `ARCHIVED_20261011.md`
> (archival reason and revival steps) in the same directory. All inputs come from existing
> main-chain artifacts (② clean reads, ① host reference/Core14, ④ detection table); outputs
> land in `08_HostExpression/`.

---

## Unified I/O Layout (v3.1)

The five pipelines now use the `standard` layout by default — **one independent output root per pipeline, numbered independently within each root**.
Existing projects (with checkpoints or running jobs) can fall back to the `legacy` layout (v3.0 directory names):

```
<project>/ 01_PublicData/ 02_Preprocessing/ 03_Discovery/ 04_Analysis/ 05_Phylo/ 06_EVE/ 90_Handoff/
           (06_EVE/ = ⑥ EVE, 08_HostExpression/ = ⑧ host-expression output roots; both ⑥⑧ archived 2026-10-11)
```

Pipeline-boundary handoffs are auto-locatable (one command for the ③→④ bridge; ①→⑤ Core14 read natively):

```bash
export MMPV_IO_LAYOUT=legacy            # existing projects only; new projects default to standard
python virome_discovery_pipeline/utils/discovery2analysis.py \
    --from-discovery <project>/03_Discovery --output_prefix <project>/90_Handoff/analysis
```

For the full project tree, the five boundary contracts, the legacy↔standard mapping table and migration rules, see
**`doc/IO_LAYOUT_DESIGN.md`**.

---

## Quick Start

### 1. One-command deployment

```bash
git clone https://github.com/zhangwenda0518/MMPV-RNA.git
cd MMPV-RNA

# Install all dependencies (pixi; every package verified to exist on its channel)
pixi install

# NOTE: penguin / refineC / flye / PhaBOX2 / CHERRY / VirHunter / VirBot / ViraLM
#       must be installed manually (no bioconda package) — see the manual section at the
#       end of pixi.toml and doc/MMPV_dependencies.md
# NOTE: BEAST/RDP5/PAML etc. for virome_phylo_pipeline are deployed separately per
#       doc/MMPV_dependencies.md A12

# Optional standalone ViraLM environment
conda env create -f envs/viralm.yaml -n viralm

# Download reference databases (see the database table below)
```

> Host genome acquisition supports multiple channels: `--download-source auto|datasets|ngd|ftp|gget` (auto = fall back from datasets to ncbi-genome-download / direct FTP; gget uses Ensembl/Ensembl Plants, suitable for non-model plants missing from NCBI). For your own genome, pass `--genome-fasta` to skip downloading.

### 1b. Data preprocessing (host reference build + QC → host depletion, two pipelines)

```bash
# QC: Fastp → Seqkit → Clumpify
python data_preprocessing_pipeline/clean-data.py -i raw_fastqs/ -o out/00a_CleanData/ -j 20 -t 16

# Host reference: NCBI multi-channel download → four indices (or --genome-fasta for your own genome)
python public_metadata_pipeline/build_host_pipeline.py --species "Solanum lycopersicum" --taxid 4081 --stage all

# Host depletion + rRNA removal: Kraken2 → alignment → Ribodetector
python data_preprocessing_pipeline/host_depletion.py -k kraken2_db/ -x bowtie2_idx/ -I out/00a_CleanData/ -O out/00b_HostDepletion/ --tool bowtie2 --seq-type rna-short --rrna -t 40 -j 10

# Report + handoff list: preprocessing_summary.tsv + assembly_ready.list + HTML
# (one-shot --stage all = clean → deplete → report)
python data_preprocessing_pipeline/data_preprocessing.py --stage report --output_dir out/

# Optional: coverage normalization (only before co-assembly when depths differ greatly; not run by default → out/00c_BBnorm/)
python data_preprocessing_pipeline/data_preprocessing.py --stage bbnorm --output_dir out/ -t 16 -j 4

# Leaner alternative: the public-data pipeline `--stage all` already includes hostref (host genome download)
# + hostdb (four-index build) and report — search species = host species, one shot:
# python public_metadata_pipeline/public_data_pipeline.py --species "..." --taxid ... --stage all

# From .sra to clean reads in one command (optional --bbnorm):
# python public_metadata_pipeline/preprocess_unified.py --sra-dir ./sra/ --outdir out/ --bbnorm \
#     --kraken2-db kraken2_db/ --step2-index bowtie2_idx/host -t 40 -j 4
```

### 2. Discovery pipeline

```bash
pixi run discovery-downstream \
    --output_dir /data/out/ \
    --input_reads /data/00b_HostDepletion/ \
    --host_db /db/hostdb/ \
    --virus_db /db/virus-db/ \
    --checkv_db /db/checkv-db-v1.7/ \
    --host-filter Plant \
    --coassembly \
    -t 120 -j 20
```

### 3. Endogenous virus (EVE) screening [ARCHIVED 2026-10-11, immature methodology]

> ⑥ moved to `archive/20261011_endogenous_virus_pipeline/` (reason and revival steps in
> its `ARCHIVED_20261011.md`). Commands below are the pre-archive usage record; not
> runnable until the module is revived.

```bash
# Single genome: host reference genome → EVE loci (database paths read from pipeline_config.yaml)
python archive/20261011_endogenous_virus_pipeline/eve_screen.py \
    -g host_reference/genome/all.genome.uniq.fasta \
    -n Solanum_lycopersicum -o out/eve/ -t 40

# Batch: batch.tsv, one NAME<TAB>genome-path per line (.tar.gz bundles accepted, auto-unpacked in Stage 1;
#       NAME must be unique — duplicates abort before start, they would share 01_Loci/<NAME>/ products)
python archive/20261011_endogenous_virus_pipeline/eve_screen.py \
    -B batch.tsv -o out/eve/ -t 40 -J 5 --fast
# -J parallel genomes, -t threads per genome; total cores = J x t (256 cores: -J5 -t40)

# Resume: stage products must be complete to count as done; re-issue the same command to resume
# Summary: kingdom_summary.tsv / family_by_genome.tsv (--merge to re-run the summary alone)
# Sequence extraction: samtools faidx (requires samtools >= 1.11; pre-flight check fails early otherwise)
#          index written to the output dir; the reference-genome dir does not need to be writable
# Anti-hang: --cmd-timeout 7200 (per-command timeout in seconds, 0 = unlimited)
```

To tell "host gene vs endogenous viral fossil", and to classify candidate contigs from virome assembly as
true DNA viruses vs degraded EVEs, two modules work together: `eve_screen.py` performs locus discovery on the
host-genome side, and `eve_distinguish/` performs candidate-side discrimination (run after the discovery pipeline):

```bash
# DNA-virus candidate vs EVE discrimination (see archive/20261011_endogenous_virus_pipeline/eve_distinguish/README.md)
bash archive/20261011_endogenous_virus_pipeline/eve_distinguish/run_all.sh \
    -D /path/to/discovery_out -A /path/to/1kp/assemblies -t 32
#   output: dna_vs_eve_filter.tsv (one action per candidate; downstream filters by column)
```

### 4. Analysis pipeline

```bash
# Rapid quantification
python virome_analysis_pipeline/auto_known_virus.py --stage detect \
    --reads_dir /data/out/00b_HostDepletion/ \
    --reference /db/ref.fasta --ref_info /db/ref_info.tsv \
    --tool salmon -t 40 -j 4

# Full workflow
python virome_analysis_pipeline/auto_known_virus.py --stage all \
    --reads_dir /data/out/00b_HostDepletion/ \
    --reference /db/ref.fasta --ref_info /db/ref_info.tsv \
    --snpeff --snpgenie -t 40 -j 4
```

### 5. Public data pipeline

```bash
# Search public SRA/GSA data for a species
python public_metadata_pipeline/public_data_pipeline.py \
    --species "Solanum lycopersicum" --taxid 4081 \
    --stage search info down plot

# Build host reference genome indices
python public_metadata_pipeline/build_host_pipeline.py \
    --species "Solanum lycopersicum" --taxid 4081 \
    --stage all --threads 30
```

### 6. Submission pipeline

```bash
# One-shot NCBI submission files from 08_Rescue
python virome_submission_pipeline/submission_pipeline.py \
    --work-dir $OUT/08_Rescue/ \
    --run-title my_plant_virome \
    --mode both \
    --suvtk-db ~/database/virus-db/suvtk_db/ \
    -t 40

# Desktop GUIs moved to the external local platform; a standalone GUI remains inside the pipeline:
# python virome_submission_pipeline/submission_gui.py
```

---

## Public Metadata Pipeline Stages

### SRA/GSA public data acquisition

| # | Stage | Script | Function |
|---|------|------|------|
| 1 | `search` | gsa_sra.search.py | Dual-engine species search (NCBI SRA + CNCB GSA) → SRA_GSA_Merged_Final.csv |
| 2 | `info` | gsa_sra.info.py | SRA XML parsing + GSA crawler + AI metadata cleaning → Global_Unified_Metadata_Core14.csv (14 columns) |
| 3 | `down` | gsa_sra.down.py | Dual-protocol download (aria2c/wget/prefetch) → raw FASTQ/SRA |
| 4 | `convert` | sra2fastx.py | Bulk SRA → FASTQ.GZ via fasterq-dump (falls back to fastq-dump) |
| 5 | `plot` | gsa_sra.plot.py | Publication-grade 6-panel visualization (time/tissue/region/institute) |
| 6 | `hostref` | download_host_genome.py | Host genome download (datasets/ngd/FTP/gget; `--host-fasta` to skip) → `host_reference/genome/` |
| 6 | `hostdb` | build_hostbase.py | Kraken2/Bowtie2/HISAT2/Minimap2 four-index build → `host_reference/hostdb/` (reuses --species/--taxid) |
| 7 | `report` | generate_report.py | Interactive HTML report + `sample_handoff.csv` handoff list (Run→FASTQ path + metadata; hands off to the cleaning pipeline) |

### Host reference database build

| # | Stage | Script | Function |
|---|------|------|------|
| 1 | `genome-down` | public_metadata_pipeline/download_host_genome.py | Multi-channel reference genome download (datasets/ngd/FTP/gget · or user --genome-fasta) + GFF3 + organelle genomes |
| 2 | `hostdb` | public_metadata_pipeline/build_hostbase.py | Kraken2 + Bowtie2 + HISAT2 + Minimap2 four-index build |

### Unified preprocessing entry (one command from .sra to clean reads)

| # | Stage | Script | Function |
|---|------|------|------|
| 1 | `convert` | sra2fastx.py | .sra → FASTQ.GZ (optional; `--skip-convert` if FASTQ exists) |
| 2 | `clean` | data_preprocessing_pipeline/clean-data.py | Fastp → Seqkit → Clumpify → `00a_CleanData/` |
| 3 | `hostref` | download_host_genome.py + build_hostbase.py | Optional: host genome multi-channel download + four-index build |
| 4 | `deplete` | data_preprocessing_pipeline/host_depletion.py | Kraken2 → alignment → Ribodetector → `00b_HostDepletion/` |
| 5 | `bbnorm` | data_preprocessing_pipeline/run_bbnorm.py | **Optional** (`--bbnorm`, off by default): coverage normalization → `00c_BBnorm/` |

Entry: `python public_metadata_pipeline/preprocess_unified.py --sra-dir ... --outdir ...`
(output names match `data_preprocessing_pipeline` / the discovery pipeline, directly usable as `--input_reads`)

---

## Discovery Pipeline Stages

| # | Stage | Script | Function |
|---|------|------|------|
| 1 | `clean` | data_preprocessing_pipeline/clean-data.py | Fastp QC + Seqkit FASTA conversion + Clumpify cluster-reorder (better compression) |
| 2 | `deplete` | data_preprocessing_pipeline/host_depletion.py | Kraken2 → Bowtie2/HISAT2/Minimap2 → rRNA removal |
| 3 | `assembly` | assembly_pipeline.py | MEGAHIT / rnaviralSPAdes / Penguin assembly |
| 4 | `identification` | virus_identification.py | 10-tool parallel identification (Genomad + Diamond BLASTX + RdrpCatch + ViraLM + VirBot + VirSorter2 + ViralVerify + VirHunter + Metabuli + viroid BLASTN) |
| 5 | `filter` | filter_virus.py + graph_propagation.py | UniProt-strict high-confidence filtering (CDD + UniProt two layers) + optional graph sub-step: second-opinion adjudication of no-evidence keepers (rescue/demote) → COBRA input; idea after GMW, 2026 (off by default, enable via graph.enabled) |
| 6 | `cobra` | cobra_pipeline.py | BWA-MEM2 → CoverM → COBRA-Meta per-sample extension (candidate source = filter main-chain mirror, restoring the filter→COBRA contract) |
| 7 | `merge` | (built-in, Flye) | Multi-sample Flye co-assembly (co-assembly mode) |
| 8 | `cluster` | cluster_pipeline.py | CD-HIT reference-guided pre-clustering + vclust Leiden clustering |
| 9 | `taxonomy` | virus_classifier.py + R | 8-tool taxonomy (Genomad + Metabuli + CAT + Diamond LCA + VITAP + MMseqs2 + ACVirus + vConTACT3) weighted voting → 8-rank taxonomy |
| 10 | `host` | run_host_prediction.py | Three-tier host prediction (ICTV > RNAVirHost > PhaBOX2) + plant-specific post-filter |
| 11 | `checkv` | (built-in) | CheckV completeness assessment |
| 12 | `rescue` | rescue_pipeline.py | Two-branch cascaded rescue (CheckV ≥90% direct pass → Virseqimprover read extension); former branches C/D archived (archive/20261011_rescue_branch_cd/) |
| 13 | `analysis` | (built-in) | Rescue-evidence integration + virome downstream analysis |
| 14 | `analysis_verify` | (built-in) | 09b analysis re-verification (HMM/CT3 evidence scan) |
| 15 | `report` | report_pipeline.py | TSV summary + Sankey diagram + interactive HTML report |

> Note: the `clean` / `deplete` stages live in the top-level pipeline `data_preprocessing_pipeline/`
> (with a clean→deplete one-shot entry); host reference genome acquisition/indexing lives in
> `public_metadata_pipeline/` (`virome_pipeline.py` can still orchestrate the whole chain).
>
> Note: **BBNorm coverage normalization is no longer a discovery-pipeline stage** (optional before
> co-assembly). It belongs to `data_preprocessing_pipeline/run_bbnorm.py`; both entries skip it by
> default: `data_preprocessing_pipeline/data_preprocessing.py --stage bbnorm` (not in `--stage all`)
> and `public_metadata_pipeline/preprocess_unified.py --bbnorm`. Output is `00c_BBnorm/`, directly
> usable as `virome_pipeline.py --input_reads`.

---

## Analysis Pipeline Stages

| # | Stage | Script | Function |
|---|------|------|------|
| 1 | `detect` | batch_virus_depth.py | Salmon/Kallisto/Bowtie2 rapid quantification + Poisson filtering |
| 2 | `filter` | utils/filter_summary.py | High-confidence filtering |
| 3 | `variants` | batch_virus_variants.py | FreeBayes/iVar/LoFreq variant calling + consensus building |
| 4 | `post` | virus_vcf_pipeline.py | VCF visualization + SnpEff macro + MAF waterfall + SnpGenie dN/dS |
| 5 | `full` | batch_virus_full.py → virus-full.py | **Full-length genome build (OmniVirusAssembler V9.0)**: 12 refinement steps — SHIVER-like Divine Fusion solid backbone → read-level iterative polishing (single-base accuracy) → gmcloser + abyss-sealer dual-engine gap closing → circularization detection (Circular=True) |
| 6 | `extract` | utils/extract_full_fasta.py | Longest-contig extraction + reference N-fill (RNA: fill only if N<5%) |
| 6 | `similarity` | virus_auto_pipeline.py | Full-length similarity heatmap + hierarchical clustering (MAFFT --auto) |
| 7 | `dvg` | batch_virema_dvg.py | ViReMa DVG/recombination detection + Circos plots |
| 8 | `report` | generate_pipeline_report.py | Interactive HTML comprehensive report + AI interpretation prompts |

> Note: positive-selection analysis (`capheine` / HyPhy FEL·MEME·BUSTED·PRIME) moved to
> `virome_phylo_pipeline/` in v3.0 (the `capheine` stage in STAGE_ORDER). The 9-stage definition
> of this pipeline is in the `auto_known_virus.py` docstring.

---

## Submission Pipeline Stages

| # | Stage | Script | Function |
|---|------|------|------|
| 1 | `topology` | viral_topology.py | BWA-MEM2 end-read alignment → genome topology (circular/linear) |
| 2 | `metadata` | unified_metadata.py | Unified metadata template (source.src + features + organism) |
| 3 | `hypothetical` | analyze_hypothetical.py | Hypothetical protein annotation (5 tools: HHsuite/Diamond/DeepLoc/PSORTb/TMHMM) |
| 4 | `sequin` | sequin_builder.py | Sequin .tbl construction (Cenote-Taker3 style) |
| 5 | `submit` | submission_pipeline.py | End-to-end orchestration: suvtk tbl2asn / Sequin tbl2asn dual mode → .sqn |
| 6 | `gui` | submission_gui.py (in-pipeline) | PySide6 desktop GUI: interactive editing/validation/export |
| 6 | `report` | report_html.py | Interactive HTML full-table editing report |

---

## Directory Structure

```
MMPV-RNA/
├── data_preprocessing_pipeline/    # Data cleaning (QC clean-data.py → host depletion host_depletion.py; one-shot entry)
├── virome_discovery_pipeline/      # Virus discovery (de novo, 15 stages)
│   ├── virome_pipeline.py          # Main orchestrator
│   ├── doc.md                      # Full pipeline documentation
│   └── utils/                      # Helpers (assembly/identification/COBRA stats, Sankey)
│
├── virome_analysis_pipeline/       # Known-virus deep analysis (9 stages)
│   ├── auto_known_virus.py         # Analysis orchestrator
│   ├── doc.md                      # Full pipeline documentation
│   └── utils/                      # Helpers
│
├── public_metadata_pipeline/       # Public data acquisition (8 stages)
│   ├── public_data_pipeline.py     # Orchestrator (search/metadata/download/visualization)
│   ├── build_host_pipeline.py      # Host reference orchestration (multi-channel download/user FASTA → four indices)
│   ├── preprocess_unified.py       # Full-chain preprocessing entry (convert→clean→hostref→deplete→[bbnorm])
│   ├── doc.md                      # Full pipeline documentation
│   └── utils/                      # Shared tool modules
│
│   ├── main.py                     # GUI entry
│   ├── controllers/                # Search bridge / AI completion / metadata control
│   ├── models/                     # Data-store models
│   ├── views/                      # Main window / search view / table / visualization / detail panel
│   └── utils/                      # Helpers
│
├── virome_phylo_pipeline/          # Phylogenetics/evolution/population genetics (18 stages, 8 modules)
│   ├── phylo_pipeline.py           # Main pipeline (BEAST dating + pypopart + positive selection)
│   └── doc/ authoritative docs (RUN_GUIDE/STAGE_REFERENCE/METHODS)
├── virome_submission_pipeline/     # Submission pipeline (GenBank/CNCB)
│   ├── submission_pipeline.py      # Main orchestrator
│   ├── sequin_builder.py           # Sequin builder
│   └── ...                         # metadata/report/topology analysis
│
├── archive/20261011_endogenous_virus_pipeline/  # [ARCHIVED 2026-10-11] Endogenous viruses (EVE): host-genome screening + candidate discrimination (was endogenous_virus_pipeline/, immature methodology)
│   ├── eve_screen.py               # Screening orchestrator (three stages: discover→verdict→annotate→merge)
│   ├── eve_genome_scan.py          # Per-genome worker (stage-product checkpointing)
│   ├── eve_scan_core.py            # Pure-logic core (coordinate restore/locus merge/verdict/summary)
│   ├── tests/                      # 142 unit tests (core logic/discrimination/invariants)
│   ├── doc.md                      # Full pipeline documentation
│   └── eve_distinguish/            # ★ post-run module: DNA-virus candidate vs EVE discrimination (v6.2)
│       ├── run_all.sh              # One-shot entry (0 panel → 1 db → 2 blastx → 3 structure/domain → 4 host → 5/6 verdict)
│       ├── s1_decay_scan.py        # Structural decay (stop-code enrichment + distribution profile)
│       ├── s2_domain_scan.py       # Domain architecture (MP/CP/AP/RT/RH component scoring)
│       ├── s2b_locus_scan.py       # Candidate blastn back to host → locus architecture
│       ├── s3_verdict.py           # decay × architecture 2-D verdict + host veto
│       ├── s4_filter.py            # verdict → action (REMOVE/MOVE_EVE/KEEP_virus/REVIEW)
│       ├── build_panel.py          # Reference panel build (ships panel.fasta/baits.fa)
│       └── README.md               # Methods & pitfalls (with measured numbers)
│
│
│
├── biosoft/                        # Third-party tools (scripts/JARs, no compilation)
│   ├── VirBot/VirBot.py
│   ├── virhunter/predict_cpu.py + weights/
│   ├── ViReMa/ViReMa.py
│   └── snpEff/snpEff.jar + config + scripts/
│
├── doc/                            # Documentation (17 docs + methods template + dependency list)
│   ├── METHODS_TEMPLATE.md         # Combined paper-grade Methods for all pipelines (ready for SCI submission)
│   └── MMPV_dependencies.md        # External dependency list (tools/databases/paths)
├── scripts/                        # Dev-time one-off script archive (non-core; see scripts/README.md)
│   ├── blacklist/                  # Host blacklist build/sync/validation (bl_layer1~4)
│   ├── db_query/                   # Result-db ad-hoc queries (q_*)
│   ├── patch/                      # One-off patches already merged into pipelines
│   ├── audit/                      # Result audits/cross-checks
│   ├── pilot_rvdb/                 # RVDB pilot analyses
│   ├── mapper_bench/               # Aligner benchmarks
│   ├── utils/                      # md2docx / server-side launchers
│   └── check_script_refs.py        # Pipeline script reference integrity check
├── paper/                          # SCI manuscript (English full-workflow draft)
├── reports/                        # Generated HTML reports (gitignored)
├── envs/viralm.yaml                # Standalone ViraLM conda environment
├── pixi.toml                       # One-command deployment (conda-installable deps; manual tools in tail comments)
├── pipeline_config.yaml            # Pipeline configuration
├── SOFTWARE_VERSIONS.txt           # All third-party software versions
└── README.md
```

---

## Key Features

- **10-tool parallel virus identification**: Genomad + Diamond BLASTX + VirSorter2 + ViralVerify + VirHunter + Metabuli + RdrpCatch + ViraLM + VirBot + viroid BLASTN
- **CD-HIT reference-guided pre-clustering**: integrates complete ICTV/NCBI genomes as references, attaching fragment contigs to known species
- **Two-branch cascaded rescue**: A: CheckV≥90% direct pass → B: Virseqimprover read extension — only evidence-backed products are kept (former C pass-through filter / D zero-yield archived: archive/20261011_rescue_branch_cd/)
- **Host filtering**: pre-filters by host category before rescue, saving 70%+ compute
- **pixi one-command deployment**: 100 conda packages with exact versions via `pixi install`
- **Checkpoint/resume**: every script supports `--resume/--force` for safe interruption and recovery at scale
- **Taxonomy-based novelty calling**: novel-virus judgment from taxonomy completeness without BLASTN
- **Complete closed loop**: from public data mining to submission-ready figures, all in one platform
- **Metadata GUI**: PySide6 desktop app for SRA/GSA metadata search, filtering, table browsing and charting
- **Submission pipeline**: GenBank/CNCB submission automation with Sequin building and hypothetical-protein analysis
- **Report template layer (2026-10)**: report HTML rendered through Jinja2 (`mmpv_common/report_template.py` + `templates/`: base skeleton + component macros); process panel / project master report / submission report migrated — auto-escaping for metadata values, fail-loud strict templates, fully offline with zero external assets

---

## Required Databases

| Database | Purpose | Source |
|--------|------|------|
| CheckV DB | Viral completeness assessment | https://bitbucket.org/berkeleylab/checkv/ |
| geNomad DB | Virus identification | https://zenodo.org/records/14828026 |
| RVDB | Viral reference sequences | https://fzer.github.io/rvdbtools/ |
| NR (diamond) | Protein filtering | NCBI nr |
| NCBI virus ref | BLAST reference | NCBI virus |
| ViralVerify HMM | HMM validation | ViralVerify |
| VirSorter2 DB | Virus taxonomy | VirSorter2 |
| ViraLM DB | DNABERT-2 identification | Google Drive (gdown) |
| Kraken2 + align DB | Host depletion | Built by this platform (build_host_pipeline.py) |
| PhaBOX2 DB | Host prediction | PhaBOX2 |
| Plant/viral split DB (hybrid) | EVE dual-side verdict | SwissProt split + id2div mapping (⑥ archived 2026-10-11) |

**Full list & deployment**: all pipelines depend on 50 databases (25 required / 25 optional, ~710 GB total).
Per-database download/indexing commands are in `DATABASE_SETUP.md`; the **reproducible centralized deployment**
(shared DB root + symlink compatibility layer + environment variables + verification) is in `doc/DB_DEPLOYMENT.md`, with the companion tool:

```bash
python scripts/db_setup/mmpv_db.py list                                 # print the 50-item dependency list
python scripts/db_setup/mmpv_db.py check                                # pre-deployment check
python scripts/db_setup/mmpv_db.py verify                               # verify existing databases
python scripts/db_setup/mmpv_db.py down --apply                          # download/build auto-fetchable databases
python scripts/db_setup/mmpv_db.py adopt --db-from ~/database --apply   # migrate to the shared DB root
python scripts/db_setup/mmpv_db.py link --apply                          # create legacy-path symlinks
```

Use `down` when moving servers; the RVDB-derived suite (diamond / mmseqs / metabuli / CAT) is built on-site
from upstream RVDB fasta by the bundled `virome_discovery_pipeline/build_virus_db.py` — see deployment §3.2.

---

## Documentation

| Document | Description |
|------|------|
| `virome_discovery_pipeline/doc.md` | Discovery pipeline, every script — parameters/inputs-outputs/result interpretation |
| `virome_analysis_pipeline/doc.md` | Analysis pipeline, every script — parameters/inputs-outputs/result interpretation |
| `public_metadata_pipeline/doc.md` | Public data pipeline + host reference build (incl. four-channel download) |
| `virome_submission_pipeline/` | Submission scripts all carry detailed docstrings |
| `doc/14-pipeline-diagrams.md` | Mermaid flow/architecture diagrams |
| `doc/DB_DEPLOYMENT.md` | Database deployment plan (shared root/symlinks/env vars/verification, reproducible) |
| `SOFTWARE_VERSIONS.txt` | All third-party software versions |
| `pixi.toml` | Dependency configuration (human-readable) |

---

## Citation

Zhang W. et al. **MMPV: Massive Meta-mining of Plant Viruses**. 2026.

> If you use MMPV in your research, please cite the above reference.
