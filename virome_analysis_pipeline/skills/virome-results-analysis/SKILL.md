---
name: virome-results-analysis
description: Analyze MMPV-RNA virome_discovery_pipeline outputs into manuscript-grade Results narratives and diagnostic reports, aligned to the MMPV Nature Communications framework (filter-first paradigm, ★/★★/★★★ novelty tiers). Use when reading pipeline_report.html stage summaries, *_summary.tsv, final_integrated_classification.tsv, ensemble_host_summary.tsv, completeness.tsv, vclust_clusters.tsv, cdd_calls.tsv, rescue_evidence_scored.tsv, or any discovery stage output directory (legacy `00a–10_Reports` / standard `03_Discovery` numbering); writing 结果段落/Results sections from virus discovery numbers; auditing funnel consistency, CheckV completeness distribution, clustering redundancy, taxonomy agreement, evidence verdicts; or diagnosing anomalies like candidate-count cliffs, low CDD pass rates, or broken host annotations.
---

# Virome Results Analysis

Analyze outputs of the MMPV-RNA virome discovery pipeline and turn them into
(a) publication-ready Results prose and (b) an audit/diagnostic report.

## When to use

- User provides a pipeline output root (legacy: `00a_CleanData` … `10_Reports`; standard: `03_Discovery/01…13`)
- User asks to "分析结果", "写Results", "解读pipeline报告", check numbers, or audit quality
- Cross-checking funnel counts between stages

## Inputs you must locate

Locate these files under the output root before doing anything:

| File | Purpose |
|---|---|
| `10_Reports/stage_summary.tsv`, `data_summary.tsv`, `assembly_summary.tsv`, `ident_summary.tsv`, `checkv_summary.tsv` | Per-stage counts and metrics |
| `02b_Filter/cdd_calls.tsv` | CDD tiering of candidates |
| `04_CLUSTER/centroids/final_centroids.fasta`, `3_vclust/vclust_clusters.tsv` | vOTU catalog size, cluster size distribution |
| `05_Taxonomy/*/final_integrated_classification.tsv` | Consensus taxonomy + per-rank agreement + confidence |
| `06_HostPrediction/ensemble_host_summary.tsv` | Host assignments by decision tree |
| `07_Checkv/completeness.tsv` | Completeness / MCP / denominator metrics |
| `08_Rescue/**` (branch_a-d, merged/all_HQ.fasta) | Rescue yield per branch |
| `09_Virome_Analysis/integrated_summary.tsv` (legacy 产物, 生成脚本已归档 — 缺失记 N/A), `suvtk_taxonomy/taxonomy.tsv`, `featuretable.tbl` | Structural annotation, R-vs-suvtk cross-check |
| `09b_Analysis_Verify/rescue_evidence_scored.tsv`, `final_judgement_table.tsv`, `class_KEEP.fasta` | Five-layer evidence chain verdicts |
| `viroid_circular_detect*` outputs | Viroid circularization support |

If a file is missing, note it explicitly as `N/A (stage not run)` — never infer the number.

## Core statistics to compute

1. **Reads funnel**: raw → post-fastp → post-Kraken2 → post-alignment → post-rRNA; report per-sample and total retention %.
2. **Assembly stats**: total bp, contig count, N50, max contig length per sample/assembler.
3. **Identification funnel**: assembled contigs → viral candidates (≥1000 bp) → post-UniProt → post-CDD → COBRA-extended → clustered centroids. Report per-tool and intersection sizes (Venn/UpSet if `comparison_plots/` exists).
4. **Cluster stats**: number of vOTUs, singletons, largest clusters, known-vs-novel split from CD-HIT reference-guided step.
5. **Taxonomy**: counts per rank/phylum/family; mean per-rank agreement; confidence distribution.
6. **Host**: assignment source fractions (ICTV / RNAVirHost / PhaBOX2); unassigned rate.
7. **CheckV**: completeness distribution (min/median/max), fraction ≥90%, branch A-D rescue yields.
8. **Evidence chain**: KEEP/REVIEW/DROP counts from `rescue_evidence_scored.tsv`; score distribution; reasons for REVIEW/DROP breakdown.
9. **Viroids**: circular vs linear candidates with identity/coverage stats.

## Standard audits (run every time)

- **Funnel monotonicity**: each downstream count ≤ upstream count. Flag any violation.
- **Threshold sensitivity notes**: report what changes at ANI 0.95/qcov 0.85 vs stricter cutoffs if vclust raw output allows re-derivation.
- **Cross-source taxonomy conflicts**: rows where R consensus disagrees with suvtk taxonomy — list them.
- **CDD pass-rate sanity**: T1/T2/T3 proportions; a very high NON_VIRAL rate after filtering suggests contamination in identification; a near-zero rate suggests the filter is not biting.
- **Length integrity**: `class_KEEP.fasta` min length ≥1000; rescued contigs with excessive N-runs (branch C gap-filling artifacts).
- **Prevalence plausibility**: samples contributing to each HQ vOTU >0 abundance (`frequency_table.tsv`); zero-prevalence catalog entries are suspicious.

## Output format

Produce two deliverables in Markdown:

### 1) `RESULTS_virus_discovery.md`
English, journal-style past-tense paragraphs following the data flow:
preprocessing → assembly → discovery funnel → extension/clustering → taxonomy → host prediction → completeness/rescue → verification/evidence integration → viroids.
Every number cites its source file in parentheses, e.g. "(Table S1; derived from ident_summary.tsv)". Use n = counts, percentages to one decimal. No speculation in Results; move interpretation to a short closing paragraph labeled *Interpretation*.

### 2) `AUDIT_pipeline_results.md`
Structured audit with a severity-tagged table: `[PASS]/[WARN]/[FAIL]` per check above, plus a "Questions for the analyst" section listing anything requiring human judgment (e.g., contradictory host calls, single-sample endemics).

## Manuscript alignment (MMPV-RNA NC framework)

The MMPV-RNA Nature Communications manuscript (`D:\桌面\文章写作\MMPV_platform_NC_v13.md` and later; RNA draft: `MMPV-RNA_NatureCommunications_revised.md`) fixes a narrative framework. Analysis output must stay consistent with it:

**Fixed benchmark anchors** — never contradict these when writing Results:
- Filter-first paradigm: MetaBuli alone under strict UniProt filtering F1 = 0.942 matches a five-tool ensemble unfiltered (F1 = 0.937); nine-tool union raw F1 = 0.827.
- Assembly: rnaviralSPAdes > MEGAHIT by ~4.3% genome fraction alone; merged contig sets gain 6.0–8.4 percentage points; chimeras <3%. metaViralSPAdes "all-or-nothing"; Penguin high-coverage/high-redundancy probe.
- Classification: ACVirus 94.1% Family/Genus accuracy vs MMseqs2 21.8%; VITAP <1% on plant viruses (phage bias).
- Host depletion default per manuscript (aligned to code, 2026-08-27): Kraken2 (confidence 0.2) + Bowtie2 two-step, RiboDetector optional and disabled by default; benchmark Evaluation 1 itself tested HISAT2 variants (99.8% host removal, 99.3% viral retention) — the numbers describe D3 config, while shipped defaults use bowtie2. When auditing a real run, read the run's `run_config.json` for the aligner actually used.

**Three-tier novelty taxonomy** (use these exact thresholds):
- ★ known: ≥90% nucleotide identity AND ≥90% coverage to reference
- ★★ potentially novel: ≥70% protein identity but <90% nucleotide identity
- ★★★ truly novel: <70% protein identity to any reference

**Terminology**: ten-stage discovery workflow; four pipelines (Public Data, Virome Discovery, Virome Analysis, Virome Submission); vOTU catalog from centroids; HQ vOTUs after cascaded rescue. Family-level claims should cite ICTV VMR MSL41-derived reference database.

**Where results data plugs into the paper**: ident_summary/clusters → Fig. survey vOTU counts and novelty tiers; integrated_classification → family composition panels; evidence chain verdicts → novel-virus characterization paragraphs; CheckV/rescue yields → extension/rescue improvement claims (paper reserves these as [TBD] slots — analysis deliverables should fill them with file-cited numbers).

If analysis reveals values inconsistent with the anchors above (e.g., a rerun used different parameters), flag it in AUDIT under `[WARN] manuscript-anchor conflict` rather than bending the narrative.

## Style constraints

- Light-color, colorblind-safe figures only if plotting is requested (Okabe-Ito palette, ≥300 dpi, English labels).
- Never fabricate a count: if a file is absent or empty, mark N/A.
- Version-drift awareness: default parameters correspond to `virome_pipeline.py` v3.x (see METHODS_Virome_Discovery.md in the pipeline repo). If the user's run used overrides recorded in `run_config.json`/`provenance.json`, prefer those values in prose.
- Chinese responses should mirror the same structure if the user writes in Chinese; keep file names and tool names in English regardless.
