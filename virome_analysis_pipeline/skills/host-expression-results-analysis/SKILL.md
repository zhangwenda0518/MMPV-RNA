---
name: host-expression-results-analysis
description: "[MODULE ARCHIVED 2026-10-11 — immature methodology; pipeline moved to archive/20261011_host_expression_pipeline/. Retained only for analyzing historical 08_HostExpression outputs.] Analyze MMPV-RNA host_expression_pipeline (platform ⑧, optional module) outputs into manuscript-grade Results narratives and audit reports. Use when reading 08_HostExpression outputs (02_Quant matrices, 03_QC qc_summary, 04_Labels virus_labels/sample_labels, 05_Coexp diff_coexp_* / network_qc.tsv / group_summary.tsv / modules/, 06_DE deseq2_* results, 07_Enrich enrich_summary, 08_LoadScan genome scans, 09_Report HTML); writing 宿主响应/差异共表达/TEA-GCN/ΔzScore Results sections; auditing coexp gates (min 30 per group, imbalance ≤3×), network z-QC verdicts, DEG-limited Louvain modules, module×load associations, DESeq2 cross-reference, or viral-load genome scans."
---

# Host-Expression Results Analysis

> **⚠️ 模块已归档（2026-10-11，方法不成熟）**：产生这些产物的
> `host_expression_pipeline/` 已移至 `archive/20261011_host_expression_pipeline/`
> （归档原因与恢复方法见其 `ARCHIVED_20261011.md`）。本技能仅保留用于分析
> **历史项目已有的** `08_HostExpression/` 输出；新分析请勿再产出基于该模块
> 的投稿级结论，除非模块恢复。

Analyze outputs of the MMPV-RNA host-expression module (⑧, optional, off the main chain)
and turn them into (a) publication-ready Results prose and (b) an audit/diagnostic report.

Key framing: the module answers "**how does the host respond to viral infection**" in three
mutually-checking layers — 一核两翼: **core** = TEA-GCN differential co-expression
(05_Coexp); **cross-reference** = DESeq2 differential expression (06_DE); **dose
validation** = viral-load × genome-wide Spearman scan (08_LoadScan). Output root
`08_HostExpression/` (same name in both layouts). Inputs come from ①②④ main-chain
products; grouping defaults to `Virus_any` positive vs negative, switchable per virus.

## When to use

- User provides an `08_HostExpression/` output root
- User asks to "分析宿主响应", "差异共表达结果", "TEA-GCN", "载量扫描", audit coexp gates, or interpret Δz hub/module tables
- Drafting manuscript host-response sections (Methods 管道七 in `doc/METHODS_TEMPLATE.md`)

## Inputs you must locate

| File | Purpose |
|---|---|
| `04_Labels/virus_labels.tsv`, `sample_labels.tsv`, `virus_columns.tsv` | Grouping actually used (Virus_any or `<VIRUS>_positive`), group sizes, Tissue column |
| `03_QC/qc_summary.tsv`, `qc_passed_samples.txt` | Samples excluded by mapping-rate/read-count thresholds |
| `02_Quant/matrix/TPM_matrix.tsv` (+NumReads/FPKM/EffectiveLength) | Quantification backing all downstream |
| `05_Coexp/group_summary.tsv` | **Gate record**: per-group n, imbalance, gate verdict (skip/pass) |
| `05_Coexp/network_qc.tsv` | Per-group network z-distribution QC: PASS/LOW/FAIL + Δz gate action |
| `05_Coexp/diff_coexp_{candidates,significant,gene_aggregate,summary}.tsv` | Δz pairs, significant pairs (0.5%/99.5% percentiles), hub aggregation |
| `05_Coexp/modules/module_table_<group>.tsv`, `module_summary.tsv`, `module_load_assoc.tsv` | DEG-limited Louvain modules (z≥1.7 both-DEG edges, res 1.0) + module×load Spearman |
| `06_DE/deseq2_all_results.tsv`, `deseq2_degs.tsv`, `deseq2_strata/` | DESeq2 (pydeseq2) cross-reference, optional tissue stratification |
| `08_LoadScan/<VIRUS>_genome_scan.tsv`, `load_scan_summary.tsv` | Dose validation (BH-FDR; `stratum` column splits all/tissue) |
| `07_Enrich/enrich_summary.tsv`, `enrich_<set>.tsv` | Optional ORA (auto-skipped without annotation) |
| `.checkpoints/`, `host_expression.log` | Which stages actually ran; coexp 7 substeps |

If a file is missing, note it explicitly as `N/A (stage not run / gate skipped / no annotation)` — never infer.

## Core statistics to compute

1. **Cohort**: QC-passed n, positive/negative group sizes, gate outcome; grouping column actually used.
2. **Coexp (core)**: significant Δz pair count + hub Top-N genes; per-group network QC verdicts; module count/sizes; module × viral-load associations (ρ, BH-FDR).
3. **DE (cross-reference)**: DEG count (with thresholds used), tissue-stratified results if present; direction agreement with Δz hubs.
4. **Load scan (dose validation)**: significant genes genome-wide and within-tissue; concordance with coexp/DE hits.
5. **Enrichment (optional)**: significant terms per gene set (DEG / dz_hub / module_*).

## Standard audits (run every time)

- **Gate first**: if `group_summary.tsv` records a gate skip (each group ≥30 AND ratio ≤3× required), 05_Coexp outputs are absent **by design** — report the downgrade path (06_DE + 08_LoadScan), never as a failure.
- **Network QC before Δz**: FAIL group ⇒ Δz not interpretable for pairs involving that group; check `network_qc.tsv` verdicts before quoting Δz numbers.
- **Counting discipline**: significant-pair counts are a ranking/truncation artifact (per-gene Top-N extremes vs percentile thresholds), NOT calibrated test counts — report hub ranking and modules, flag raw pair counts if quoted.
- **Label/quant consistency**: QC-excluded samples must not appear in label groups; label n per group must match `group_summary.tsv`.
- **Module composition**: every module edge must be z≥1.7 AND both endpoints DEG; module member sets ⊆ quantified genes; module×load uses TPM module means (rank-equivalent), DE uses raw counts — never mix the two models in one claim.
- **Grouping口径**: single-virus mode (`--virus GCVA` → `<VIRUS>_positive`) overwrites prior 05_Coexp products; check `group_summary.tsv` for which grouping produced the current tables.
- **Stratified vs pooled**: pooled-cohort summaries are tissue-confounded (ERAD lesson); conclusions prefer tissue-stratified + dose-response cross-validation when Tissue exists.

## Output format

Two Markdown deliverables:

### 1) `RESULTS_host_response.md`
English, journal past tense: cohort & gates → DE overview → coexp core (Δz hubs, modules, load associations) → dose validation → enrichment. Every number cites its source file; interpretation quarantined to a closing *Interpretation* paragraph.

### 2) `AUDIT_host_expression.md`
Severity-tagged `[PASS]/[WARN]/[FAIL]` table for the checks above + "Questions for the analyst" (e.g. modules driven by one tissue, hubs without DE support).

## Style constraints

- Terminology: module ⑧ (optional, off main chain); 一核两翼 core/cross-reference/dose-validation; TEA-GCN; Δz(=z_pos−z_neg); DEG-limited Louvain modules; output root `08_HostExpression/`.
- Limitations to state verbatim when writing manuscript prose: Δz percentiles are a ranking device, not a calibrated test; single-resolution Louvain (res 1.0, seed 42) without bootstrap; DE on counts vs TPM for load scans; tissue confounding guarded by stratification + dose cross-validation.
- Never fabricate counts; absent file → N/A.
- Chinese responses mirror the same structure; keep file names, tool names, and statistic names in English.
- Manuscript anchor: Methods 管道七 in `doc/METHODS_TEMPLATE.md`; module detail in `archive/20261011_host_expression_pipeline/README.md` (module archived 2026-10-11).
