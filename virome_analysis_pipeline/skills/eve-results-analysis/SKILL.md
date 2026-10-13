---
name: eve-results-analysis
description: "[MODULE ARCHIVED 2026-10-11 — immature methodology; pipeline moved to archive/20261011_endogenous_virus_pipeline/. Retained only for analyzing historical 06_EVE outputs.] Analyze MMPV-RNA endogenous_virus_pipeline (platform ⑥, EVE) outputs into manuscript-grade Results narratives and audit reports. Use when reading 06_EVE outputs (01_Loci, 02_Verdict s2_verdict.tsv, 03_RVDB, 04_Summary, kingdom_summary.tsv, kingdom_loci_all.tsv.gz, family_by_genome.tsv), 05_Elements element tables (element_table.tsv, element_regions.tsv, eve_consensus_library.fa, family_quant.tsv), or eve_distinguish outputs (eve_distinguish_verdict.tsv, dna_vs_eve_filter.tsv, standard 06_EVE/eve_distinguish vs legacy 08b_EVE_Distinguish); writing EVE/内源性病毒元件 Results sections; auditing verdict distributions, cross-genome element clustering, or DNA-virus-vs-degenerate-EVE filter consistency with the discovery pipeline."
---

# EVE Results Analysis

> **⚠️ 模块已归档（2026-10-11，方法不成熟）**：产生这些产物的
> `endogenous_virus_pipeline/` 已移至 `archive/20261011_endogenous_virus_pipeline/`
> （归档原因与恢复方法见其 `ARCHIVED_20261011.md`）。本技能仅保留用于分析
> **历史项目已有的** `06_EVE/` 输出；新分析请勿再产出基于该模块的投稿级结论，
> 除非模块恢复。

Analyze outputs of the MMPV-RNA endogenous virus element (EVE) pipeline and turn them into
(a) publication-ready Results prose and (b) an audit/diagnostic report.

Key framing: this module's input is **host reference genomes, not sequencing reads** — it
answers "where in the host genome are viral insertions and what are they", complementing
(not duplicating) the reads-based discovery pipeline. Platform module ⑥; standard output
root `06_EVE/` (legacy runs keep everything under the discovery output root).

## When to use

- User provides a `06_EVE/` output root (or legacy `08b_EVE_Distinguish/` for eve_distinguish)
- User asks to "分析EVE结果", "解读内源性病毒元件", "EVE 筛查报告", audit verdicts, or cross-check eve_distinguish actions against discovery-pipeline verdict tables
- Drafting manuscript EVE sections (Methods 管道六 in `doc/METHODS_TEMPLATE.md`)

## Inputs you must locate

| File | Purpose |
|---|---|
| `01_Loci/<NAME>/<NAME>.s1_best.tsv`, `loci.bed`, `loci.fa` | Stage1 loci + best-hit coarse family (`ref_family`) |
| `02_Verdict/<NAME>/<NAME>.s2_verdict.tsv` | Per-locus verdict: `viral_supported` / `host_like` / `undetermined` + both-side best hits |
| `03_RVDB/<NAME>/<NAME>.s3_rvdb.tsv` | RVDB deep family/genus annotation (empty file = no candidates, valid checkpoint state) |
| `04_Summary/<NAME>_eve_summary.tsv` (+ `_viroid.tsv` optional) | Per-genome 12-column rollup; completion markers |
| `kingdom_summary.tsv`, `kingdom_loci_all.tsv.gz`, `family_by_genome.tsv` | Cross-genome merge: verdict counts / all loci (13 col) / ref_family × genome matrix |
| `05_Elements/element_table.tsv`, `element_regions.tsv`, `eve_consensus_library.fa`, `family_quant.tsv`, `consensus_domains.tsv`, `phylo/` | Element-level layer (only if eve_element_db.py was run) |
| `eve_distinguish/eve_distinguish_verdict.tsv`, `dna_vs_eve_filter.tsv` | Contig-level true-DNA-virus vs degenerate-EVE call + action/action_reason (consumed by discovery evidence chain) |
| `logs/` + batch log | Which stages actually ran; `--fast` / `--overlap` / viroid layer on/off |

If a file is missing, note it explicitly as `N/A (stage not run / element layer not run)` — never infer.

## Core statistics to compute

1. **Loci funnel**: genomes screened → loci found → verdict split (viral_supported / host_like / undetermined). Per-genome and pooled.
2. **Family composition**: `ref_family` distribution (coarse, keyword-derived) corroborated by RVDB `rvdb_stitle` where available.
3. **Elements** (if 05_Elements exists): element count, sizes, cross-genome element families (`n_genomes`), consensus library size, copy-number matrix highlights, chimera flags.
4. **eve_distinguish**: candidate counts by verdict and by action (keep / review / drop-out semantics per `action_reason`); how many discovery candidates were reclassified as EVE-like.
5. **Viroid layer** (optional): circular-vs-linear support stats.

## Standard audits (run every time)

- **Verdict sanity**: a near-100% viral_supported rate suggests the Stage2 split database path was wrong; a near-100% host_like rate suggests Stage1 db drift. Both are `[FAIL]`-worthy.
- **ref_family is NOT taxonomy**: Stage1 `ref_family` is a lightweight keyword mapping — never report it as precise classification; use RVDB annotation for family claims.
- **Merge consistency**: `kingdom_summary.tsv` counts must equal sums over per-genome `*_eve_summary.tsv`; `family_by_genome` row sums must equal loci with that ref_family.
- **Run-mode awareness**: `--fast` cuts loci by ~half (diamond --fast); `--overlap` densifies borders; both change counts — read the batch log before quoting numbers.
- **Cross-pipeline consistency**: `dna_vs_eve_filter.tsv` actions vs discovery `final_judgement_table.tsv` / host-contamination columns — flag any candidate KEEPed upstream but actioned out here (or vice versa) as `[WARN] cross-pipeline conflict`.
- **Naming**: loci IDs are `contig:start-end` — quote them verbatim.

## Output format

Two Markdown deliverables:

### 1) `RESULTS_eve.md`
English, journal past tense: screening scope → loci funnel → verdict split → family composition (RVDB-backed) → elements (if run) → eve_distinguish reconciliation. Every number cites its source file; percentages one decimal; interpretation quarantined to a closing *Interpretation* paragraph.

### 2) `AUDIT_eve.md`
Severity-tagged `[PASS]/[WARN]/[FAIL]` table for the checks above + "Questions for the analyst" (e.g. undetermined-heavy genomes worth manual review, single-genome element families).

## Style constraints

- Terminology: EVE; verdict terms exactly `viral_supported` / `host_like` / `undetermined`; element terms element / element region / consensus library; module ⑥; never call ref_family a "taxonomic assignment".
- Never fabricate counts; absent file → N/A.
- Chinese responses mirror the same structure; keep file names, verdict values, and tool names in English.
- Manuscript anchor: Methods 管道六 in `doc/METHODS_TEMPLATE.md` (merged canonical); per-stage detail in `archive/20261011_endogenous_virus_pipeline/doc.md` (module archived 2026-10-11).
