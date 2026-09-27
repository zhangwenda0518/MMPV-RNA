---
name: mmpv-paper-public-data
description: Draft and audit the public-data mining sections of the MMPV SCI manuscript (Methods 2.1-2.2 public-data retrieval + host-reference construction + read cleaning/depletion; Results 3.2 study-queue composition; Fig 2; Table S2). Use when writing or revising the MMPV/Lycium manuscript sections about SRA/GSA retrieval, metadata unification, AI-assisted cleaning, host database building, or the 1,348-sample queue; or when auditing queue-composition claims against pipeline outputs.
---

# MMPV Paper — Public-Data Mining Sections

You own these manuscript sections: **Methods §2.1 (public data retrieval + host reference construction), §2.2 (read cleaning and host depletion), Results §3.2 (study queue), Fig. 2 legend, Table S2**.

## Authoritative sources (read before writing; never write from memory)

| Source | What to take from it |
|---|---|
| `doc/METHODS_TEMPLATE.md` 管道一 | Canonical methods paragraph: five-stage metadata pipeline, dual SRA+GSA retrieval ("biomol rna"), DeepSeek/Kimi arbitration, 14-column core table, dual-channel download (NGDC CRR/aria2c; NCBI SRR/prefetch), four host indices (Kraken2/Bowtie2/HISAT2/Minimap2) |
| `public_metadata_pipeline/doc.md` + script argparse defaults | Exact stage names and parameter values |
| `pipeline_config.yaml` deplete block | Kraken2 confidence 0.2, keep_taxids 10239, vmtouch off, Bowtie2 two-step default |
| `MMPV-RNA-paper/MMPV_Lycium_virome_manuscript_EN_v1.md` §2.1-2.2, §3.2 | Established prose and queue numbers to revise, not contradict |
| `virome_analysis_pipeline/skills/virome-results-analysis/SKILL.md` | Style rules + fixed benchmark anchors |

## Anchors (fixed; never contradict)

- Queue = **1,348 metatranscriptomic samples**, dominated by *L. barbarum* (incl. cv. Ningxia No.1–No.7), *L. ruthenicum*, *L. chinense*, plus solanaceous/non-solanaceous comparators.
- Host-depletion shipped default: **Kraken2 (confidence 0.2, viral taxid 10239 preserved) + Bowtie2 two-step; RiboDetector optional and disabled by default**. The benchmark's HISAT2 numbers (99.8% host removal, 99.3% viral retention) describe Evaluation-1/D3 config — attribute them to the benchmark, never to the shipped default.
- BBNorm normalization is **optional** (target 70×, min 2) and only before co-assembly; it is **not a stage of the discovery pipeline** — it lives in the preprocessing pipeline (`data_preprocessing_pipeline/run_bbnorm.py`, output `00c_BBnorm/`), reachable via `data_preprocessing.py --stage bbnorm` or `preprocess_unified.py --bbnorm`, and is disabled by default in both.
- Manuscript frames MMPV as **five pipelines** (Public Data, Discovery, Known-virus Analysis, Phylodynamics, Submission); older repo docs say four — the manuscript anchor is five.

## Anti-fabrication rules

- Any queue statistic not present in v1 or pipeline outputs (e.g., candidate runs before QC, % runs needing AI cleaning, species×sample breakdown) stays a bracketed placeholder: `[N — TO CONFIRM]` — never estimate.
- Every number cites its source: "(Table S2)", "(derived from `Global_Unified_Metadata_Core14.tsv`)".
- Tool names in English even in Chinese prose; past tense; no interpretation in Results.

## Style

Journal past tense, "we retrieved/we built". Methods paragraph must mention: RNA-library restriction, LLM-assisted cleaning as *arbitration of rule-based extraction* (not as source of truth), resume-capable downloads, and that host DBs are built **per host species** from NCBI RefSeq assemblies.
