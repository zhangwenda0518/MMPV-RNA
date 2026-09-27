---
name: mmpv-paper-discovery
description: Draft and audit the virus-discovery sections of the MMPV SCI manuscript (Methods 2.3 discovery pipeline; Results 3.5 novel-virus candidates; Table S3; part of Fig 1). Use when writing/revising manuscript text about the ten-tool identification panel, CDD two-layer filtering, COBRA extension, CD-HIT+vclust clustering, eight-tool taxonomy ensemble, host-prediction decision tree, four-branch rescue cascade, five-layer evidence chain, or the Lycium novel-species candidate counts.
---

# MMPV Paper — Virus-Discovery Sections

You own: **Methods §2.3, Results §3.5, Table S3 caption, Fig. 1 (discovery part)**.

## Authoritative sources

| Source | Take from it |
|---|---|
| `doc/METHODS_TEMPLATE.md` 管道二 | Canonical methods: 16 stages, 10-tool panel, E-value ≤ 1e−5, two-layer FP filtering (UniRef90 E ≤ 1e−3 majority vote; MMseqs2-CDD tiers 1/2/3), COBRA (k-mer 21–141, linkage mismatch 2), metaFlye scaffolding, CD-HIT (ANI ≥ 0.95, qcov ≥ 0.85) + vclust Leiden, 8-tool taxonomy weighted voting, host decision tree, 4-branch rescue (A ≥ 90%; B Virseqimprover; C dc-megablast; D ±15% genus length), 5-layer evidence (score = 0.30·BLASTN + 0.30·BLASTX + 0.40·domain, domain = max(CDD, CT3-HMM viral-profile evidence) → KEEP/REVIEW/DROP), suvtk structural annotation, viroid circularity (terminal self-BLASTN ≥ 95% id, ≥ 90% qcov) |
| `virome_discovery_pipeline/METHODS_Virome_Discovery.md` | Per-parameter detail if the merged template is insufficient |
| `virome_analysis_pipeline/skills/virome-results-analysis/SKILL.md` | Novelty tiers + anchors |
| `MMPV-RNA-paper/MMPV_Lycium_virome_manuscript_EN_v1.md` §2.3, §3.5 | Established prose/numbers |

## Lycium result anchors (fixed)

- *L. barbarum*: 284 rescued contigs → **24 novel-species candidates**, 38 known, 77 distant-relative, 40 host-contamination, 125 unclassified/distant.
- *L. ruthenicum*: 294 contigs → **21 novel-species candidates**, 51 known. *L. chinense*: 4 confirmed known viruses from 130 contigs.
- Representative novel candidate: SRR8191340, 14,939 nt, depth 528.6×, **79.3% nt identity** to *Cytorhabdovirus* sp. 'lycii' OR489165.1 at 78% qcov, 89.7% aa identity over RdRp (pfam00946) → distinct cytorhabdovirus species candidate.
- Novelty tiers (exact thresholds): ★ known ≥ 90% nt id AND ≥ 90% coverage; ★★ ≥ 70% aa id but < 90% nt id; ★★★ < 70% aa id to any reference.

## Anchors / terminology

- "ten-stage discovery workflow"; tools named exactly as in METHODS_TEMPLATE (geNomad, DIAMOND BLASTX, RdRp-Catch, ViraLM, VirBot, VirSorter2, ViralVerify, VirHunter, Metabuli, viroid BLASTN).
- vConTACT3 is **opt-in** (`--tools all` runs 7 taxonomy tools) — if the manuscript lists 8-tool ensemble, keep v1's phrasing "eight-tool ensemble" only when the run actually enabled it; otherwise say "seven-tool default ensemble with optional vConTACT3".
- Verdict terms KEEP/REVIEW/DROP; catalog terms vOTU / HQ vOTU from centroids.

## Anti-fabrication rules

Never invent candidate counts, CDD pass rates, or rescue-branch yields — take them from `final_judgement_table.tsv` / `rescue_evidence_scored.tsv` or keep v1's numbers; anything else stays `[N — TO CONFIRM]`.
