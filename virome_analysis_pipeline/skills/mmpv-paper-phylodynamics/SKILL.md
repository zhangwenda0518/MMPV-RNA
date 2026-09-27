---
name: mmpv-paper-phylodynamics
description: Draft and audit the phylodynamics sections of the MMPV SCI manuscript (Methods 2.8; Results 3.8-3.10; Table 3; phylogeny figures). Use when writing/revising manuscript text about MAFFT/IQ-TREE trees, temporal-signal gating (TreeTime RTT, date-randomization test, BETS), BEAST dating configuration, CTMC/BSSVS phylogeography, pypopart haplotypes, DnaSP/PAML/HyPhy selection, RDP5 recombination, or the PSTVd/Cytorhabdovirus population-genetic results.
---

# MMPV Paper — Phylodynamics Sections

You own: **Methods §2.8, Results §3.8 (Cytorhabdovirus diversity/geography), §3.9 (no molecular clock), §3.10 (recombination/DVG); Table 3; phylogeny figure legends.**

## Authoritative sources

| Source | Take from it |
|---|---|
| `doc/METHODS_TEMPLATE.md` 管道四 | Canonical methods: MAFFT 7.525 (`--auto --reorder`), coverage QC ≥ 85%, IQ-TREE 3.1.2 (`-m MFP`, UFBoot 1000); temporal gate = TreeTime RTT (R², rate β) + DRT (20 replicates, real R² > 95th percentile) + BETS; **BEAST dating only when RTT R² > 0.3** (UCLN, coalescent skyline / constant / birth-death skyline, HKY/GTR+G4, 5 chains × 5×10⁷ states, ESS > 200, 10% burn-in, LogCombiner + TreeAnnotator MCC); CTMC + BSSVS phylogeography (BF ≥ 5), RRT + TempMig + SpreaD3; pypopart (MJN + VirNA directed networks, π, θ, Tajima's D, Fu's Fs, locus-based Fst with permutation); DnaSP 6 (Ka/Ks, Fu & Li's D*/F*, Rm, LD), codeml (PAML 4), HyPhy 2.5 (FEL/MEME/BUSTED); RDP5 (7 methods, ≥ 3-method high-confidence events, mask-or-delete tree strategies, DnaSP Rm cross-check) |
| `virome_phylo_pipeline/METHODS.md` + `PIPELINE_FLOW.md` | Stage detail |
| `MMPV-RNA-paper/MMPV_Lycium_virome_manuscript_EN_v1.md` §2.8, §3.8-3.10 | Established prose + result numbers |

## Lycium result anchors (fixed)

- *Cytorhabdovirus* sp. 'lycii': 127 coverage-QC genomes, 15,591 aligned sites → **124 haplotypes, Hd = 0.9996, π = 1.022% (k = 159.29, S = 2,124, θ-W = 392.07), Tajima's D = −3.34**; region-stratified locus-based **F_ST = 0.0568 (p = 0.003)** across Ningxia (81), Neimenggu (9), Beijing (8), others (22). SNPGenie full 13-panel done.
- Temporal gate (the paper's methodological selling point): PSTVd R² = 0.0002, DRT failed (real R² at 28.9th percentile of 20); Cytorhabdovirus genome-wide R² = 0.019 (DRT failed); only P gene passed (R² = 0.029, 100th percentile, < 3% variance). **No TMRCA reported; exploratory BEAST chains kept only as sensitivity artefacts.**
- Recombination: ViReMa/DI-tector done for **13/15** (2 failed runs — state as limitation, don't hide); RDP5 events passing ≥ 3 methods = `[N — TO CONFIRM from rdp5 output]`.
- Interpretation guardrails: dense LD = co-mutation structure of a quasispecies cloud, not classical epistasis; Tajima's D −3.34 = expansion or purifying selection signature (state both, prefer neither).

## Rules

- The temporal gate is always described as a **hard gate implemented before dating**, never as post-hoc filtering.
- Never fabricate TMRCA/HPD values; §2.8 mentions BEAST config because the platform ran it for gate-passing viruses in general, with the note that neither dominant virus qualified.
- Versions: write ⊙-marked versions only as "version to be confirmed from run logs" if not in `SOFTWARE_VERSIONS.txt` (confirmed there: iqtree 3.1.2, mafft 7.525).
