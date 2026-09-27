---
name: mmpv-paper-known-virus
description: Draft and audit the known-virus analysis sections of the MMPV SCI manuscript (Methods 2.4-2.7 and 2.9; Results 3.3-3.4, 3.6-3.7, 3.11; Tables 1-2; Figs 3-7). Use when writing/revising manuscript text about Salmon detection thresholds, high-confidence filtering, iVar/FreeBayes/LoFreq variant calling, SnpEff/SNPGenie annotation, AF-matrix reconstruction, QST haplotypes, Fisher variant-host associations, full-length assembly, SDT similarity, or the 15-virus Lycium tables.
---

# MMPV Paper — Known-Virus Analysis Sections

You own: **Methods §2.4 (detection/quantification), §2.5 (variants), §2.6 (population-genetic reconstruction), §2.7 (full-length assembly + similarity), §2.9 (statistics); Results §3.3 (detection funnel), §3.4 (virus–host associations), §3.6 (variant complexity), §3.7 (PSTVd quasispecies), §3.11 (genome recovery/similarity); Tables 1-2; Figs 3-7.**

## Authoritative sources

| Source | Take from it |
|---|---|
| `doc/METHODS_TEMPLATE.md` 管道三 | Canonical thresholds: Salmon 2.5.1; detection = coverage ≥ 10%, Poisson ratio ≥ 0.3, depth ≥ 0.5×, TPM ≥ 1.0, ≥ 10 unique reads, ANI ≥ 95%; high-confidence = coverage ≥ 50%, depth ≥ 5×; iVar 1.4.4 pipeline (`samtools mpileup -aa -A -d 0 -B -Q 0 | ivar variants -q 20 -t 0.01`); FreeBayes `-p 1 --pooled-continuous`; LoFreq 2.1.5; dynamic DP (10/20/100), AF ≥ 0.05, bcftools `QUAL>20 && DP && AF`; SnpEff 5.4c; SNPGenie; AF-matrix reconstruction rule (AF ≥ 0.05 → 1); LD r² ≥ 0.5 with same-position artefact removal; OmniVirusAssembler 12 steps, `--min_covered 10`; SDT-style identity + SciPy/Seaborn clustering; Fisher one-sided + BH; QST Hamming ≤ 0.15 union-find; Hd formula |
| `virome_analysis_pipeline/METHODS.md` + `default_profile.yaml` | Parameter fallback |
| `MMPV-RNA-paper/MMPV_Lycium_virome_manuscript_EN_v1.md` §2.4-2.7, §2.9, §3.3-3.4, §3.6-3.7, §3.11 | Established prose, Tables 1-2 (verbatim numbers) |

## Lycium result anchors (fixed; verbatim from Table 1/2 of v1)

- Funnel: 1,348 processed → 386 virus-positive → 47 taxa → **442 high-confidence associations (15 viruses, 285 samples)** → full-length 442/442 → **388** after extraction QC → recombination/DVG done for 13/15.
- Top prevalence: *Cytorhabdovirus* sp. 'lycii' 17.0% (186 HC; cov 86.2%; depth 875.8×; 5,385 sites; 210,487 raw LD pairs), PSTVd 16.4% (172 HC; cov 99.3%; depth 2,248.2×; 71 sites; 39→37 LD pairs after artefact removal, mapping to the RdRp-encoding region), TCDVd third (47 HC).
- Variant complexity split: Potyvirus sacchari 404 sites/42,664 LD; PVH 209/11,567; PSTVd/TCDVd extremely conserved.
- Reference database: **15,202 accessions** (`final.cluster.ref.fasta` / `ref_info.tsv`).
- §3.7: PCA separates PSTVd quasispecies by host (*L. barbarum* n = 127 vs *L. ruthenicum* n = 34); host-associated site count is still `[N — TO CONFIRM from variant_species_association output]`.

## Style / rules

- Past tense; every threshold stated as value, not "default".
- Never round Table 2 numbers; one-decimal percentages.
- Don't upgrade `TO CONFIRM` slots into numbers; flag any anchor conflict as `[WARN] manuscript-anchor conflict`.
- Term fixes: "predicted support ×0.1" is an internal scoring artefact — never report raw predicted-support as abundance.
