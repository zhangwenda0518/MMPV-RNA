# Pipeline Audit — OneKP Run `/home/zhangwenda/data-test/out`

Audit per `virome-results-analysis` skill (MMPV-RNA v2.3 run, reports dated 2026-09-12, provenance 2026-09-10). Companion narrative: `RESULTS_virus_discovery.md`.

## Verdict table

| # | Check | Verdict | Evidence |
|---|-------|---------|----------|
| 1 | Reads-funnel monotonicity (global) | **[PASS]** | 66.13 B reads → 64.67 B (fastp, 97.8%) → 12.52 B pairs (dedup unit change, see #14) → 10.86 B pairs after Kraken2 (86.8%). data_summary.tsv / hostdep_summary.tsv |
| 2 | Identification-funnel monotonicity (global) | **[PASS]** | 6,967,670 candidates → 1,476,114 filtered (21.2%) → 761,295 clustered → 610,535 vOTUs → 610,425 centroids. ident/filter_summary.tsv, cluster_pipeline_reduction.tsv |
| 3 | Identification-funnel monotonicity (per-sample) | **[WARN]** | ≥1 sample violates: ERR2040117 has `All_Candidate=0` in ident_summary.tsv yet `Passed=466` in filter_summary.tsv, producing `Retained(%)=46600.0`. Exactly 1 bad-row flagged in each filter mode. Fix the ÷0 guard and re-sync the two tables. |
| 4 | Filter-mode differentiation | **[WARN]** | `comb`, `strict`, `filter` modes all pass exactly 1,476,114 — the mode switch is not biting. Confirm which mode the manuscript numbers should quote. filter_summary.tsv |
| 5 | Tool panel sanity | **[WARN]** | VirSorter2 / VirHunter / VirBot / ViralML are zero in all 1,343 ident rows; BLAST/Metabuli/ViralVerify/geNomad/RdRp-Catch are active. `provenance.json` tool registry lists neither group (no metabuli/rdrpcatch, no vs2/virhunter…), so the registry is not authoritative. Confirm disabled-vs-broken before claiming a "ten-tool panel" for this run. |
| 6 | COBRA accounting | **[WARN]** | 115,180 queries, `Extended_Circular=0`, `Extended_Partial=0`, `Extended_Failed=1`, yet `Total_Gain(bp)=274,009,544` and median orphan rate 0.0%. Gain with zero extensions is contradictory — column semantics need re-derivation. cobra_summary.tsv |
| 7 | Cluster structure | **[PASS]** | 88.8% singletons, max cluster 64, mean 1.25 — extreme but plausible for transcriptome data; reduction 19.8% recorded consistently in 3 files. cluster_size_distribution.tsv, global_summary.txt |
| 8 | Reference-guided (CD-HIT) propagation | **[WARN]** | 448 contigs link to known refs (known_association.tsv) but `ref_cov` is 0.00 for all 20,054 rows of plant_virus_cluster_summary.tsv — ref-guided coverage was never propagated to the plant cluster table. |
| 9 | Cross-source taxonomy conflicts (R vs suvtk) | **[N/A (stage output absent)]** | No `suvtk_taxonomy/` under 09_Virome_Analysis; no structural-annotation cross-check table. Per skill rule the number is not inferred. |
| 10 | Taxonomy fill / agreement | **[PASS]** | Fill declines Realm 91.8% → Family 81.3% → Genus 56.1%; species agreement mmseqs 0.940 / VITAP 0.217 — VITAP weak-on-plant-viruses matches the manuscript anchor. taxonomy_fill_stats.tsv, taxonomy_agreement_stats.tsv |
| 11 | Host-table consistency | **[WARN]** | Two host universes coexist: host_ictv_classification_summary.tsv (Sep 04): Plant 19,210 / total 535,103; host_distribution.tsv + checkv_summary.tsv (Sep 12): Plant 20,054 / total 534,438. taxonomy_fill_stats denominator = 535,103. Stale-file drift (≈665 rows) — regenerate or annotate which is canonical. |
| 12 | CheckV distribution | **[PASS]** | 534,438 assessed: Complete 364 / HQ 1,237 / Medium 8,473 / Low 257,049 / ND 267,315; Plant: 172 / 447 / 1,308 / 2,286 / 15,841. Internally sums to totals. checkv_summary.tsv |
| 13 | Rescue-cascade bookkeeping | **[WARN]** | rescue_summary.md: inputs 20,517 → branch products 9,392 (A 130, B 50, C 9,212, D 0) → final 9,215; the file itself admits the product-vs-centroid 口径 mismatch. plant_virus_summary.tsv `branch` column contradicts it ('fail' 8,669, 'B' 48, '' 198) — branch attribution of final vOTUs is not reconstructable from that column. |
| 14 | Host-depletion units & no-op aligner | **[WARN]** | hostdep input (12.52 B pairs, post-Clumpify) ≠ fastp clean reads (64.67 B reads) — units differ and must be stated in Methods. `After_Kraken2 == After_Host` exactly: the Bowtie2 host-alignment step removed 0 reads — expected for OneKP (no per-sample host genome) but verify the config intentionally skipped it. hostdep_summary.tsv |
| 15 | rRNA step | **[WARN]** | `Total_rRNA=0` for all samples. Matches shipped default (RiboDetector/rRNA step off), but the empty column invites misreading as "0 rRNA found". hostdep_summary.tsv |
| 16 | Length integrity (≥1,000 bp rule) | **[WARN]** | 997 of 8,915 (11.2%) scored plant-virus records are <1,000 bp (min 500) in plant_virus_summary.tsv. Likely rescue products admitted below the candidate threshold — either enforce the rule post-rescue or document the exception. class_KEEP.fasta min length not yet re-measured (Q5). |
| 17 | Evidence-chain vs final-judgement consistency | **[FAIL]** | In rescue_evidence_scored.tsv, 1,336 rows carry `verdict=KEEP` while their `category` says Host gene contamination (963), Insect virus/Contamination (357) or Fungal/Contamination (16). KEEP-verdict and final category disagree for ~half of the KEEP set; class_KEEP.fasta holds 2,322 seqs vs 2,592 KEEP verdicts. This is the highest-priority defect before quoting KEEP counts. |
| 18 | final_judgement_table host_species column | **[FAIL]** | Constant value `out` in all 9,414 rows — the host-species column was never populated. Plant-species-level virus attribution currently impossible from this table. |
| 19 | CDD pass-rate sanity | **[PASS]** | PASS tiers 2,778/8,915 (31.2%), AMBIGUOUS 42.6%, REVIEW 26.2% — filter is biting (neither saturated nor empty). cdd_evidence_report.tsv |
| 20 | Prevalence plausibility | **[PASS]** | Only 1 zero-prevalence row; median 11 samples; 52.5% of records in ≥10 samples; max 1,296/1,341 samples (near-universal — likely a plant endophytic/shared virus such as the CEVd/viroid complex; worth a manual look, see Q6). frequency_table.tsv |
| 21 | Viroid metrics | **[WARN]** | `coverage_pct` aggregates over contigs (up to 9,349%) — not per-genome coverage; mean aligned length 100–300 nt vs 342–410 nt refs shows mostly partial hits. Re-report per-contig coverage and circularization support before claiming full viroid genomes. Viroid.per_virus.tsv |
| 22 | Threshold sensitivity (ANI 0.95/qcov 0.85 vs stricter) | **[SKIP at this scale]** | vclust raw output exists under 04_CLUSTER/3_vclust; re-derivation over 761,295 inputs deferred. |
| 23 | Manuscript-anchor conflicts | **[WARN]** | (a) Ten-tool panel claim vs 4 silent tools (#5); (b) host-depletion anchor says Kraken2+bowtie2 two-step — bowtie2 leg is a no-op here (#14); (c) rRNA optional-off as anchored (#15); (d) taxonomy ensemble agreement values consistent with anchors (VITAP weakest). Resolve (a)/(b) wording before submission. |

## Questions for the analyst

1. **Q1 (blocks manuscript numbers)** — KEEP-verdict rows categorized as contamination (#17): should the 963 "Host gene contamination + KEEP" rows be downgraded to REVIEW/DROP, or does the contamination category only annotate the taxonomy tool's call while the 5-layer evidence stands? Decide, then regenerate `class_KEEP.fasta` (2,322 vs 2,592 KEEP).
2. **Q2** — Which filter mode (`comb`/`strict`/`filter`, currently identical) is the manuscript's "two-layer filtering"? (#4)
3. **Q3** — VirSorter2/VirHunter/VirBot/ViralML: intentionally disabled for OneKP, or failed silently? Check the 02a logs; if disabled, the paper must say "5 of 9 tools produced calls in this run". (#5)
4. **Q4** — Split ★★ vs ★★★ among the 461 blast-protein + 2,095 domain-only novel candidates using the 70% protein-identity threshold (`aa_pident` in rescue_evidence_scored.tsv) so Results can quote exact tier counts.
5. **Q5** — Measure min/median length of class_KEEP.fasta and confirm ≥1,000 bp, and decide the fate of the 997 sub-1 kb rescue records (#16).
6. **Q6** — Inspect the 1,296-sample vOTU (frequency_table.tsv): genuine ultra-prevalent plant virus/viroid vs index-hopping or reference-reagent carryover.
7. **Q7** — Regenerate host tables from one snapshot (#11) and populate `host_species` in final_judgement_table.tsv from OneKP sample metadata (the mapping to plant species currently exists only upstream in onekp_analysis).
8. **Q8** — COBRA: re-derive what `Extended_*` columns mean vs `Total_Gain(bp)` (#6); if 274 Mbp of flanks are real, state per-contig median gain in Methods.
9. **Q9** — Suvtk structural-annotation stage absent → run it or drop the "R-vs-suvtk cross-check" claim from the methods-to-results chain (#9).
10. **Q10** — Confirm the plant universe canonical count (20,054 host-assigned vs 20,517 rescue inputs vs 9,215 post-rescue vs 8,915 scored vs 9,414 judged) and publish a one-line accounting chain in the paper's supplementary methods.

## Bottom line

The global funnel, taxonomy, host prediction, prevalence and CDD tiering are consistent and quotable. Two [FAIL] bookkeeping defects (#17, #18) and the disabled-tool/filter-mode ambiguities (#4, #5) must be resolved before the KEEP counts, the "ten-tool panel" claim, or any per-species host attribution enter the manuscript.
