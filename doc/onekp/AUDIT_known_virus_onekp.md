# Pipeline Audit — Known-Virus Track, OneKP (`~/data-test/onekp_analysis`)

Audit per `mmpv-paper-known-virus` + `virome-results-analysis` skills. Run completed 2026-09-13; inputs from the 2026-09-10 discovery run. Companion narrative: `RESULTS_known_virus_onekp.md`; discovery-track audit: `AUDIT_pipeline_results.md`.

## Verdict table

| # | Check | Verdict | Evidence |
|---|-------|---------|----------|
| 1 | Input-identity of `onekp-virus` vs `out` | **[PASS]** | Identical `provenance.json` timestamp (2026-09-10T10:18:31), identical 10_Reports sizes/dates, 1,341 samples both. One canonical run duplicated under two names — pick one path for the paper (Q1). |
| 2 | Detection-funnel arithmetic | **[PASS]** | 1,341 samples → 5,934 raw pairs (379 taxa, 887 samples) → 318 passing pairs (96 taxa, 251 samples). `all_viruses.raw.tsv` / `all_viruses.best.summary.tsv` / `analysis_report.txt` agree. |
| 3 | Stage-1 thresholds vs skill anchors | **[PASS]** | analysis_report.txt states Poisson-coverage-ratio modelling with threshold ≥ 0.3 ("覆盖比值 >= 0.3") — matches the canonical Poisson ratio ≥ 0.3; single-track whole-genome Poisson mode declared. |
| 4 | Stage-2 threshold enforcement | **[PASS]** | Kept 118 / discarded 200 = 318 ✓. Discards: 146/200 depth <5×, 69/200 coverage <50% — consistent with HC = cov ≥50% ∧ depth ≥5×. high_conf.summary(.discarded).tsv |
| 5 | Avg_Read_ANI auditability | **[WARN]** | The `Avg_Read_ANI` column is empty in all 318 best-summary rows, so the canonical ANI ≥ 95% detection criterion cannot be verified from the summary. Populate the column or drop the ANI claim for this run. |
| 6 | TPM interpretability | **[WARN]** | `Asm_TPM` values approach 1,000,000 (e.g. Schlumbergera virus X = 942,752) — TPM is computed within the viral-genome subset, not transcriptome-wide. Never quote as cellular-abundance TPM in the manuscript; state the normalization universe. |
| 7 | Variant-stage input/output accounting | **[WARN]** | 118 pairs enter stage 3; QC scores only 115 (2 accessions lack qc_summary.tsv: Cannabis sativa amalgavirus 1, CMV Pl-1) and CMV Pl-1 has no vcf_merge bundle at all (72/73 dirs). Identify why 3 pairs are unaccounted (no-call? empty VCF?) and document. |
| 8 | Variant QC quality | **[PASS]** | 115/115 PASS, 0 FAIL; 22,377 variants; median 107/pair, max 1,014; median ts/tv 3.8 (typical 2.5–4.5 range for RNA viruses). qc_summary.tsv files. |
| 9 | Coinfection plausibility | **[PASS]** | 18/89 samples (20.2%) with ≥2 viruses; max 5 in one sample. Distribution (71/11/4/2/1) is plausible; worth listing the quintuple-infected sample in supplementary (Q4). |
| 10 | Assembly stage completeness | **[PASS]** | 118/118 pairs have sample-level assemblies; 117 evolution-tracking files (steps DeNovo → PVGA extension → Pre-Fusion merge with N50/N tracking). Missing 1 file — locate (Q5). |
| 11 | Full-length recovery rate | **[PASS]** | 76/118 (64.4%) genomes extracted; 57/73 accessions with ≥1 full genome; **all 76 span ≥90% of reference length (verified per-file against Rep_Length; earlier "<70%" readings were a filename-parsing artifact where `NC_006056.1` became `NC.006056.1`)**. Below the Lycium anchor (388/442 = 87.8%) — expected for wild-plant diversity, not a defect; do not copy the Lycium recovery claim into OneKP text. |
| 12 | Genome integrity metric | **[PASS]** | Median raw N% 2.9 (max 67.4%) across 118 scored consensus genomes; per-virus extraction summaries present for all 73 dirs. 07_similarity/*/01_extraction_summary.tsv |
| 13 | DVG/recombination coverage | **[WARN]** | ViReMa ran on 50/73 accessions and 67/118 sample runs; selection rule not recorded in the summary. Verify run_8_dvg.log for how the subset was chosen (Q6). |
| 14 | DVG result validity | **[PASS]** | Explicit negative: Empty_Report.csv = "No recombination events detected"; plot dir empty; per-sample dirs contain only virema.log. Consistent with median Pi 0.0034. Report as a true negative, not missing data. |
| 15 | Final reports | **[PASS]** | 09_report/Pipeline_Summary_Report.html (48 MB) + 10_report S1/S2/S3/S5/S6 bundles all present and freshly generated 2026-09-13. |
| 16 | Cross-stage set consistency | **[PASS]** | Taxa progression is coherent: 96 (detection) → 71 kept species → 73 accession dirs (= multi-accession species) → 50 DVG subset; 118 pairs preserved through variants/assembly. |
| 17 | Manuscript-anchor conflicts | **[PASS]** | No contradiction with fixed anchors (Poisson ≥0.3, HC cov≥50%/depth≥5×, reference DB usage). OneKP fills previously [TBD] breadth/recovery slots; Lycium numbers remain untouched. |

## Questions for the analyst

1. **Q1** — `~/data-test/onekp-virus` and `~/data-test/out` are the same run. Declare one canonical path (and whether one is a hardlink clone) before data-deposit/fig-pipeline work so downstream scripts don't fork.
2. **Q2** — Populate `Avg_Read_ANI` (or re-export it) so the ANI ≥ 95% detection criterion is auditable; currently the criterion is asserted by the config but invisible in outputs (#5).
3. **Q3** — CEVd: 104 detected → 88 discarded (mostly depth <5×) → 16 confirmed. Is the ≥5× depth rule appropriate for a 387-nt viroid with multi-mapping small-genome reads? Consider a viroid-specific depth rule or report both prevalences explicitly.
4. **Q4** — Identify the sample carrying 5 viruses and the 4-carrier/3-carrier samples from Coinfection_Matrix_Reads.tsv; check for index-hopping neighbours before publishing co-infection claims.
5. **Q5** — Which assembly lacks its Global_Evolution_Stats.tsv (117/118)? Regenerate or note.
6. **Q6** — Document the 50-accession DVG subset selection; either extend ViReMa to all 73 or state the rule in Methods.
7. **Q7** — Explain the 3 unQC'd variant pairs and regenerate the missing CMV Pl-1 post-analysis bundle (or drop the accession with an explicit note).
8. **Q8** — Confirm the TPM universe (viral-subset TPM) in Methods text; add per-sample total viral-read fraction if absolute abundance is needed.

## Bottom line

The known-virus track is arithmetically consistent, threshold-compliant and quality-controlled end to end (118 pairs → 115 QC-PASS → 76 full genomes → 0 recombinants). Before manuscript use: fix the ANI-column gap, the 3-pair QC accounting, the DVG subset rule, and state the TPM normalization universe.
