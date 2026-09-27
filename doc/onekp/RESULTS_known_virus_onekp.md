# Known-Virus Analysis Results — OneKP Cohort (`~/data-test/onekp_analysis`)

**Input**: the OneKP discovery output (`~/data-test/onekp-virus`, byte-identical to `/home/zhangwenda/data-test/out`; same provenance, MMPV-RNA v2.3, 2026-09-10) and the 1,341 cleaned host-depleted libraries. **Analysis run**: stages 01–10 of the known-virus pipeline completed 2026-09-13 (`known_virus.log`; final report `09_report/Pipeline_Summary_Report.html`, 48 MB). Every figure cites its source file.

## Detection funnel (stage 1, Salmon-based quantification)

Salmon-based detection (engine "V42.4", single-track whole-genome Poisson filtering, coverage ratio ≥ 0.3) processed all 1,341 samples and retained 5,934 raw sample × virus pairs spanning 379 virus taxa in 887 samples before thresholding (01_detection/summary/all_viruses.raw.tsv; analysis_report.txt). Applying the detection thresholds left **318 sample × virus pairs from 96 virus taxa in 251 samples (18.7% of the cohort)** on the best-summary whitelist, with 0 pairs downgraded as unclassified (all_viruses.best.summary.tsv). Median replicate coverage of passing pairs was 67.7%, median depth 6.0×, median Poisson ratio 0.79 (all_viruses.best.summary.tsv). Genome-type composition of the 318 pairs: ssRNA(+) 155, ssRNA 110, dsRNA 22, ssRNA(±) 10, ssRNA(−) 9, ssDNA(±) 5, dsDNA-RT 3, other/DNA 4 (Molecule_type column).

The most prevalent taxa were Citrus exocortis Yucatan viroid (104 samples; 7.8% of the cohort), Passiflora latent virus (19), Cymbidium mosaic virus (13; mean depth 3,999×), Narcissus symptomless virus (10), Impatiens necrotic spot virus (10), Alfalfa-associated picorna-like virus 2 (8) and Schlumbergera virus X (6) (filter_stats.per_virus.tsv). Fifty-six of the 96 taxa were detected in a single sample.

## High-confidence filtering (stage 2)

The high-confidence rule (coverage ≥ 50% and depth ≥ 5×) kept **118 pairs (71 taxa, 89 samples; 6.6% of the cohort)** and discarded 200 pairs (49 taxa) (high_conf.summary.tsv; high_conf.summary.discarded.tsv). Kept pairs had median coverage 85.8% and median depth 78×, with 67/118 at ≥80% coverage. Of the 200 discards, 146 fell below 5× depth and 69 below 50% coverage. Citrus exocortis Yucatan viroid accounts for the largest single discard block (88 of its 104 detections), so its high-confidence prevalence (16 samples) is 6.5-fold lower than its detection prevalence — both figures should be quoted when contrasting cohorts.

## Variant calling and population genetics (stage 3)

 iVar-based variant calling was run on all 118 high-confidence pairs (73 reference accessions, 89 samples; all_summary.tsv). Mapping quality before calling was high: median coverage 98.3% and median depth 115.9× (03_variants/summary/all_summary.tsv). Within-sample diversity was low overall — median nucleotide diversity Pi = 0.0034 and median Shannon entropy = 0.0158.

Variant QC (per-accession `vcf_merge/stats/qc_summary.tsv`) scored 115 of the 118 pairs; **all 115 passed** (0 FAIL), yielding **22,377 variants** in total with a median of 107 variants per pair (maximum 1,014) and median transition/transversion ratio 3.8. Viruses contributing the most variants: Passiflora latent virus (2,212), Schlumbergera virus X (1,977), Cymbidium mosaic virus (1,536), Potato virus S (1,450), Garlic virus B (829), Arabis mosaic virus (823), Garlic virus X (757) and Verbena latent virus (694) (qc_summary.tsv files).

The co-infection matrix (89 samples × 71 viruses) shows 71 samples with a single virus, 11 with two, 4 with three, 2 with four and 1 sample with five concurrent viruses — 18 of 89 high-confidence-positive samples (20.2%) carry mixed infections (Coinfection_Matrix_Reads.tsv).

## Post-hoc variant analysis (stage 4)

Of 73 per-accession working directories, 72 contain the full post-analysis bundle — merged VCFs, allele-frequency matrices, phylogenetic trees, QC statistics and Figure-1A/1B variant-landscape plots — and one (Cucumber mosaic virus strain Pl-1, AM183116.1) has none (04_post_analysis/). Two accessions (Cannabis sativa amalgavirus 1; CMV Pl-1) lack a qc_summary.tsv, which explains the 115-of-118 QC coverage above.

## Read-mapping funnel and virus-count tiers

Reads mapping to the 15,202-accession known-virus reference form a three-tier funnel: 7,949,421 virus-mapped reads across 5,934 raw sample × virus pairs (379 species, 887 samples; `all_viruses.raw.tsv`) → 7,237,092 reads across 318 threshold-passing pairs (96 species, 100 accessions, 251 samples; `all_viruses.summary.tsv`, identical to the best-summary whitelist this run with 0 unclassified downgrades) → 7,225,617 reads (99.8%) retained in the 118 high-confidence pairs (71 species, 89 samples), while the 200 discarded pairs carry only 11,475 reads in total (~57 reads/pair) — i.e. the discarded detections are trace-level signals, which also explains the 104-detected-vs-16-confirmed gap for Citrus exocortis Yucatan viroid. Virus reads represent ~0.04% of the ~21.7-billion-read post-depletion pool. At the variant-calling stage, re-mapping of the same 118 pairs against individual references yields 18,995,828 reads (samtools/bowtie2 counting, a different unit from the Salmon EM estimates above).

## Full-length genome recovery (stages 5–6)

Per-sample de novo assembly with stepwise extension produced 118 sample-level assemblies across all 73 accessions (maximum 16 assemblies for Citrus exocortis Yucatan viroid), each tracked through assembly-evolution steps (DeNovo_Cleaned → PVGA_Extension → Pre_Fusion_Merge) with length, non-N bases and contig N50 recorded per step (117 Global_Evolution_Stats.tsv files; 05_assembly/). Full-length genome extraction succeeded for **76 of 118 assemblies (64.4%)**, covering 57 of 73 accessions (78.1%); **all 76 recovered genomes span ≥90% of their reference length (most at ~100%)**, with sequence integrity captured by N-content instead (median raw N 2.9%, max 67.4%; 07_similarity/*/01_extraction_summary.tsv). The deepest recoveries were Cymbidium mosaic virus (6 full genomes), followed by Citrus exocortis Yucatan viroid, Narcissus symptomless virus and Schlumbergera virus X (3 each) (06_extraction/*.full.fasta counts).

## Recombination / DVG screening (stage 8)

ViReMa-based defective-viral-genome and recombination screening ran on 50 of the 73 accessions (67 sample-level runs). The consolidated report records **"No recombination events detected"** across all runs, and the per-virus plot directory is empty (08_dvg/Summary_Analysis_Report/Empty_Report.csv; Virus_Specific_Plots/). This clean negative is consistent with the low within-sample diversity (median Pi 0.0034) of the cohort.

## Reporting (stages 9–10)

The pipeline closed with a 48 MB summary report (09_report/Pipeline_Summary_Report.html) and bundled result sets S1_Detection, S2_Filter, S3_Variants (+ coinfection matrix), S5_Assembly_Stats and S6_Post_hoc per-virus variant bundles (10_report/).

## Interpretation

The OneKP known-virus track behaves as a stringent second pass over the discovery output: Salmon detection nominates 318 sample × virus pairs (96 taxa, 18.7% of samples virus-positive), the ≥50%-coverage/≥5×-depth rule confirms 118 pairs (71 taxa, 6.6% of samples), and the confirmed set is quality-controlled end to end — 115/118 pairs pass variant QC with 22,377 variants at median ts/tv 3.8, 76 full-length genomes are recovered, and recombination screening returns a clean negative. Compared with the domesticated-cohort anchor (Lycium: 442 high-confidence pairs of 15 viruses, 87.8% full-length recovery), the wild-plant OneKP cohort shows far higher taxonomic breadth (71 vs 15 taxa) but lower per-virus depth and lower full-length recovery (64.4%), which is the expected trade-off when moving from a focused host to 1,341 wild-plant transcriptomes. The detection-vs-confirmation gap for the flagship viroid (104 detected vs 16 confirmed) is the single most manuscript-relevant number to verify before quoting prevalence.
