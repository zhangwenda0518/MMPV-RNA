# Virus Discovery Results — OneKP Cohort (MMPV-RNA v2.3)

**Run**: `/home/zhangwenda/data-test/out` · pipeline MMPV-RNA v2.3 · provenance timestamp 2026-09-10 · ICTV MSL41 reference database (`provenance.json`). Reports regenerated 2026-09-12 (`10_Reports/` file dates). Cohort: the OneKP (1,000 Plants) transcriptome dataset — 1,341 plant samples plus one co-assembly (1,342 rows in `data_summary.tsv`; 1,341 rows in `hostdep_summary.tsv` / per-sample assembly table).

All counts are derived from the staged output files; source file for every figure is given in parentheses. Percentages to one decimal, n = counts. N/A marks a stage output that is absent — no number was inferred.

## Sequence preprocessing and host depletion

Preprocessing retained 64,672,764,260 of 66,126,541,256 raw reads (97.8%; median per-sample clean Q20 = 97.7%, clean Q30 = 92.0%, duplicate rate = 20.2%) from 5,826.0 Gbp of raw sequence (data_summary.tsv). Host depletion (Kraken2 confidence 0.2 + Bowtie2 two-step) was run on 12,517,071,378 read pairs (post-Clumpify pairs, a different unit from the fastp read counts above); Kraken2 removed 13.2% of pairs, leaving 10,864,586,837 pairs, and the host-genome alignment step removed no additional pairs (hostdep_summary.tsv). The rRNA-depletion accounting columns report zero across all samples (hostdep_summary.tsv), consistent with the shipped default in which the optional rRNA step is disabled.

## Assembly

rnaviralSPAdes assembled each of the 1,341 samples individually, with one additional co-assembly row ("all") in the assembly table (assembly_summary.tsv). The cohort produced 179.95 M contigs totalling 92.24 Gbp, with per-sample median N50 = 853 bp (max 3,286 bp) and a global maximum contig of 52,036 bp; 11.34 M contigs (6.3%) exceeded 1,000 bp, and the median per-sample fraction above 1,000 bp was 13.8% (assembly_summary.tsv). These short contigs are typical of viral reconstruction from transcriptomes and motivated the rescue cascade described below.

## Viral identification funnel

Summed over samples, 6,967,670 contigs were flagged by at least one identification tool (ident_summary.tsv). Tool-level positive sums were: DIAMOND BLASTX 4,726,082; Metabuli 2,911,776; ViralVerify 599,938; geNomad 378,486; RdRp-Catch 55,838. VirSorter2, VirHunter, VirBot and ViralML recorded zero calls in every sample (ident_summary.tsv; see AUDIT — tools either disabled or silent in this run). Two-layer false-positive filtering passed 1,476,114 candidates (21.2% of tool-flagged contigs; filter_summary.tsv). COBRA extension processed 115,180 queries and appended 274.0 Mbp of sequence in total, with zero fully circular extensions recorded (cobra_summary.tsv; accounting semantics flagged in AUDIT).

## Clustering and the vOTU catalog

Clustering of 761,295 filtered inputs yielded 610,535 vOTUs reduced to 610,425 final centroids, a 19.8% redundancy reduction (cluster_pipeline_reduction.tsv; vclust global_summary.txt). The catalog is dominated by singletons: 542,349 vOTUs (88.8%) contain a single contig, 68,186 clusters have ≥2 members, and the largest cluster holds 64 contigs with a mean cluster size of 1.25 (cluster_size_distribution.tsv; global_summary.txt). The reference-guided CD-HIT step associated 448 contigs with known reference genomes (04_CLUSTER/2_cdhit/known_association.tsv). For the plant-hosted subset, 20,054 cluster records comprise 18,377 singletons (91.6%) and 1,677 multi-member clusters, the largest containing 53 members (plant_virus_cluster_summary.tsv).

## Taxonomic composition

Across the 534,438 host-assigned vOTUs, classification is dominated by Nucleocytoviria (290,397 vOTUs; chiefly Nucleocytoviricota, 271,083; Mimiviridae, 120,124; Chlorovirus, 17,267), followed by Duplodnaviria (101,806) and Riboviria (98,142) (taxonomy_composition.tsv; host_distribution.tsv). This composition mirrors the algal/protist biomass that dominates OneKP transcriptomes. Rank fill rates decline from Realm (91.8%) to Genus (56.1%) and Species (66.2%) (taxonomy_fill_stats.tsv); species-level cross-tool agreement ranged from 0.940 (MMseqs2) and 0.863 (DIAMOND LCA) down to 0.545 (geNomad) and 0.217 (VITAP) (taxonomy_agreement_stats.tsv).

Host prediction assigned 81.7% of vOTUs by the ICTV-preferred rule (436,801 of 534,438), with 13,671 (2.6%) left unassigned; Plant-hosted vOTUs number 20,054 (3.75%) (host_decision_method.tsv; host_distribution.tsv).

Within the plant subset (n = 20,054; plant_final_taxonomy.tsv), 75 virus families were detected, led by Partitiviridae (2,608), Caulimoviridae (2,117), Potyviridae (2,033), Tombusviridae (1,642), Tymoviridae (1,000), Secoviridae (837), Rhabdoviridae (803), Endornaviridae (655) and Closteroviridae (631); the leading genera were Potyvirus (1,712), Deltapartitivirus (1,608), Badnavirus (1,061), Luteovirus (943) and Maculavirus (643). Species-level labels were filled for 16,538 contigs (82.5%). The ACVirus module independently classified 1,969 nucleotide records with per-rank confidence (09b_Analysis_Verify/acvirus_classify/final_result_with_confidence.tsv).

## Completeness, rescue cascade and HQ vOTUs

CheckV assessment of all host-assigned vOTUs classified 364 as complete, 1,237 as high-quality, 8,473 as medium-quality, 257,049 as low-quality and 267,315 as not-determined (checkv_summary.tsv); the plant subset held 172 complete, 447 high-quality and 1,308 medium-quality genomes (619 HQ-col total) with 15,841 not-determined (checkv_summary.tsv, row Plant). Within the scored plant-virus evidence table, 3,901 contigs carried a completeness estimate (median 8.8%; 185 contigs ≥90%) and AAI confidence was high for 979, medium for 1,430, low for 1,492 and unavailable for 5,014 (plant_virus_summary.tsv; rescue_evidence_scored.tsv).

The four-branch rescue cascade started from 20,517 plant vOTU centroids: branch A (CheckV ≥90% protein completeness) rescued 130, branch B (Virseqimprover read extension) 50, branch C (dc-megablast + RGA/ragtag genus-guided scaffold extension) 9,212, and branch D (genus-length fallback) 0 — 9,392 branch products in total, condensing to 9,215 final non-redundant plant vOTUs; 20,467 inputs remained unrescued, chiefly because they were too short (<2,000 bp; 59.9% of unsaved) (08_Rescue/Plant/rescue_summary.md). Accordingly, the "rescued" source label covers 8,717 of the 8,915 scored plant-virus records versus 198 flagged as needing no rescue (plant_virus_summary.tsv, source column).

## Five-layer evidence chain and novelty verdicts

Five-layer evidence scoring (0.30·BLASTN + 0.30·BLASTX + 0.40·CDD) of the 8,915 plant-virus records returned 2,592 KEEP (29.1%), 2,614 REVIEW (29.3%) and 3,709 DROP (41.6%), with a median total score of 28.0/100 (rescue_evidence_scored.tsv). CDD domain tiers were PASS_CORE 590, PASS_CORE_FAM 139, PASS_VIRAL 2,049, AMBIGUOUS 3,800 and REVIEW 2,337, i.e. 31.2% of records carry a passing CDD tier (cdd_evidence_report.tsv). Novelty labels were: known virus by BLAST protein 222; novel candidate by BLAST protein 461; novel candidate with domain-only evidence 2,095; novel candidate with BLAST but no viral domain 1,291; short fragment pending 1,137; host/non-viral 3,709 (rescue_evidence_scored.tsv, novelty column). Against the three-tier scheme these correspond to ★ known (≥90% nt identity and coverage) for the 222 blast-confirmed knowns, with ★★/★★★ separation of the remaining candidates pending the protein-identity threshold split (see AUDIT, question Q4).

The final judgement table (n = 9,414, which includes the scored records plus viroid/HQ additions) categorised: Known plant virus 286; novel-species candidates (新种候选) 184; distantly-related virus candidate 843; distantly-related undetermined 2,173; undetermined 1,396; host-gene contamination 3,793 (40.3%); insect-virus contamination 586; fungal-virus contamination 153 (final_judgement_table.tsv). The curated keep set class_KEEP.fasta contains 2,322 sequences (09b_Analysis_Verify/class_KEEP.fasta).

## Prevalence

Median sample prevalence of the 8,915 plant-virus records was 11 samples; 7,217 (80.9%) were detected in ≥2 samples, 4,679 (52.5%) in ≥10 samples, and 1,698 (19.0%) were single-sample; one record reached 1,296 samples (~97% of the cohort) (frequency_table.tsv). A single zero-prevalence row was found (frequency_table.tsv).

## Viroids

Viroid screening detected 4 species across 22 reference accessions, 573 contigs and 404 samples (30.1% of the cohort; Viroid.per_virus.tsv; Viroid.per_sample.tsv). Citrus exocortis viroid (Pospiviroid exocortiscitri) dominated with 563 contigs in 398 samples at 94.3% mean identity; Pospiviroid latenscolumneae (7 contigs, 4 samples), P. alphairesinis (2 contigs) and Coleviroid zetacolei (1 contig) were rare (Viroid.species_info.tsv). Per-contig mean aligned lengths of 100–300 nt against 342–410 nt references indicate partial-genome hits predominate, and the reported coverage_pct column aggregates across contigs (values up to 9,349%) and should not be read as per-genome coverage (Viroid.per_virus.tsv; flagged in AUDIT).

## OneKP high-confidence call set

The dedicated OneKP detection/filtering layer consolidated assembly+read evidence into a per-sample × per-virus high-confidence table: 118 sample-virus pairs retained and 200 discarded, covering 96 unique viruses; the most prevalent call is Citrus exocortis Yucatan viroid in 104 samples (onekp_analysis/02_filtering/high_conf.summary.tsv; filter_stats.per_virus.tsv). Read-based rescue detection produced a per-sample species matrix (rescue_detection_summary.tsv).

## Interpretation

The OneKP transcriptome cohort behaves as a plant-associated meta-virome in which the discovery signal is carried by RNA viruses of plant, fungal and insect origin plus an overwhelming algal giant-virus background. The funnel is internally coherent at the global level (6.97 M tool-flagged contigs → 1.48 M filtered → 761,295 clustered → 610,535 vOTUs), and the 20,054 plant-hosted vOTUs condense through the rescue cascade into 9,215 workable genomes from which the evidence chain retains 2,322 sequences, 286 known plant viruses and 184 novel-species candidates. Low completeness (median 8.8%) and 88.8% singleton clustering are expected consequences of transcriptomic, fragmentary viral DNA/RNA rather than pipeline faults, but several bookkeeping defects documented in the companion audit (KEEP-vs-contamination verdict conflicts, a constant host_species column, zero-hit tools, a 1,000 bp length-rule exception list) must be resolved before these counts are quoted in the manuscript.
