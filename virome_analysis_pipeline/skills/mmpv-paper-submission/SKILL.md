---
name: mmpv-paper-submission
description: Draft and audit the submission/data-availability sections of the MMPV SCI manuscript (Methods submission paragraph; Declarations — data availability, code availability, author contributions; Table S1; the .sqn/GenBank workflow description). Use when writing/revising manuscript text about suvtk taxonomy/features/table2asn, MIUVIG compliance, GenBank deposition plans, or the platform citation and provenance statements.
---

# MMPV Paper — Submission & Declarations Sections

You own: **Methods submission paragraph (inside §2.3 or as §2.10 per journal style), Declarations (data/code availability, author contributions skeleton), Table S1 caption, the platform-citation statement.**

## Authoritative sources

| Source | Take from it |
|---|---|
| `doc/METHODS_TEMPLATE.md` 管道五 | Canonical methods: suvtk v0.1.1 five steps — (i) `suvtk taxonomy` MMseqs2 LCA + genome-type prediction (`-s 0.7`); (ii) `suvtk features` direction correction + pyrodigal-gv ORFs (coding-complete, CDS > 50% genome) + BFVD annotation, five-column .tbl, hypothetical-protein labelling; (iii) `lcl|` ID normalization; (iv) source.src / miuvig.tsv / assembly.tsv / template.sbt per GenBank+MIUVIG; (v) `suvtk comments` + `suvtk table2asn` → submission.sqn + .val; pre-submission checklist (taxonomy non-NA, no internal stops, placeholders replaced, segmented-isolate consistency, metagenomic = TRUE, no ERROR in .val) |
| `virome_submission_pipeline/METHODS.md` | Detail fallback |
| `doc/MMPV_dependencies.md` D3 | suvtk DB path/deployment (`suvtk download-database`, ~5 GB) |
| `MMPV-RNA-paper/MMPV_Lycium_virome_manuscript_EN_v1.md` Declarations | Existing declarations skeleton |

## Fixed conventions

- Platform citation: "MMPV (Massive Meta-mining of Plant Viruses)" with GitHub/DOI placeholder `[TO CONFIRM]`; code version **v2.3** (matches `virome_pipeline.py`/`pixi.toml`) while README displays v3.0 as the platform name — in the manuscript use "MMPV v3.x" only if release alignment is confirmed, otherwise cite the repo without version or use v2.3 with a note. Flag the mismatch, don't paper over it.
- Data availability must name: NCBI SRA + CNCB GSA source accessions (Table S2), GenBank deposition `[accessions — TO CONFIRM]`, output archives `known_virus_all_v2/` and `virome_phylo_pipeline/phylo_results/` (PSTVd `pstvd_full_20260831`, Cytorhabdovirus `gcva_full_20260831`).
- Novel genomes enter GenBank **through this pipeline's own output** (.sqn from rescue-stage centroids and analysis-stage full genomes) — keep that loop explicit; it is a platform selling point.
- Software without papers get GitHub URL + access date (suvtk, CAPHEINE, VirNA, RdRp-Catch, ViralVerify, VITAP, ACVirus, RNAVirHost, vConTACT3→cite vConTACT2 paper, BBTools, BWA-MEM2, CoverM, Penguin, PhaBOX2, CHERRY, RVDB); NCBI tools cite URL.

## Anti-fabrication rules

- Accession numbers, author lists, funding, e-mail: always bracketed placeholders; never keep a template example e-mail.
- Table S1 = "Software and versions" must match `SOFTWARE_VERSIONS.txt` where versions are confirmed (samtools 1.21, diamond 2.2.0, iqtree 3.1.2, mafft 7.525, salmon 2.5.1, megahit 1.2.9, mmseqs bd01c22); others stay ⊙/TO CONFIRM.
