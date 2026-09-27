# Changelog

All notable changes to PyPopART will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed (breaking)

- **All network algorithms now follow the C++ PopART semantics by
  default**; outputs differ from earlier PyPopART releases:
  - MSN uses PopART's ascending level sweep: equal-distance ties are
    kept, and the old automatic redundant-edge removal is now the
    opt-in `prune_redundant=True`. `max_connections` was removed.
  - TCS no longer applies a connection limit by default (networks
    fully connect, as in PopART); the Poisson-derived cap is opt-in
    via `connection_limit='auto'`/an int, and the CLI's
    `-p/--parsimony-limit` only activates it when given. Inferred
    intermediates are unlabelled (no fake `N...` sequences).
  - MJN is a full MedJoinNet port: feasible-link threshold-graph
    pruning, vertex-incident triplets, and quasi-medians computed in
    condensed site-pattern space with weighted costs. The degree-2
    median smoothing is now the opt-in `simplify=True` (default off).
  - TSW is a full TightSpanWalker port (auxiliary-graph 2-colouring,
    dT >= d invariant, memoised internal vertices with no invented
    sequences).
  - Parsimony Network is a full AncestralSeqNet port with randomized
    Fitch ancestral sampling over stepwise-addition parsimony trees;
    `min_edge_frequency` was replaced by `alpha` (default 0.95, with
    a connectivity guard).
- Haplotypes deduplicate by exact aligned sequence (gaps included),
  matching PopART's condenseSeqs; gap-position-differing sequences are
  no longer merged.
- Distance calculations skip the full IUPAC ambiguity set
  (`N Y R M S V W K D H B ? X` and gaps), matching PopART.
- Distance-method aliases normalised (`tamura_nei` -> `tn`); unknown
  algorithm parameters now raise `TypeError` instead of being silently
  ignored.
- Core install slimmed: `dash`/`matplotlib`/`plotly`/`numba` moved to
  the `gui`/`viz`/`speed` extras; cartopy, folium, scikit-learn, lxml
  and pandas dropped. `scipy` remains (required by networkx layouts).
- The broken `geo-visualize` command and remaining geographic
  visualization remnants were removed.

### Added

- Top-level API: `pypopart.load_alignment`, `pypopart.build_network`,
  `pypopart.build`, and friends.
- Algorithm registry (`pypopart.algorithms.build`) shared by CLI/GUI.
- Site-pattern condensation (`pypopart.core.site_patterns`), site
  weights and character masks in `pairwise_distance_matrix`.
- Validation suite with hand-traced golden networks and property tests
  (`tests/validation/`), plus a pytest-benchmark harness
  (`tests/benchmarks/`).
- Readers parse from in-memory text (`from_string`), used by the GUI
  (no more leaked temp files); uploads detect format by content.
- GitHub Actions CI (tests on Python 3.11-3.14, lint, build + wheel
  smoke test).

### Fixed

- `pypopart analyze`/`visualize` crashed on import errors; network
  saving from the CLI failed for every format; `load_network` returned
  a raw graph downstream code could not use.
- The GUI's algorithm parameters were silently ignored due to a
  misspelled keyword; median counts always displayed as 0.
- `GenBankReader` was missing `read_alignment`.
- Large speedups: TCS ~11x (cached all-pairs paths), TSW ~10x
  (vectorised geodesic), PN ~4x (bitmask Fitch), plus a parallel
  numba kernel for Hamming distance matrices.

## [0.1.0] - Initial Release

### Added

- Core haplotype network algorithms (MST, MSN, TCS, MJN, PN)
- Distance calculation with multiple evolutionary models
- FASTA, NEXUS, PHYLIP, and GenBank file format support
- Static and interactive visualization
- Command-line interface
- Python API for programmatic access
- Network statistics and analysis tools

[Unreleased]: https://github.com/Adamtaranto/pypopart/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/Adamtaranto/pypopart/releases/tag/v0.1.0
