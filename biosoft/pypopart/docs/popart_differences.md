# Differences from PopART

PyPopART's network algorithms are ports of the C++ implementations in
[PopART](https://github.com/jessicawleigh/popart-current) and follow
their semantics by default. This page lists everything that differs.

## Shared machinery, matching PopART

- **Sequence condensation**: identical aligned sequences (gaps
  included) collapse into haplotypes, as in `HapNet::condenseSeqs`.
- **Site-pattern condensation**: MJN, TSW, and PN operate on condensed
  site patterns with per-column weights (`HapNet::condenseSitePats`);
  their inferred median/ancestral vertices carry condensed-space
  sequences.
- **Ambiguity handling**: distance calculations skip gaps and the full
  IUPAC ambiguity set (`N Y R M S V W K D H B`, plus `?` and `X`).
- **Weighted paths**: shortest-path lengths use genetic distances, like
  PopART's Floyd–Warshall.

## PyPopART additions (not in PopART)

- **MST** (`-a mst`): a plain minimum spanning tree; PopART only
  exposes the minimum spanning _network_.
- **Corrected distance metrics** (`-d jc|k2p|tn`): PopART's network
  algorithms use integer Hamming-style distances only. MJN/TSW/PN
  always use weighted Hamming distances regardless of `-d`, as the C++
  does.
- **Opt-in extras**, all off by default:
  - `MinimumSpanningNetwork(prune_redundant=True)` removes edges with
    equal-or-shorter alternative paths.
  - `TCS(connection_limit=...)` caps connections (an int, or
    `'auto'` to derive one from `confidence`; the CLI enables this
    via `-p/--parsimony-limit`). PopART's TCS has no connection limit
    and always yields a fully connected network.
  - `MedianJoiningNetwork(simplify=True)` smooths degree-2 medians;
    `max_median_vectors=` caps median count.

## Known divergences

- **Parsimony Network tree source**: PopART reads parsimony trees from
  a Nexus `TREES` block supplied by the user. PyPopART has no tree
  input path yet, so it generates trees by random-order stepwise
  addition with best-scoring insertion. The consensus sampling,
  ancestral-state resolution, alpha-threshold pruning, and
  connectivity guard all follow `AncestralSeqNet`.
- **Integer Neighbor-Joining (IntNJ)** is not ported.
- **AMOVA** is a simplified, non-permutational variant; PopART's
  nested/permutation AMOVA is not ported.
- **Traits**: PopART's Nexus `TRAITS`/`GEOTAGS` blocks and trait-driven
  frequency recomputation are not ported; PyPopART uses CSV metadata
  for populations instead.
- **Character masking** exists at the distance-calculation API level
  (`pairwise_distance_matrix(mask=...)`) but is not yet threaded
  through the CLI/GUI.

## Validation

`tests/validation/` holds golden networks hand-traced from the C++
algorithms (MST, MSN, TCS, MJN, TSW) plus property tests (MSN contains
every MST edge, TCS fully connects, sampled haplotypes always appear,
and more). Run them with `pytest tests/validation`.
