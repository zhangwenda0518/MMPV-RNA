"""
Statistics and population genetics analysis for PyPopART.

This package provides functions for:
- Network statistics and diversity metrics
- Population genetics measures (Tajima's D, Fu's Fs, FST, AMOVA)
- Network topology analysis
"""

from .popgen import (
    AMOVAResult,
    FstResult,
    FuFsResult,
    TajimaDResult,
    calculate_amova,
    calculate_fst_matrix,
    calculate_fu_fs,
    calculate_mismatch_distribution,
    calculate_pairwise_fst,
    calculate_tajimas_d,
)
from .statistics import (
    DiversityMetrics,
    NetworkMetrics,
    calculate_diversity_metrics,
    calculate_haplotype_frequencies,
    calculate_network_metrics,
    calculate_reticulation_index,
    calculate_summary_statistics,
    get_frequency_distribution,
    identify_central_haplotypes,
)
from .topology import (
    AncestralNode,
    Partition,
    StarPattern,
    calculate_node_centrality,
    calculate_topology_summary,
    detect_bridges,
    detect_network_partitions,
    find_central_hub_nodes,
    identify_ancestral_nodes,
    identify_bottleneck_nodes,
    identify_star_patterns,
)

__all__ = [
    # Statistics
    'calculate_haplotype_frequencies',
    'calculate_diversity_metrics',
    'calculate_network_metrics',
    'identify_central_haplotypes',
    'calculate_reticulation_index',
    'get_frequency_distribution',
    'calculate_summary_statistics',
    'DiversityMetrics',
    'NetworkMetrics',
    # Population genetics
    'calculate_tajimas_d',
    'calculate_fu_fs',
    'calculate_pairwise_fst',
    'calculate_fst_matrix',
    'calculate_amova',
    'calculate_mismatch_distribution',
    'TajimaDResult',
    'FuFsResult',
    'FstResult',
    'AMOVAResult',
    # Topology
    'identify_star_patterns',
    'detect_network_partitions',
    'calculate_node_centrality',
    'identify_ancestral_nodes',
    'calculate_topology_summary',
    'find_central_hub_nodes',
    'detect_bridges',
    'identify_bottleneck_nodes',
    'StarPattern',
    'Partition',
    'AncestralNode',
]
