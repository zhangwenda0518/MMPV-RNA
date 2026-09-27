"""
Network statistics and diversity metrics for PyPopART.

This module provides functions for calculating various statistics about
haplotype networks, including diversity measures, network topology metrics,
and frequency distributions.
"""

from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import networkx as nx
import numpy as np

from ..core.alignment import Alignment
from ..core.graph import HaplotypeNetwork


@dataclass
class DiversityMetrics:
    """Container for diversity metrics."""

    haplotype_diversity: float
    nucleotide_diversity: float
    shannon_index: float
    num_haplotypes: int
    num_samples: int


@dataclass
class NetworkMetrics:
    """Container for network topology metrics."""

    diameter: int
    clustering_coefficient: float
    average_path_length: float
    reticulation_index: float
    num_components: int
    central_haplotypes: List[str]


def calculate_haplotype_frequencies(
    network: HaplotypeNetwork, normalize: bool = False
) -> Dict[str, Dict[str, float]]:
    """
    Calculate haplotype frequencies overall and by population.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.
    normalize : bool, default=False
        If True, return frequencies as proportions (0-1).

    Returns
    -------
    Dict[str, Dict[str, float]]
        Dictionary with an 'overall' entry mapping haplotype_id to
        frequency or count, and a 'by_population' entry mapping each
        population to its own haplotype_id mapping.
    """
    overall_counts: Dict[str, int] = defaultdict(int)
    pop_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))

    total_samples = 0
    pop_totals: Dict[str, int] = defaultdict(int)

    # Collect counts
    for hap_id in network.nodes:
        haplotype = network.get_haplotype(hap_id)
        if haplotype is None:
            continue

        count = haplotype.frequency
        overall_counts[hap_id] = count
        total_samples += count

        # By population
        freq_info = haplotype.get_frequency_info()
        for pop, pop_count in freq_info.by_population.items():
            pop_counts[pop][hap_id] = pop_count
            pop_totals[pop] += pop_count

    # Normalize if requested
    if normalize and total_samples > 0:
        overall_freq = {
            hap_id: count / total_samples for hap_id, count in overall_counts.items()
        }
    else:
        overall_freq = dict(overall_counts)

    by_pop_freq = {}
    for pop, counts in pop_counts.items():
        if normalize and pop_totals[pop] > 0:
            by_pop_freq[pop] = {
                hap_id: count / pop_totals[pop] for hap_id, count in counts.items()
            }
        else:
            by_pop_freq[pop] = dict(counts)

    return {
        'overall': overall_freq,
        'by_population': by_pop_freq,
        'total_samples': total_samples,
    }


def calculate_diversity_metrics(
    network: HaplotypeNetwork, alignment: Optional[Alignment] = None
) -> DiversityMetrics:
    """
    Calculate diversity metrics for a haplotype network.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.
    alignment : Alignment, optional
        Optional alignment for nucleotide diversity calculation.

    Returns
    -------
    DiversityMetrics
        DiversityMetrics object with calculated metrics.
    """
    # Get frequencies
    freq_info = calculate_haplotype_frequencies(network, normalize=True)
    overall_freq = freq_info['overall']
    total_samples = freq_info['total_samples']

    num_haplotypes = len(overall_freq)

    # Haplotype diversity (Nei's gene diversity)
    # H = (n / (n-1)) * (1 - sum(p_i^2))
    if total_samples <= 1:
        haplotype_diversity = 0.0
    else:
        sum_p_squared = sum(p**2 for p in overall_freq.values())
        haplotype_diversity = (total_samples / (total_samples - 1)) * (
            1 - sum_p_squared
        )

    # Shannon diversity index
    # H = -sum(p_i * ln(p_i))
    shannon_index = 0.0
    for p in overall_freq.values():
        if p > 0:
            shannon_index -= p * np.log(p)

    # Nucleotide diversity (π)
    nucleotide_diversity = 0.0
    if alignment is not None and total_samples > 1:
        nucleotide_diversity = _calculate_nucleotide_diversity(
            network, alignment, overall_freq
        )

    return DiversityMetrics(
        haplotype_diversity=haplotype_diversity,
        nucleotide_diversity=nucleotide_diversity,
        shannon_index=shannon_index,
        num_haplotypes=num_haplotypes,
        num_samples=total_samples,
    )


def _calculate_nucleotide_diversity(
    network: HaplotypeNetwork, alignment: Alignment, frequencies: Dict[str, float]
) -> float:
    """
    Calculate nucleotide diversity (π).

    π = sum_ij(p_i * p_j * d_ij)

    where p_i is the frequency of haplotype i and d_ij is the
    number of differences between haplotypes i and j.

    Parameters
    ----------
    network : HaplotypeNetwork
        Network holding the haplotypes.
    alignment : Alignment
        Alignment the haplotypes came from.
    frequencies : Dict[str, float]
        Relative frequency per haplotype id.

    Returns
    -------
    float
        The nucleotide diversity (π).
    """
    haplotype_ids = list(frequencies.keys())
    pi = 0.0

    for i, hap_i in enumerate(haplotype_ids):
        for j, hap_j in enumerate(haplotype_ids):
            if i >= j:
                continue

            # Get sequences
            hap_obj_i = network.get_haplotype(hap_i)
            hap_obj_j = network.get_haplotype(hap_j)

            if hap_obj_i is None or hap_obj_j is None:
                continue

            seq_i = hap_obj_i.sequence.data
            seq_j = hap_obj_j.sequence.data

            # Count differences
            differences = sum(
                1 for a, b in zip(seq_i, seq_j) if a != b and a != '-' and b != '-'
            )

            # Add to diversity
            pi += 2 * frequencies[hap_i] * frequencies[hap_j] * differences

    # Normalize by sequence length (excluding gaps)
    if haplotype_ids and alignment.length > 0:
        pi /= alignment.length

    return pi


def calculate_network_metrics(network: HaplotypeNetwork) -> NetworkMetrics:
    """
    Calculate comprehensive network topology metrics.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.

    Returns
    -------
    NetworkMetrics
        NetworkMetrics object with calculated metrics.
    """
    G = network.to_networkx()

    # Basic metrics
    if G.number_of_nodes() == 0:
        return NetworkMetrics(
            diameter=0,
            clustering_coefficient=0.0,
            average_path_length=0.0,
            reticulation_index=0.0,
            num_components=0,
            central_haplotypes=[],
        )

    # Diameter (max shortest path in largest component)
    if nx.is_connected(G):
        diameter = nx.diameter(G)
        avg_path_length = nx.average_shortest_path_length(G)
    else:
        # For disconnected graphs, use largest component
        largest_cc = max(nx.connected_components(G), key=len)
        subgraph = G.subgraph(largest_cc)
        diameter = nx.diameter(subgraph) if len(largest_cc) > 1 else 0
        avg_path_length = (
            nx.average_shortest_path_length(subgraph) if len(largest_cc) > 1 else 0.0
        )

    # Clustering coefficient
    clustering_coefficient = nx.average_clustering(G)

    # Reticulation index (proportion of alternative connections)
    # For a tree: edges = nodes - components
    # Reticulation = (actual_edges - tree_edges) / tree_edges
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()
    num_components = nx.number_connected_components(G)
    tree_edges = num_nodes - num_components

    if tree_edges > 0:
        reticulation_index = (num_edges - tree_edges) / tree_edges
    else:
        reticulation_index = 0.0

    # Identify central haplotypes (top 3 by degree centrality)
    degree_centrality = nx.degree_centrality(G)
    sorted_nodes = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)
    central_haplotypes = [node for node, _ in sorted_nodes[:3]]

    return NetworkMetrics(
        diameter=diameter,
        clustering_coefficient=clustering_coefficient,
        average_path_length=avg_path_length,
        reticulation_index=reticulation_index,
        num_components=num_components,
        central_haplotypes=central_haplotypes,
    )


def identify_central_haplotypes(
    network: HaplotypeNetwork, method: str = 'degree'
) -> List[Tuple[str, float]]:
    """
    Identify central haplotypes using various centrality measures.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.
    method : str, default='degree'
        Centrality method ('degree', 'betweenness', 'closeness', 'eigenvector').

    Returns
    -------
    List[Tuple[str, float]]
        List of (haplotype_id, centrality_score) tuples, sorted by score (descending).
    """
    G = network.to_networkx()

    if G.number_of_nodes() == 0:
        return []

    # Calculate centrality
    if method == 'degree':
        centrality = nx.degree_centrality(G)
    elif method == 'betweenness':
        centrality = nx.betweenness_centrality(G)
    elif method == 'closeness':
        if nx.is_connected(G):
            centrality = nx.closeness_centrality(G)
        else:
            centrality = {}
            for component in nx.connected_components(G):
                subgraph = G.subgraph(component)
                comp_centrality = nx.closeness_centrality(subgraph)
                centrality.update(comp_centrality)
    elif method == 'eigenvector':
        try:
            centrality = nx.eigenvector_centrality(G, max_iter=1000)
        except nx.PowerIterationFailedConvergence:
            # Fall back to degree centrality
            centrality = nx.degree_centrality(G)
    else:
        raise ValueError(f'Unknown centrality method: {method}')

    # Sort by centrality score
    sorted_centrality = sorted(centrality.items(), key=lambda x: x[1], reverse=True)

    return sorted_centrality


def calculate_reticulation_index(network: HaplotypeNetwork) -> float:
    """
    Calculate the reticulation index of the network.

    The reticulation index measures the proportion of reticulations
    (alternative connections) in the network compared to a simple tree.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.

    Returns
    -------
    float
        Reticulation index (0 for a tree, >0 for networks with reticulations).
    """
    G = network.to_networkx()

    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()
    num_components = nx.number_connected_components(G)

    if num_nodes == 0:
        return 0.0

    # For a tree: edges = nodes - components
    tree_edges = num_nodes - num_components

    if tree_edges == 0:
        return 0.0

    # Reticulation index
    return (num_edges - tree_edges) / tree_edges


def get_frequency_distribution(network: HaplotypeNetwork) -> Dict[int, int]:
    """
    Get the frequency distribution of haplotypes.

    Returns a dictionary mapping frequency values to the number of
    haplotypes with that frequency.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.

    Returns
    -------
    Dict[int, int]
        Dictionary mapping frequency -> count of haplotypes.
    """
    frequencies = []

    for hap_id in network.nodes:
        haplotype = network.get_haplotype(hap_id)
        if haplotype is not None:
            frequencies.append(haplotype.frequency)

    return dict(Counter(frequencies))


def calculate_summary_statistics(
    network: HaplotypeNetwork, alignment: Optional[Alignment] = None
) -> Dict[str, any]:
    """
    Calculate comprehensive summary statistics for a network.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.
    alignment : Alignment, optional
        Optional alignment for additional metrics.

    Returns
    -------
    Dict[str, any]
        Dictionary with all calculated statistics.
    """
    # Get diversity metrics
    diversity = calculate_diversity_metrics(network, alignment)

    # Get network metrics
    net_metrics = calculate_network_metrics(network)

    # Get frequency information
    freq_info = calculate_haplotype_frequencies(network, normalize=False)

    # Get frequency distribution
    freq_dist = get_frequency_distribution(network)

    return {
        'diversity': {
            'haplotype_diversity': diversity.haplotype_diversity,
            'nucleotide_diversity': diversity.nucleotide_diversity,
            'shannon_index': diversity.shannon_index,
            'num_haplotypes': diversity.num_haplotypes,
            'num_samples': diversity.num_samples,
        },
        'network': {
            'diameter': net_metrics.diameter,
            'clustering_coefficient': net_metrics.clustering_coefficient,
            'average_path_length': net_metrics.average_path_length,
            'reticulation_index': net_metrics.reticulation_index,
            'num_components': net_metrics.num_components,
            'central_haplotypes': net_metrics.central_haplotypes,
        },
        'frequencies': freq_info,
        'frequency_distribution': freq_dist,
    }
