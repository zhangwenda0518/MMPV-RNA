"""
Network topology analysis for PyPopART.

This module provides functions for analyzing the topology of haplotype
networks, including identifying star-like patterns, partitions, and
ancestral nodes.
"""

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import networkx as nx

from ..core.graph import HaplotypeNetwork


@dataclass
class StarPattern:
    """Represents a star-like pattern in the network."""

    center: str
    leaves: List[str]
    center_frequency: int
    total_frequency: int
    is_perfect_star: bool


@dataclass
class Partition:
    """Represents a partition (connected component) in the network."""

    nodes: List[str]
    num_nodes: int
    total_samples: int
    diameter: int
    is_connected: bool


@dataclass
class AncestralNode:
    """Represents a potential ancestral node."""

    node_id: str
    score: float
    degree: int
    frequency: int
    centrality: float
    is_median_vector: bool


def identify_star_patterns(
    network: HaplotypeNetwork, min_leaves: int = 3
) -> List[StarPattern]:
    """
    Identify star-like patterns in the network.

    A star pattern has a central node connected to multiple leaves
    (nodes with degree 1).

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.
    min_leaves : int, default=3
        Minimum number of leaves for a pattern to be considered.

    Returns
    -------
    List[StarPattern]
        List of StarPattern objects.
    """
    G = network.to_networkx()
    patterns = []

    for node in G.nodes():
        neighbors = list(G.neighbors(node))

        # Check if this node could be a star center
        if len(neighbors) < min_leaves:
            continue

        # Count how many neighbors are leaves (degree 1)
        leaves = [n for n in neighbors if G.degree(n) == 1]

        if len(leaves) >= min_leaves:
            # Get frequencies
            center_hap = network.get_haplotype(node)
            center_freq = center_hap.frequency if center_hap else 0

            total_freq = center_freq
            for leaf in leaves:
                leaf_hap = network.get_haplotype(leaf)
                if leaf_hap:
                    total_freq += leaf_hap.frequency

            # Check if it's a perfect star (all neighbors are leaves)
            is_perfect = len(leaves) == len(neighbors)

            pattern = StarPattern(
                center=node,
                leaves=leaves,
                center_frequency=center_freq,
                total_frequency=total_freq,
                is_perfect_star=is_perfect,
            )
            patterns.append(pattern)

    # Sort by number of leaves (descending)
    patterns.sort(key=lambda p: len(p.leaves), reverse=True)

    return patterns


def detect_network_partitions(network: HaplotypeNetwork) -> List[Partition]:
    """
    Detect partitions (connected components) in the network.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.

    Returns
    -------
    List[Partition]
        List of Partition objects, sorted by size (descending).
    """
    G = network.to_networkx()
    components = list(nx.connected_components(G))

    partitions = []
    for component in components:
        nodes = list(component)
        num_nodes = len(nodes)

        # Calculate total samples
        total_samples = 0
        for node in nodes:
            hap = network.get_haplotype(node)
            if hap:
                total_samples += hap.frequency

        # Calculate diameter
        subgraph = G.subgraph(component)
        if num_nodes > 1:
            diameter = nx.diameter(subgraph)
        else:
            diameter = 0

        partition = Partition(
            nodes=nodes,
            num_nodes=num_nodes,
            total_samples=total_samples,
            diameter=diameter,
            is_connected=True,
        )
        partitions.append(partition)

    # Sort by number of nodes (descending)
    partitions.sort(key=lambda p: p.num_nodes, reverse=True)

    return partitions


def calculate_node_centrality(
    network: HaplotypeNetwork, methods: Optional[List[str]] = None
) -> Dict[str, Dict[str, float]]:
    """
    Calculate various centrality measures for all nodes.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.
    methods : List[str], optional
        Centrality methods to calculate: 'degree', 'betweenness',
        'closeness', 'eigenvector'. If None, all methods are used.

    Returns
    -------
    Dict[str, Dict[str, float]]
        Dictionary mapping node_id -> {method -> centrality_score}.
    """
    if methods is None:
        methods = ['degree', 'betweenness', 'closeness', 'eigenvector']

    G = network.to_networkx()
    centrality_scores: Dict[str, Dict[str, float]] = defaultdict(dict)

    if G.number_of_nodes() == 0:
        return {}

    # Calculate each centrality measure
    for method in methods:
        if method == 'degree':
            scores = nx.degree_centrality(G)
        elif method == 'betweenness':
            scores = nx.betweenness_centrality(G)
        elif method == 'closeness':
            if nx.is_connected(G):
                scores = nx.closeness_centrality(G)
            else:
                # Calculate per component
                scores = {}
                for component in nx.connected_components(G):
                    subgraph = G.subgraph(component)
                    comp_scores = nx.closeness_centrality(subgraph)
                    scores.update(comp_scores)
        elif method == 'eigenvector':
            try:
                scores = nx.eigenvector_centrality(G, max_iter=1000)
            except (nx.PowerIterationFailedConvergence, nx.NetworkXError):
                # Fall back to degree centrality
                scores = nx.degree_centrality(G)
        else:
            continue

        for node, score in scores.items():
            centrality_scores[node][method] = score

    return dict(centrality_scores)


def identify_ancestral_nodes(
    network: HaplotypeNetwork, top_n: int = 5
) -> List[AncestralNode]:
    """
    Identify potential ancestral nodes in the network.

    Ancestral nodes are typically characterized by:
    - High frequency
    - High degree (many connections)
    - Central position in the network
    - High betweenness centrality

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.
    top_n : int, default=5
        Number of top candidates to return.

    Returns
    -------
    List[AncestralNode]
        List of AncestralNode objects, sorted by score (descending).
    """
    G = network.to_networkx()

    if G.number_of_nodes() == 0:
        return []

    # Calculate centrality measures
    nx.degree_centrality(G)
    betweenness_cent = nx.betweenness_centrality(G)

    # Score each node
    candidates = []

    for node in G.nodes():
        haplotype = network.get_haplotype(node)
        if haplotype is None:
            continue

        frequency = haplotype.frequency
        degree = G.degree(node)

        # Calculate composite score
        # Normalize each component
        max_freq = max(
            network.get_haplotype(n).frequency
            for n in G.nodes()
            if network.get_haplotype(n) is not None
        )
        max_degree = max(G.degree(n) for n in G.nodes())

        freq_score = frequency / max_freq if max_freq > 0 else 0
        degree_score = degree / max_degree if max_degree > 0 else 0
        betweenness_score = betweenness_cent[node]

        # Composite score (weighted average)
        score = 0.4 * freq_score + 0.3 * degree_score + 0.3 * betweenness_score

        is_median = node in network._median_vectors

        ancestral = AncestralNode(
            node_id=node,
            score=score,
            degree=degree,
            frequency=frequency,
            centrality=betweenness_cent[node],
            is_median_vector=is_median,
        )
        candidates.append(ancestral)

    # Sort by score (descending) and return top N
    candidates.sort(key=lambda x: x.score, reverse=True)

    return candidates[:top_n]


def calculate_topology_summary(network: HaplotypeNetwork) -> Dict[str, any]:
    """
    Create a comprehensive topology summary report.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.

    Returns
    -------
    Dict[str, any]
        Dictionary with topology analysis results.
    """
    G = network.to_networkx()

    # Basic topology
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()

    if num_nodes == 0:
        return {
            'num_nodes': 0,
            'num_edges': 0,
            'star_patterns': [],
            'partitions': [],
            'ancestral_candidates': [],
            'centrality': {},
        }

    # Identify patterns
    star_patterns = identify_star_patterns(network, min_leaves=2)
    partitions = detect_network_partitions(network)
    ancestral_candidates = identify_ancestral_nodes(network, top_n=5)

    # Calculate centrality
    centrality = calculate_node_centrality(network)

    # Network connectivity
    is_connected = nx.is_connected(G)
    num_components = nx.number_connected_components(G)

    # Degree distribution
    degrees = [G.degree(n) for n in G.nodes()]
    degree_dist = {
        'min': min(degrees) if degrees else 0,
        'max': max(degrees) if degrees else 0,
        'mean': sum(degrees) / len(degrees) if degrees else 0,
        'median': sorted(degrees)[len(degrees) // 2] if degrees else 0,
    }

    return {
        'num_nodes': num_nodes,
        'num_edges': num_edges,
        'is_connected': is_connected,
        'num_components': num_components,
        'degree_distribution': degree_dist,
        'star_patterns': [
            {
                'center': p.center,
                'num_leaves': len(p.leaves),
                'center_frequency': p.center_frequency,
                'is_perfect_star': p.is_perfect_star,
            }
            for p in star_patterns
        ],
        'partitions': [
            {
                'num_nodes': p.num_nodes,
                'total_samples': p.total_samples,
                'diameter': p.diameter,
            }
            for p in partitions
        ],
        'ancestral_candidates': [
            {
                'node_id': a.node_id,
                'score': a.score,
                'degree': a.degree,
                'frequency': a.frequency,
                'is_median_vector': a.is_median_vector,
            }
            for a in ancestral_candidates
        ],
        'top_central_nodes': {
            method: sorted(scores.items(), key=lambda x: x[1], reverse=True)[:5]
            for method, scores in (
                (
                    'degree',
                    {n: c['degree'] for n, c in centrality.items() if 'degree' in c},
                ),
                (
                    'betweenness',
                    {
                        n: c['betweenness']
                        for n, c in centrality.items()
                        if 'betweenness' in c
                    },
                ),
            )
        },
    }


def find_central_hub_nodes(
    network: HaplotypeNetwork, degree_threshold: Optional[int] = None
) -> List[Tuple[str, int]]:
    """
    Find hub nodes (nodes with high degree).

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.
    degree_threshold : int, optional
        Minimum degree for a node to count as a hub. If None, the
        mean degree plus one standard deviation is used.

    Returns
    -------
    List[Tuple[str, int]]
        List of (node_id, degree) tuples, sorted by degree (descending).
    """
    G = network.to_networkx()

    if G.number_of_nodes() == 0:
        return []

    degrees = dict(G.degree())

    if degree_threshold is None:
        # Calculate threshold as mean + 1 std dev
        degree_values = list(degrees.values())
        mean_degree = sum(degree_values) / len(degree_values)

        if len(degree_values) > 1:
            variance = sum((d - mean_degree) ** 2 for d in degree_values) / len(
                degree_values
            )
            std_dev = variance**0.5
            degree_threshold = int(mean_degree + std_dev)
        else:
            degree_threshold = 0

    # Find hubs
    hubs = [
        (node, degree) for node, degree in degrees.items() if degree >= degree_threshold
    ]

    # Sort by degree (descending)
    hubs.sort(key=lambda x: x[1], reverse=True)

    return hubs


def detect_bridges(network: HaplotypeNetwork) -> List[Tuple[str, str]]:
    """
    Detect bridge edges in the network.

    A bridge is an edge whose removal would disconnect the network
    or increase the number of connected components.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.

    Returns
    -------
    List[Tuple[str, str]]
        List of (node1, node2) tuples representing bridge edges.
    """
    G = network.to_networkx()

    if G.number_of_nodes() < 2:
        return []

    # Find bridges using NetworkX
    bridges = list(nx.bridges(G))

    return bridges


def identify_bottleneck_nodes(network: HaplotypeNetwork) -> List[Tuple[str, float]]:
    """
    Identify bottleneck nodes (articulation points with high betweenness).

    Bottleneck nodes are those whose removal would significantly
    disrupt information flow in the network.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.

    Returns
    -------
    List[Tuple[str, float]]
        List of (node_id, betweenness_score) tuples, sorted by score (descending).
    """
    G = network.to_networkx()

    if G.number_of_nodes() < 3:
        return []

    # Find articulation points (cut vertices)
    articulation_points = set(nx.articulation_points(G))

    if not articulation_points:
        return []

    # Calculate betweenness centrality for articulation points
    betweenness = nx.betweenness_centrality(G)

    bottlenecks = [(node, betweenness[node]) for node in articulation_points]

    # Sort by betweenness (descending)
    bottlenecks.sort(key=lambda x: x[1], reverse=True)

    return bottlenecks
