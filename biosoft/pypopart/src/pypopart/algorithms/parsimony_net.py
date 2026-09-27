"""
Parsimony Network algorithm, a port of PopART's AncestralSeqNet.

Builds a consensus network by repeatedly sampling ancestral-state
reconstructions on parsimony trees: each iteration picks a random tree,
samples one randomised Fitch reconstruction, and records every tree
edge as a pair of sequences. Network vertices are keyed by sequence
string, so identical ancestral states across trees collapse onto the
same vertex. Edges whose sampling frequency falls below alpha are
removed lowest-frequency-first - but only when their endpoints remain
connected without them - and ancestral vertices left with degree <= 1
are dropped.

PopART reads its trees from a Nexus TREES block; PyPopART generates
them by random-order stepwise addition (see algorithms.ancestral).

References
----------
.. [1] Excoffier, L. & Smouse, P. E. (1994). Using allele frequencies
       and geographic subdivision to reconstruct gene trees within a
       species: molecular variance parsimony. Genetics 136(1), 343-359.
"""

import random
from typing import Dict, List, Optional, Tuple

import networkx as nx

from ..core.alignment import Alignment
from ..core.distance import DistanceMatrix
from ..core.graph import HaplotypeNetwork
from ..core.haplotype import Haplotype, identify_haplotypes_from_alignment
from ..core.sequence import Sequence
from ..core.site_patterns import condense_site_patterns, is_ambiguous
from .ancestral import sample_parsimony_trees
from .base import NetworkAlgorithm


class ParsimonyNetwork(NetworkAlgorithm):
    """
    Construct a consensus haplotype network from sampled parsimony trees.

    Parameters
    ----------
    distance_method : str, default='hamming'
        Kept for interface compatibility; edge weights use PopART's
        site-weighted Hamming distances over condensed site patterns.
    n_trees : int, default=20
        Number of stepwise-addition parsimony trees to sample.
    alpha : float, default=0.95
        Edge-frequency threshold: edges sampled in fewer than alpha of
        iterations are candidates for removal (lowest first), kept only
        when removing them would disconnect their endpoints.
    n_iterations : int, optional
        Ancestral-sampling iterations. Defaults to 100 * n_trees, as in
        PopART.
    random_seed : int, optional
        Random seed for reproducibility.
    **kwargs : dict
        Additional parameters passed to the base class.
    """

    def __init__(
        self,
        distance_method: str = 'hamming',
        n_trees: int = 20,
        alpha: float = 0.95,
        n_iterations: Optional[int] = None,
        random_seed: Optional[int] = None,
        **kwargs,
    ):
        """
        Initialize Parsimony Network algorithm.

        Parameters
        ----------
        distance_method : str, default='hamming'
            Method for calculating distances.
        n_trees : int, default=20
            Number of parsimony trees to sample.
        alpha : float, default=0.95
            Edge-frequency threshold for pruning.
        n_iterations : int, optional
            Sampling iterations (default 100 * n_trees).
        random_seed : int, optional
            Random seed for reproducibility.
        **kwargs : dict
            Additional parameters.
        """
        super().__init__(distance_method, **kwargs)
        self.n_trees = n_trees
        self.alpha = alpha
        self.n_iterations = n_iterations
        self.random_seed = random_seed

    def construct_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
        Construct the parsimony consensus network.

        Parameters
        ----------
        alignment : Alignment
            Multiple sequence alignment.
        distance_matrix : DistanceMatrix, optional
            Ignored; distances are computed in condensed space.

        Returns
        -------
        HaplotypeNetwork
            Consensus network; inferred ancestral vertices are flagged
            as median vectors and carry their condensed-space sequence.
        """
        haplotypes = identify_haplotypes_from_alignment(alignment)

        if len(haplotypes) == 0:
            return HaplotypeNetwork()
        if len(haplotypes) == 1:
            network = HaplotypeNetwork()
            network.add_haplotype(haplotypes[0])
            return network

        rng = random.Random(self.random_seed)
        condensed, weights = condense_site_patterns([h.data for h in haplotypes])

        trees = sample_parsimony_trees(condensed, weights, self.n_trees, rng)
        n_iterations = (
            self.n_iterations if self.n_iterations is not None else 100 * self.n_trees
        )

        graph, seq_of_vertex, n_samples = self._sample_consensus(
            condensed, weights, trees, n_iterations, rng
        )
        self._prune_by_frequency(graph, n_iterations)
        self._drop_isolated_ancestors(graph, n_samples)

        # Build the HaplotypeNetwork
        network = HaplotypeNetwork()
        for idx, haplotype in enumerate(haplotypes):
            if idx in graph:
                network.add_haplotype(haplotype)
        label_of = {idx: haplotypes[idx].id for idx in range(n_samples)}
        ancestor_counter = 0
        for vertex in sorted(graph.nodes):
            if vertex < n_samples:
                continue
            label = f'Median_{ancestor_counter}'
            ancestor_counter += 1
            label_of[vertex] = label
            seq = Sequence(
                id=label,
                data=seq_of_vertex[vertex],
                description='Sampled ancestral sequence (condensed sites)',
            )
            network.add_haplotype(
                Haplotype(sequence=seq, sample_ids=[]), median_vector=True
            )
            network.graph.nodes[label]['is_median'] = True
        for u, v, attrs in graph.edges(data=True):
            network.add_edge(label_of[u], label_of[v], distance=attrs['distance'])

        return network

    def _sample_consensus(
        self,
        condensed: List[str],
        weights: List[int],
        trees,
        n_iterations: int,
        rng: random.Random,
    ) -> Tuple[nx.Graph, Dict[int, str], int]:
        """
        Sample tree edges into a frequency-annotated consensus graph.

        Port of AncestralSeqNet::computeGraph's sampling loop: vertices
        are keyed by sequence string; per iteration, each distinct
        vertex and edge seen is counted once.

        Parameters
        ----------
        condensed : list of str
            Condensed sample sequences.
        weights : list of int
            Site weights.
        trees : list of ParsimonyTree
            Sampled parsimony trees.
        n_iterations : int
            Number of ancestral samplings.
        rng : random.Random
            Random source.

        Returns
        -------
        tuple
            (graph with 'count' and 'distance' edge attrs, vertex id to
            sequence mapping, number of sample vertices).
        """
        graph = nx.Graph()
        seq_to_vertex: Dict[str, int] = {}
        seq_of_vertex: Dict[int, str] = {}

        for idx, seq in enumerate(condensed):
            seq_to_vertex.setdefault(seq, idx)
            seq_of_vertex[idx] = seq
            graph.add_node(idx)
        n_samples = len(condensed)
        next_vertex = n_samples

        for _ in range(n_iterations):
            tree = trees[rng.randrange(len(trees))]
            ancestors = tree.sample_ancestors(rng)
            edge_list = tree.edge_sequences(ancestors)

            vertices_seen = set()
            edges_seen = set()
            for seq_from, seq_to in edge_list:
                u = seq_to_vertex.get(seq_from)
                if u is None:
                    u = next_vertex
                    next_vertex += 1
                    seq_to_vertex[seq_from] = u
                    seq_of_vertex[u] = seq_from
                    graph.add_node(u)
                v = seq_to_vertex.get(seq_to)
                if v is None:
                    v = next_vertex
                    next_vertex += 1
                    seq_to_vertex[seq_to] = v
                    seq_of_vertex[v] = seq_to
                    graph.add_node(v)

                vertices_seen.update((u, v))
                if u != v:
                    if not graph.has_edge(u, v):
                        distance = self._weighted_distance(seq_from, seq_to, weights)
                        graph.add_edge(u, v, distance=distance, count=0)
                    edges_seen.add(frozenset((u, v)))

            for edge in edges_seen:
                u, v = tuple(edge)
                graph[u][v]['count'] += 1

        return graph, seq_of_vertex, n_samples

    def _prune_by_frequency(self, graph: nx.Graph, n_iterations: int) -> None:
        """
        Remove low-frequency edges unless removal disconnects endpoints.

        Edges with sampling frequency below alpha are processed lowest
        frequency first; each is removed only when a path between its
        endpoints survives, matching the C++ areConnected guard.

        Parameters
        ----------
        graph : nx.Graph
            Consensus graph (modified in place).
        n_iterations : int
            Total sampling iterations (frequency denominator).
        """
        candidates = sorted(
            (attrs['count'] / n_iterations, u, v)
            for u, v, attrs in graph.edges(data=True)
            if attrs['count'] / n_iterations < self.alpha
        )
        for _freq, u, v in candidates:
            attrs = dict(graph[u][v])
            graph.remove_edge(u, v)
            if not nx.has_path(graph, u, v):
                graph.add_edge(u, v, **attrs)

    @staticmethod
    def _drop_isolated_ancestors(graph: nx.Graph, n_samples: int) -> None:
        """
        Drop ancestral vertices of degree <= 1 (descending pass).

        Parameters
        ----------
        graph : nx.Graph
            Consensus graph (modified in place).
        n_samples : int
            Number of sample vertices (never removed).
        """
        for vertex in sorted((v for v in graph.nodes if v >= n_samples), reverse=True):
            if graph.degree(vertex) <= 1:
                graph.remove_node(vertex)

    @staticmethod
    def _weighted_distance(s1: str, s2: str, weights: List[int]) -> float:
        """
        Site-weighted Hamming distance over condensed sequences.

        Parameters
        ----------
        s1 : str
            First condensed sequence.
        s2 : str
            Second condensed sequence.
        weights : list of int
            Site weights per condensed column.

        Returns
        -------
        float
            Weighted distance (ambiguous positions skipped).
        """
        total = 0
        for c1, c2, w in zip(s1, s2, weights):
            if is_ambiguous(c1) or is_ambiguous(c2):
                continue
            if c1 != c2:
                total += w
        return float(total)

    def get_parameters(self) -> dict:
        """
        Get algorithm parameters.

        Returns
        -------
        dict
            Parameters including n_trees, alpha, n_iterations and seed.
        """
        params = super().get_parameters()
        params['n_trees'] = self.n_trees
        params['alpha'] = self.alpha
        params['n_iterations'] = self.n_iterations
        params['random_seed'] = self.random_seed
        return params
