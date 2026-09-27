"""Minimum Spanning Network (MSN) algorithm for haplotype networks."""

from collections import defaultdict
from typing import List, Optional, Tuple

from ..core.alignment import Alignment
from ..core.distance import DistanceMatrix
from ..core.graph import HaplotypeNetwork
from ..core.haplotype import identify_haplotypes_from_alignment as identify_haplotypes
from .mst import MinimumSpanningTree


class MinimumSpanningNetwork(MinimumSpanningTree):
    """
    Construct a Minimum Spanning Network (MSN) from haplotype data.

    Follows PopART's AbstractMSN::computeMSN: pairs are processed in
    ascending distance levels; at each level, edges are added for every
    pair joining two distinct components (strict mode, epsilon == 0) or
    for every pair at that level (relaxed mode, epsilon > 0). Component
    merging happens only after a whole level is added, so ties at the
    connection threshold produce the MSN's alternative paths. Once the
    network first becomes connected at threshold T, levels are processed
    up to T + epsilon and then construction stops.

    Parameters
    ----------
    distance_method : str, default='hamming'
        Method for calculating distances.
    epsilon : float, default=0.0
        Extension of the connection threshold: levels up to
        (first-connection threshold + epsilon) are included, and any
        epsilon > 0 relaxes the cross-component requirement, as in
        PopART.
    prune_redundant : bool, default=False
        PyPopART-specific extra (no PopART analogue, off by default):
        after construction, remove edges for which an alternative path
        of equal or shorter total distance exists.
    **kwargs : dict
        Additional parameters passed to the base class.
    """

    def __init__(
        self,
        distance_method: str = 'hamming',
        epsilon: float = 0.0,
        prune_redundant: bool = False,
        **kwargs,
    ):
        """
        Initialize MSN algorithm.

        Parameters
        ----------
        distance_method : str, default='hamming'
            Method for calculating distances.
        epsilon : float, default=0.0
            Connection-threshold extension (PopART semantics; a value
            greater than zero also relaxes the cross-component rule).
        prune_redundant : bool, default=False
            Opt-in removal of edges with equal-or-shorter alternative
            paths. Not part of PopART's MSN.
        **kwargs : dict
            Additional parameters passed to the base class.
        """
        super().__init__(distance_method, algorithm='prim', **kwargs)
        self.epsilon = epsilon
        self.prune_redundant = prune_redundant

    def construct_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
        Construct MSN from sequence alignment.

        Parameters
        ----------
        alignment : Alignment
            Multiple sequence alignment.
        distance_matrix : DistanceMatrix, optional
            Optional pre-computed distance matrix.

        Returns
        -------
        HaplotypeNetwork
            Haplotype network representing the MSN.
        """
        haplotypes = identify_haplotypes(alignment)

        if len(haplotypes) <= 1:
            return super().construct_network(alignment, distance_matrix)

        haplotype_dist_matrix = self.calculate_haplotype_distances(haplotypes)
        self._distance_matrix = haplotype_dist_matrix

        labels = [h.id for h in haplotypes]
        edges = self._msn_edges(labels, haplotype_dist_matrix)

        network = HaplotypeNetwork()
        for haplotype in haplotypes:
            network.add_haplotype(haplotype)
        for u, v, dist in edges:
            network.add_edge(u, v, distance=dist)

        if self.prune_redundant:
            self._prune_redundant_edges(network)

        return network

    def _msn_edges(
        self, labels: List[str], distance_matrix: DistanceMatrix
    ) -> List[Tuple[str, str, float]]:
        """
        Compute MSN edges with PopART's level-sweep algorithm.

        Parameters
        ----------
        labels : list of str
            Node labels, in distance-matrix order.
        distance_matrix : DistanceMatrix
            Pairwise distances between the labelled nodes.

        Returns
        -------
        list of tuple
            Edges as (label_u, label_v, distance).

        Raises
        ------
        RuntimeError
            If the pair queue empties before the graph is connected
            (cannot happen with a complete finite distance matrix).
        """
        n = len(labels)
        matrix = distance_matrix.matrix
        strict = self.epsilon == 0

        levels = defaultdict(list)
        for i in range(n):
            for j in range(i):
                levels[float(matrix[i, j])].append((i, j))

        component = list(range(n))
        ncomps = n
        max_value = float('inf')
        edges: List[Tuple[str, str, float]] = []

        for dist in sorted(levels):
            if dist > max_value:
                break

            # Add edges for the whole level before merging components,
            # so equal-distance ties between the same components all
            # make it into the network.
            new_pairs = []
            for i, j in levels[dist]:
                if not strict or component[i] != component[j]:
                    new_pairs.append((i, j))
                    edges.append((labels[i], labels[j], dist))

            for i, j in new_pairs:
                comp_i, comp_j = component[i], component[j]
                if comp_i != comp_j:
                    low, high = min(comp_i, comp_j), max(comp_i, comp_j)
                    component = [
                        low if c == high else (c - 1 if c > high else c)
                        for c in component
                    ]
                    ncomps -= 1
                if ncomps == 1 and max_value == float('inf'):
                    max_value = dist + self.epsilon

        if ncomps > 1:
            raise RuntimeError('Pair queue empty before the graph is connected')

        return edges

    def _prune_redundant_edges(self, network: HaplotypeNetwork) -> None:
        """
        Remove edges that have an equal-or-shorter alternative path.

        PyPopART-specific opt-in behaviour (PopART keeps all tied edges).
        Edges are examined longest-first; an edge is dropped when the
        remaining network still offers a path between its endpoints of
        no greater total distance.

        Parameters
        ----------
        network : HaplotypeNetwork
            Network to prune in place.
        """
        import networkx as nx

        graph = network.graph
        edges_by_length = sorted(
            graph.edges(data=True),
            key=lambda e: e[2].get('distance', 0),
            reverse=True,
        )
        for u, v, attrs in edges_by_length:
            dist = attrs.get('distance', 0)
            graph.remove_edge(u, v)
            try:
                alternative = nx.shortest_path_length(
                    graph, u, v, weight=lambda a, b, d: d.get('distance') or 1
                )
            except nx.NetworkXNoPath:
                alternative = float('inf')
            if alternative > dist:
                graph.add_edge(u, v, **attrs)

    def get_parameters(self) -> dict:
        """
        Get algorithm parameters.

        Returns
        -------
        dict
            Parameters including epsilon and prune_redundant.
        """
        params = super().get_parameters()
        params['epsilon'] = self.epsilon
        params['prune_redundant'] = self.prune_redundant
        return params
