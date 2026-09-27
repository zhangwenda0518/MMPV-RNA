"""
Minimum Spanning Tree (MST) algorithm for haplotype network construction.

This module implements both Prim's and Kruskal's algorithms for constructing
minimum spanning trees from genetic sequence data. The MST forms the foundation
for more complex network algorithms including MSN and median-joining networks.

References
----------
.. [1] Excoffier, L. & Smouse, P. E. (1994). Using allele frequencies and
       geographic subdivision to reconstruct gene trees within a species:
       molecular variance parsimony. Genetics, 136(1), 343-359.
"""

import heapq
from typing import List, Optional, Tuple

from ..core.alignment import Alignment
from ..core.distance import DistanceMatrix
from ..core.graph import HaplotypeNetwork
from ..core.haplotype import identify_haplotypes_from_alignment as identify_haplotypes
from .base import NetworkAlgorithm


class MinimumSpanningTree(NetworkAlgorithm):
    """
    Construct a Minimum Spanning Tree from haplotype data.

    A minimum spanning tree (MST) connects all haplotypes with the minimum total
    genetic distance. This is the simplest haplotype network algorithm and
    forms the basis for more complex methods like MSN (Minimum Spanning Network).

    The MST is guaranteed to be a tree (no cycles) and provides the most
    parsimonious representation of relationships between haplotypes.

    Supports both Prim's and Kruskal's algorithms for MST construction:

    - **Prim's algorithm**: Grows the tree from a single starting node,
    always adding the minimum-weight edge that connects a new node.
    Time complexity: O(E log V) with binary heap.

    - **Kruskal's algorithm**: Sorts all edges and adds them in order of
    increasing weight, skipping edges that would create cycles.
    Time complexity: O(E log E) with union-find.

    Parameters
    ----------
    distance_method : str, default='hamming'
        Method for calculating pairwise distances between sequences:
        'hamming', 'jukes_cantor', 'kimura_2p' or 'tamura_nei'.
    algorithm : str, default='prim'
        MST construction algorithm to use: 'prim' or 'kruskal'.
    **kwargs : dict
        Additional parameters passed to base NetworkAlgorithm.

    Attributes
    ----------
    algorithm : str
        The selected MST algorithm.
    _distance_matrix : DistanceMatrix
        Cached distance matrix from the last construction.

    See Also
    --------
    MinimumSpanningNetwork : Extension of MST allowing alternative connections.
    TCS : Statistical parsimony network construction.

    Notes
    -----
    For most applications, Prim's algorithm is preferred as it's typically
    faster and uses less memory. Kruskal's algorithm can be advantageous
    when the graph is sparse or when edges are already sorted.

    Examples
    --------
    >>> from pypopart.algorithms import MinimumSpanningTree
    >>> from pypopart.io import load_alignment
    >>>
    >>> # Load alignment
    >>> alignment = load_alignment('sequences.fasta')
    >>>
    >>> # Construct MST using Prim's algorithm
    >>> mst = MinimumSpanningTree(algorithm='prim')
    >>> network = mst.build_network(alignment)
    >>>
    >>> # Construct using Kruskal's algorithm
    >>> mst = MinimumSpanningTree(algorithm='kruskal')
    >>> network = mst.build_network(alignment)
    """

    def __init__(
        self, distance_method: str = 'hamming', algorithm: str = 'prim', **kwargs
    ):
        """
        Initialize MST algorithm.

        Parameters
        ----------
        distance_method : str, default='hamming'
            Method for calculating distances.
        algorithm : str, default='prim'
            MST algorithm to use ('prim' or 'kruskal').
        **kwargs : dict
            Additional parameters.
        """
        super().__init__(distance_method, **kwargs)
        self.algorithm = algorithm.lower()
        if self.algorithm not in ('prim', 'kruskal'):
            raise ValueError(
                f"Unknown MST algorithm: {algorithm}. Use 'prim' or 'kruskal'"
            )

    def construct_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
        Construct MST from sequence alignment.

        Parameters
        ----------
        alignment : Alignment
            Multiple sequence alignment.
        distance_matrix : DistanceMatrix, optional
            Optional pre-computed distance matrix.

        Returns
        -------
        HaplotypeNetwork
            Haplotype network representing the MST.
        """
        # Identify unique haplotypes
        haplotypes = identify_haplotypes(alignment)

        if len(haplotypes) == 0:
            return HaplotypeNetwork()

        if len(haplotypes) == 1:
            # Single haplotype - just add it to network
            network = HaplotypeNetwork()
            network.add_haplotype(haplotypes[0])
            return network

        # Calculate distances between haplotypes
        haplotype_dist_matrix = self.calculate_haplotype_distances(haplotypes)
        self._distance_matrix = haplotype_dist_matrix

        # Build MST using selected algorithm
        if self.algorithm == 'prim':
            edges = self._prim_mst(haplotypes, haplotype_dist_matrix)
        else:  # kruskal
            edges = self._kruskal_mst(haplotypes, haplotype_dist_matrix)

        # Construct network from MST edges
        network = self._build_network(haplotypes, edges)

        return network

    def _prim_mst(
        self, haplotypes: List, distance_matrix: DistanceMatrix
    ) -> List[Tuple[str, str, float]]:
        """
        Construct MST using Prim's algorithm.

        Parameters
        ----------
        haplotypes : List
            List of Haplotype objects.
        distance_matrix : DistanceMatrix
            Distance matrix between haplotypes.

        Returns
        -------
        List[Tuple[str, str, float]]
            List of edges (id1, id2, distance).
        """
        if len(haplotypes) == 0:
            return []

        # Map haplotype IDs to haplotype objects for quick lookup
        hap_dict = {h.id: h for h in haplotypes}
        hap_ids = list(hap_dict.keys())

        # Track which nodes are in the tree
        in_tree = {hap_ids[0]}  # Start with first haplotype
        not_in_tree = set(hap_ids[1:])

        edges = []

        # Priority queue: (distance, from_id, to_id)
        pq = []

        # Add all edges from first haplotype
        for other_id in not_in_tree:
            dist = distance_matrix.get_distance(hap_ids[0], other_id)
            heapq.heappush(pq, (dist, hap_ids[0], other_id))

        # Build MST
        while not_in_tree and pq:
            dist, from_id, to_id = heapq.heappop(pq)

            # Skip if to_id already in tree (edge is outdated)
            if to_id not in not_in_tree:
                continue

            # Add edge to MST
            edges.append((from_id, to_id, dist))
            in_tree.add(to_id)
            not_in_tree.remove(to_id)

            # Add new edges from newly added node
            for other_id in not_in_tree:
                new_dist = distance_matrix.get_distance(to_id, other_id)
                heapq.heappush(pq, (new_dist, to_id, other_id))

        return edges

    def _kruskal_mst(
        self, haplotypes: List, distance_matrix: DistanceMatrix
    ) -> List[Tuple[str, str, float]]:
        """
        Construct MST using Kruskal's algorithm with Union-Find.

        Parameters
        ----------
        haplotypes : List
            List of Haplotype objects.
        distance_matrix : DistanceMatrix
            Distance matrix between haplotypes.

        Returns
        -------
        List[Tuple[str, str, float]]
            List of edges (id1, id2, distance).
        """
        if len(haplotypes) == 0:
            return []

        hap_ids = [h.id for h in haplotypes]

        # Get all edges sorted by distance
        all_edges = []
        for i, id1 in enumerate(hap_ids):
            for id2 in hap_ids[i + 1 :]:
                dist = distance_matrix.get_distance(id1, id2)
                all_edges.append((dist, id1, id2))

        all_edges.sort()  # Sort by distance

        # Union-Find data structure
        parent = {hap_id: hap_id for hap_id in hap_ids}
        rank = dict.fromkeys(hap_ids, 0)

        def find(x):
            """
            Find root of x with path compression.

            Parameters
            ----------
            x : str
                Element to look up.

            Returns
            -------
            str
                The representative element of x's set.
            """
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]

        def union(x, y):
            """
            Union sets containing x and y.

            Parameters
            ----------
            x : str
                Element to look up.
            y : str
                Second element to union with x.

            Returns
            -------
            bool
                True if the two sets were distinct and got merged.
            """
            root_x = find(x)
            root_y = find(y)

            if root_x == root_y:
                return False

            # Union by rank
            if rank[root_x] < rank[root_y]:
                parent[root_x] = root_y
            elif rank[root_x] > rank[root_y]:
                parent[root_y] = root_x
            else:
                parent[root_y] = root_x
                rank[root_x] += 1

            return True

        # Build MST
        mst_edges = []
        for dist, id1, id2 in all_edges:
            if union(id1, id2):
                mst_edges.append((id1, id2, dist))
                # Stop when we have n-1 edges
                if len(mst_edges) == len(hap_ids) - 1:
                    break

        return mst_edges

    def _build_network(
        self, haplotypes: List, edges: List[Tuple[str, str, float]]
    ) -> HaplotypeNetwork:
        """
        Build HaplotypeNetwork from haplotypes and MST edges.

        Parameters
        ----------
        haplotypes : List
            List of Haplotype objects.
        edges : List[Tuple[str, str, float]]
            List of edges (id1, id2, distance).

        Returns
        -------
        HaplotypeNetwork
            Constructed haplotype network.
        """
        network = HaplotypeNetwork()

        # Add all haplotypes
        for haplotype in haplotypes:
            network.add_haplotype(haplotype)

        # Add MST edges
        for id1, id2, dist in edges:
            network.add_edge(id1, id2, distance=int(round(dist)))

        return network

    def get_parameters(self) -> dict:
        """
        Get algorithm parameters including MST algorithm type.

        Returns
        -------
        dict
            Algorithm parameters, including the MST algorithm type.
        """
        params = super().get_parameters()
        params['algorithm'] = self.algorithm
        return params
