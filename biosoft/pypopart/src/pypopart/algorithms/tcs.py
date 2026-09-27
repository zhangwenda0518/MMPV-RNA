"""
TCS (Statistical Parsimony) algorithm - Enhanced to match C++ PopART implementation.

This implementation closely follows the C++ PopART TCS.cpp logic:
- Component-based connection algorithm
- Intermediate sequence inference with scoring
- Post-processing vertex collapse
- Floyd-Warshall for path length calculations

Implements the method from Clement, Posada & Crandall (2000):
"TCS: a computer program to estimate gene genealogies"
Molecular Ecology 9: 1657-1659
"""

import math
from typing import Dict, List, Optional, Tuple

import numpy as np

from ..core.alignment import Alignment
from ..core.distance import DistanceMatrix
from ..core.graph import HaplotypeNetwork
from ..core.haplotype import Haplotype
from ..core.haplotype import identify_haplotypes_from_alignment as identify_haplotypes
from ..core.sequence import Sequence
from .base import NetworkAlgorithm


class TCS(NetworkAlgorithm):
    """
    Construct haplotype network using Statistical Parsimony (TCS algorithm).

    This implementation accurately reproduces the C++ PopART TCS algorithm:
    1. Groups haplotypes by pairwise distances
    2. Iteratively connects components at increasing distance levels
    3. Infers intermediate sequences when distance > 1
    4. Uses scoring system to find optimal intermediates
    5. Post-processes to collapse degree-2 vertices

    Parameters
    ----------
    distance_method : str, default='hamming'
        Method for calculating distances (should be hamming).
    confidence : float, default=0.95
        Confidence level for the 'auto' connection limit.
    connection_limit : int or 'auto', optional
        Optional connection cap; None (default) matches PopART's
        unlimited behaviour.
    infer_intermediates : bool, default=True
        Whether to infer intermediate sequences.
    collapse_vertices : bool, default=True
        Whether to collapse degree-2 intermediate vertices.
    **kwargs : dict
        Additional parameters passed to the base class.
    """

    # Scoring constants from C++ implementation
    BONUS = 20
    SHORTCUTPENALTY = 10
    LONGPENALTY = 5

    def __init__(
        self,
        distance_method: str = 'hamming',
        confidence: float = 0.95,
        connection_limit=None,
        infer_intermediates: bool = True,
        collapse_vertices: bool = True,
        **kwargs,
    ):
        """
        Initialize TCS algorithm.

        Parameters
        ----------
        distance_method : str, default='hamming'
            Method for calculating distances (should be hamming).
        confidence : float, default=0.95
            Confidence level for the 'auto' connection limit. Only used
            when connection_limit='auto'.
        connection_limit : int or 'auto', optional
            PopART's TCS has no connection limit and always produces a
            fully connected network; None (the default) matches that.
            Pass an int to cap connections at that distance, or 'auto'
            to derive a cap from `confidence` (PyPopART-specific extras,
            both off by default).
        infer_intermediates : bool, default=True
            Whether to infer intermediate sequences.
        collapse_vertices : bool, default=True
            Whether to collapse degree-2 intermediate vertices.
        **kwargs : dict
            Additional parameters.
        """
        super().__init__(distance_method, **kwargs)
        self.confidence = confidence
        self.connection_limit = connection_limit
        self.infer_intermediates = infer_intermediates
        self.collapse_vertices = collapse_vertices
        self._intermediate_counter = 0
        self._path_cache: Dict[str, Dict[str, float]] = {}
        self._path_cache_version: Optional[int] = None

        if distance_method not in ('hamming',):
            import warnings

            warnings.warn(
                f"TCS is designed for hamming distance. Using '{distance_method}' "
                'may not produce theoretically correct results.',
                stacklevel=2,
            )

    def construct_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
        Construct TCS network from sequence alignment.

        Implements the component-based connection algorithm from C++ TCS.cpp.

        Parameters
        ----------
        alignment : Alignment
            Multiple sequence alignment.
        distance_matrix : DistanceMatrix, optional
            Optional pre-computed distance matrix.

        Returns
        -------
        HaplotypeNetwork
            Haplotype network constructed using statistical parsimony.
        """
        # Identify unique haplotypes
        haplotypes = identify_haplotypes(alignment)

        if len(haplotypes) == 0:
            return HaplotypeNetwork()

        if len(haplotypes) == 1:
            network = HaplotypeNetwork()
            network.add_haplotype(haplotypes[0])
            return network

        # Calculate distances between haplotypes
        haplotype_dist_matrix = self.calculate_haplotype_distances(haplotypes)
        self._distance_matrix = haplotype_dist_matrix

        # PopART's TCS has no connection limit; a cap is opt-in
        if self.connection_limit == 'auto':
            limit = self._calculate_connection_limit(alignment.length, len(haplotypes))
        else:
            limit = self.connection_limit

        # Build network using component-based algorithm (matches C++ TCS.cpp)
        network = self._build_network_with_components(
            haplotypes, haplotype_dist_matrix, alignment.length, limit
        )

        # Collapse degree-2 vertices (post-processing simplification)
        if self.collapse_vertices:
            network = self._collapse_degree2_vertices(network)

        return network

    def _calculate_connection_limit(
        self, sequence_length: int, num_haplotypes: int
    ) -> int:
        """
        Calculate maximum parsimony connection limit.

        Uses the formula from Templeton et al. (1992) to estimate the
        maximum number of mutational differences that can be explained
        by parsimony at the specified confidence level.

        Parameters
        ----------
        sequence_length : int
            Length of aligned sequences.
        num_haplotypes : int
            Number of unique haplotypes.

        Returns
        -------
        int
            Maximum connection distance for parsimony criterion.
        """
        connection_limit = 1
        cumulative_prob = 0.0

        for k in range(1, sequence_length + 1):
            prob_k = self._poisson_probability(k, sequence_length, num_haplotypes)
            cumulative_prob += prob_k

            if cumulative_prob >= (1 - self.confidence):
                connection_limit = k
                break

        return max(1, connection_limit)

    def _poisson_probability(self, k: int, seq_length: int, sample_size: int) -> float:
        """
        Calculate probability using Poisson distribution.

        Parameters
        ----------
        k : int
            Number of mutations.
        seq_length : int
            Sequence length.
        sample_size : int
            Number of sequences.

        Returns
        -------
        float
            Probability.
        """
        lambda_param = 2.0 * math.log(sample_size) if sample_size > 1 else 1.0

        try:
            prob = (lambda_param**k) * math.exp(-lambda_param) / math.factorial(k)
        except (OverflowError, ValueError):
            prob = 0.0

        return prob

    def _build_network_with_components(
        self,
        haplotypes: List,
        distance_matrix: DistanceMatrix,
        sequence_length: int,
        connection_limit: Optional[int] = None,
    ) -> HaplotypeNetwork:
        """
        Build network using component-based algorithm from C++ TCS.cpp.

        This accurately reproduces the C++ logic:
        1. Create all haplotypes as vertices, each in its own component
        2. Group pairs by distance
        3. For each distance level (1 to connection_limit):
           - Process pairs at that distance
           - Connect vertices from different components
           - Infer intermediates for distance > 1
           - Merge components when connected
        4. Clean up intermediate vertices in post-processing

        Parameters
        ----------
        haplotypes : list of Haplotype
            List of Haplotype objects.
        distance_matrix : DistanceMatrix, optional
            Distance matrix.
        sequence_length : int
            Length of sequences for creating intermediates.
        connection_limit : int, optional
            Optional maximum connection distance (None = no limit,
            matching PopART).

        Returns
        -------
        HaplotypeNetwork
            Network with inferred intermediate sequences.
        """
        # Initialize network with all haplotypes
        network = HaplotypeNetwork()
        for haplotype in haplotypes:
            network.add_haplotype(haplotype)

        # Component tracking: maps haplotype ID to component ID
        # Each haplotype starts in its own component
        component_ids: Dict[str, int] = {h.id: i for i, h in enumerate(haplotypes)}

        # Group pairs by distance (like C++ VertContainer and priority queue)
        pairs_by_distance: Dict[int, List[Tuple[str, str]]] = {}

        for i, h1 in enumerate(haplotypes):
            for j in range(i + 1, len(haplotypes)):
                h2 = haplotypes[j]
                dist = int(round(distance_matrix.get_distance(h1.id, h2.id)))

                if connection_limit is None or dist <= connection_limit:
                    if dist not in pairs_by_distance:
                        pairs_by_distance[dist] = []
                    pairs_by_distance[dist].append((h1.id, h2.id))

        # Process pairs in order of increasing distance (like C++ priority queue)
        for M in sorted(pairs_by_distance.keys()):
            # Keep processing this distance level until no more pairs remain
            while M in pairs_by_distance and len(pairs_by_distance[M]) > 0:
                pairs = pairs_by_distance[M]

                # Track which component pair we're working on
                comp_a = -1
                comp_b = -1
                other_pairs = []

                for u_id, v_id in pairs:
                    comp_u = component_ids.get(u_id, -1)
                    comp_v = component_ids.get(v_id, -1)

                    # Skip if already in same component
                    if comp_u == comp_v:
                        continue

                    # Ensure comp_u < comp_v
                    if comp_u > comp_v:
                        comp_u, comp_v = comp_v, comp_u
                        u_id, v_id = v_id, u_id

                    # Set component pair on first distinct pair
                    if comp_a < 0:
                        comp_a = comp_u
                        comp_b = comp_v

                    # Process pairs for the current component pair
                    if comp_u == comp_a and comp_v == comp_b:
                        if M == 1:
                            # Direct connection for distance 1
                            network.add_edge(u_id, v_id, distance=1)
                        else:
                            # Infer intermediates for distance > 1
                            if self.infer_intermediates:
                                self._add_connection_with_intermediates(
                                    network,
                                    u_id,
                                    v_id,
                                    M,
                                    sequence_length,
                                    component_ids,
                                    comp_u,
                                    comp_v,
                                )
                            else:
                                # Simple connection without intermediates
                                network.add_edge(u_id, v_id, distance=M)
                    else:
                        # Save for next iteration of this distance level
                        other_pairs.append((u_id, v_id))

                # Merge components (like C++ component renumbering)
                if comp_a >= 0:
                    for hap_id in list(component_ids.keys()):
                        if component_ids[hap_id] < 0 or component_ids[hap_id] == comp_b:
                            component_ids[hap_id] = comp_a
                        elif component_ids[hap_id] > comp_b:
                            component_ids[hap_id] -= 1

                # Update pairs for this distance level
                if other_pairs:
                    pairs_by_distance[M] = other_pairs
                else:
                    # No more pairs at this distance, remove from dict
                    del pairs_by_distance[M]

        return network

    @staticmethod
    def _all_path_lengths(
        network: HaplotypeNetwork, sources: List[str]
    ) -> Dict[str, Dict[str, float]]:
        """
        Compute weighted shortest-path lengths from the given sources.

        One Dijkstra sweep per relevant source, computed once per
        network state; the per-pair path queries in findIntermediates
        and computeScore read from this table (C++ caches a
        Floyd-Warshall the same way). Restricting sources to the two
        components being joined (plus no-man's-land intermediates)
        keeps this cheap while components are small.

        Parameters
        ----------
        network : HaplotypeNetwork
            Current network.
        sources : list of str
            Vertices to run Dijkstra from.

        Returns
        -------
        dict
            Mapping source -> {target: total distance}; missing targets
            are unreachable.
        """
        import networkx as nx

        graph = network.graph
        return {
            source: nx.single_source_dijkstra_path_length(
                graph, source, weight='distance'
            )
            for source in sources
            if source in graph
        }

    def _add_connection_with_intermediates(
        self,
        network: HaplotypeNetwork,
        u_id: str,
        v_id: str,
        distance: int,
        sequence_length: int,
        component_ids: Dict[str, int],
        comp_u: int,
        comp_v: int,
    ) -> None:
        """
        Add connection between u and v, inferring intermediates if needed.

        Implements findIntermediates() and newCompositePath() from C++ TCS.cpp.

        Parameters
        ----------
        network : HaplotypeNetwork
            Current network.
        u_id : str
            Source haplotype ID.
        v_id : str
            Target haplotype ID.
        distance : int
            Distance between u and v.
        sequence_length : int
            Length of sequences.
        component_ids : dict
            Component membership tracker.
        comp_u : int
            Component of u.
        comp_v : int
            Component of v.
        """
        # Find optimal intermediate vertices (like C++ findIntermediates)
        # Only vertices in the two components being joined (or in
        # no-man's-land) can appear as path endpoints. Sweeps are cached
        # per network version and recomputed only after topology changes.
        version = network._version
        if self._path_cache_version != version:
            self._path_cache = {}
            self._path_cache_version = version
        # Only paths *from* the connection endpoints and the sampled
        # members of the two components are ever read (candidate paths
        # are read from member rows by symmetry), so those are the only
        # Dijkstra sources needed.
        index_of = self._distance_matrix._label_index
        sources = {u_id, v_id}
        sources.update(
            node_id
            for node_id, comp in component_ids.items()
            if comp in (comp_u, comp_v) and node_id in index_of
        )
        missing = [s for s in sources if s not in self._path_cache]
        self._path_cache.update(self._all_path_lengths(network, missing))
        path_lengths = self._path_cache
        int_u, int_v, min_path_length = self._find_intermediates(
            network, u_id, v_id, distance, component_ids, comp_u, comp_v, path_lengths
        )

        # Check if path already exists
        existing_path_length = path_lengths.get(int_u, {}).get(int_v, float('inf'))

        if existing_path_length < min_path_length:
            # C++ TCS treats this as an internal inconsistency
            raise RuntimeError('Shorter path already exists between these vertices!')

        # Only add the new path when none of the right length exists
        if existing_path_length > min_path_length:
            self._create_composite_path(
                network, int_u, int_v, min_path_length, sequence_length, component_ids
            )

    def _find_intermediates(
        self,
        network: HaplotypeNetwork,
        u_id: str,
        v_id: str,
        dist: int,
        component_ids: Dict[str, int],
        comp_u: int,
        comp_v: int,
        path_lengths: Dict[str, Dict[str, float]],
    ) -> Tuple[str, str, int]:
        """
        Find optimal intermediate vertices to connect two components.

        Implements C++ TCS::findIntermediates() with scoring system.

        Parameters
        ----------
        network : HaplotypeNetwork
            Current network.
        u_id : str
            Source haplotype ID from component comp_u.
        v_id : str
            Target haplotype ID from component comp_v.
        dist : int
            Distance between components.
        component_ids : dict
            Component membership tracker.
        comp_u : int
            Source component.
        comp_v : int
            Target component.
        path_lengths : dict
            All-pairs weighted path lengths from _all_path_lengths.

        Returns
        -------
        tuple
            Tuple of (intermediate_u_id, intermediate_v_id, path_length).
        """
        max_score = float('-inf')
        min_path_length = dist
        best_u = u_id
        best_v = v_id

        score_context = self._build_score_context(
            comp_u, comp_v, component_ids, path_lengths
        )

        # Try all vertices in comp_u (or "no man's land" with comp_id < 0)
        for i_id in list(component_ids.keys()):
            if (
                component_ids.get(i_id, -1) != comp_u
                and component_ids.get(i_id, -1) >= 0
            ):
                continue

            # Check if connected to u
            path_ui = path_lengths.get(u_id, {}).get(i_id)
            if path_ui is None:
                continue

            if path_ui >= dist:
                continue

            # Try all vertices in comp_v
            for j_id in list(component_ids.keys()):
                if (
                    component_ids.get(j_id, -1) != comp_v
                    and component_ids.get(j_id, -1) >= 0
                ):
                    continue

                # Check if connected to v
                path_vj = path_lengths.get(v_id, {}).get(j_id)
                if path_vj is None:
                    continue

                if path_vj + path_ui >= dist:
                    continue

                dP = dist - path_vj - path_ui
                score = self._compute_score(i_id, j_id, dP, dist, score_context)

                # Select best scoring pair (or shortest path if tied)
                if score > max_score or (score == max_score and dP < min_path_length):
                    min_path_length = dP
                    max_score = score
                    best_u = i_id
                    best_v = j_id

        return best_u, best_v, min_path_length

    @staticmethod
    def _build_score_context(
        comp_u: int,
        comp_v: int,
        component_ids: Dict[str, int],
        path_lengths: Dict[str, Dict[str, float]],
    ) -> Dict:
        """
        Precompute the shared inputs of computeScore for one merge.

        Member lists, the original-distance submatrix, and per-candidate
        path vectors (cached lazily) are hoisted out of the candidate
        loop, which calls _compute_score once per candidate pair.

        Parameters
        ----------
        comp_u : int
            Source component.
        comp_v : int
            Target component.
        component_ids : dict
            Component membership tracker.
        path_lengths : dict
            Weighted path lengths keyed by source.

        Returns
        -------
        dict
            Context consumed by _compute_score.
        """
        return {
            'members_u': [h for h, c in component_ids.items() if c == comp_u],
            'members_v': [h for h, c in component_ids.items() if c == comp_v],
            'path_lengths': path_lengths,
            'orig': None,
            'pu_cache': {},
            'pv_cache': {},
        }

    def _compute_score(
        self, u_id: str, v_id: str, dP: int, clust_dist: int, context: Dict
    ) -> float:
        """
        Compute the score for one intermediate pair (C++ computeScore).

        Vectorised over all (i, j) pairs of original haplotypes in the
        two components: +BONUS when the path through the intermediates
        matches the original distance, -LONGPENALTY when longer,
        -SHORTCUTPENALTY when shorter, and -inf when a shortcut
        undercuts the cluster distance.

        Parameters
        ----------
        u_id : str
            Candidate intermediate on the comp_u side.
        v_id : str
            Candidate intermediate on the comp_v side.
        dP : int
            Path length between the candidate intermediates.
        clust_dist : int
            Distance between the original components.
        context : dict
            Precomputed inputs from _build_score_context.

        Returns
        -------
        float
            Score for this intermediate pair.
        """
        index_of = self._distance_matrix._label_index

        if context['orig'] is None:
            context['members_u'] = [h for h in context['members_u'] if h in index_of]
            context['members_v'] = [h for h in context['members_v'] if h in index_of]
            context['orig'] = self._distance_matrix.matrix[
                np.ix_(
                    [index_of[h] for h in context['members_u']],
                    [index_of[h] for h in context['members_v']],
                )
            ]

        members_u = context['members_u']
        members_v = context['members_v']
        if not members_u or not members_v:
            return 0.0

        path_lengths = context['path_lengths']
        pu = context['pu_cache'].get(u_id)
        if pu is None:
            # Undirected graph: path(candidate -> member) is read from
            # the member's row, so candidates need no Dijkstra sweep
            pu = np.array(
                [path_lengths.get(h, {}).get(u_id, np.inf) for h in members_u]
            )
            context['pu_cache'][u_id] = pu
        pv = context['pv_cache'].get(v_id)
        if pv is None:
            pv = np.array(
                [path_lengths.get(h, {}).get(v_id, np.inf) for h in members_v]
            )
            context['pv_cache'][v_id] = pv

        total = dP + pu[:, None] + pv[None, :]
        orig = context['orig']

        valid = np.isfinite(total)
        close = valid & (np.abs(total - orig) < 0.5)
        longer = valid & ~close & (total > orig)
        shorter = valid & ~close & ~longer

        if bool((shorter & (total < clust_dist)).any()):
            return float('-inf')  # Invalid shortcut

        return float(
            self.BONUS * close.sum()
            - self.LONGPENALTY * longer.sum()
            - self.SHORTCUTPENALTY * shorter.sum()
        )

    def _create_composite_path(
        self,
        network: HaplotypeNetwork,
        start_id: str,
        end_id: str,
        distance: int,
        sequence_length: int,
        component_ids: Dict[str, int],
    ) -> None:
        """
        Create path of intermediate vertices connecting start to end.

        Implements C++ TCS::newCompositePath().

        Parameters
        ----------
        network : HaplotypeNetwork
            Current network.
        start_id : str
            Starting vertex ID.
        end_id : str
            Ending vertex ID.
        distance : int
            Number of intermediates to create.
        sequence_length : int
            Sequence length for intermediates.
        component_ids : dict
            Component membership tracker.
        """
        current_id = start_id

        # Create distance-1 intermediate vertices
        for _i in range(1, distance):
            intermediate_id = f'intermediate_{self._intermediate_counter}'
            self._intermediate_counter += 1

            # Unlabelled vertex: no sequence, as in C++ newCompositePath
            intermediate_seq = Sequence(
                id=intermediate_id,
                data='',
                description='Inferred intermediate vertex',
            )

            intermediate_hap = Haplotype(
                sequence=intermediate_seq,
                sample_ids=[],
            )

            # Add to network as a median vector, so the GUI draws it as an
            # unlabelled vertex like PopART does, and marked for collapse.
            network.add_haplotype(intermediate_hap, median_vector=True)
            network.graph.nodes[intermediate_id]['is_intermediate'] = True
            network.add_edge(current_id, intermediate_id, distance=1)

            # Mark as "no man's land" (component ID = -1)
            component_ids[intermediate_id] = -1

            current_id = intermediate_id

        # Connect last intermediate to end
        network.add_edge(current_id, end_id, distance=1)

    def _collapse_degree2_vertices(self, network: HaplotypeNetwork) -> HaplotypeNetwork:
        """
        Collapse vertices with degree 2 (post-processing simplification).

        Implements C++ TCS post-processing that removes degree-2 vertices
        connecting only two other vertices. This matches lines 161-194 of TCS.cpp.

        Parameters
        ----------
        network : HaplotypeNetwork
            Network with potential degree-2 vertices.

        Returns
        -------
        HaplotypeNetwork
            Simplified network.
        """
        # Process one vertex at a time (like C++ implementation)
        # Keep looping until no more collapses possible
        changed = True

        while changed:
            changed = False

            # Find and collapse ONE degree-2 intermediate vertex
            for hap_id in list(network.nodes):
                # Skip if already removed
                if not network.has_node(hap_id):
                    continue

                degree = network.get_degree(hap_id)

                # Only collapse degree-2 vertices
                if degree != 2:
                    continue

                neighbors = network.get_neighbors(hap_id)

                if len(neighbors) != 2:
                    continue

                n1, n2 = neighbors

                # Get edge weights
                try:
                    w1 = network.get_edge_distance(hap_id, n1)
                    w2 = network.get_edge_distance(hap_id, n2)
                except Exception:
                    continue

                # Only collapse intermediates (not original haplotypes);
                # C++ checks vertex index >= nseqs, we use the marker attr
                if network.graph.nodes[hap_id].get('is_intermediate'):
                    try:
                        combined_weight = w1 + w2

                        # Remove the intermediate vertex and its edges
                        network.remove_haplotype(hap_id)

                        # Add direct edge if it doesn't exist
                        if not network.has_edge(n1, n2):
                            network.add_edge(n1, n2, distance=combined_weight)

                        # Mark that we made a change - restart loop
                        changed = True
                        break  # Exit for loop and restart while loop
                    except Exception:
                        # Skip if removal fails
                        pass

        return network

    def get_parameters(self) -> dict:
        """
        Get algorithm parameters.

        Returns
        -------
        dict
            Dictionary of algorithm parameters.
        """
        params = super().get_parameters()
        params['confidence'] = self.confidence
        params['connection_limit'] = self.connection_limit
        params['infer_intermediates'] = self.infer_intermediates
        params['collapse_vertices'] = self.collapse_vertices
        return params
