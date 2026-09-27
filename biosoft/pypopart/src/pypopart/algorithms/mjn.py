"""
Median-Joining Network (MJN) algorithm, ported from PopART's MedJoinNet.

Implements Bandelt, Forster & Röhl (1999) as realised in the C++ code:
the algorithm runs in condensed site-pattern space (see
core.site_patterns), repeatedly builds a relaxed MSN over samples plus
inferred medians, tags *feasible links* (edges joining distinct
components of the threshold graph, whose components merge for pairs at
distance < threshold - epsilon), prunes everything else, drops obsolete
medians (degree < 2), and adds quasi-median vertices for triplets formed
by feasible links sharing an endpoint whose cost is within epsilon of
the global minimum.
"""

from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import networkx as nx

from ..core.alignment import Alignment
from ..core.distance import DistanceMatrix
from ..core.graph import HaplotypeNetwork
from ..core.haplotype import Haplotype
from ..core.haplotype import identify_haplotypes_from_alignment as identify_haplotypes
from ..core.sequence import Sequence
from ..core.site_patterns import condense_site_patterns, is_ambiguous
from .msn import MinimumSpanningNetwork

#: Safety cap on refinement iterations; PopART loops purely on `changed`
#: and terminates in practice, this guards a port bug from hanging.
_MAX_ITERATIONS = 10_000


class MedianJoiningNetwork(MinimumSpanningNetwork):
    """
    Construct a Median-Joining Network from haplotype data.

    Follows PopART's MedJoinNet: quasi-median vertices are inferred in
    condensed site-pattern space with site-weighted costs, and network
    pruning is driven by feasible links of a threshold graph.

    Parameters
    ----------
    distance_method : str, default='hamming'
        Method for calculating distances (MJN itself works on weighted
        Hamming distances over condensed sites, as PopART does).
    epsilon : float, default=0.0
        Bandelt's epsilon: widens both the feasible-link threshold graph
        (components merge below threshold - epsilon) and the accepted
        median cost (cost <= min cost + epsilon).
    max_median_vectors : int, optional
        PyPopART-specific cap on the number of inferred medians
        (default None = unlimited, matching PopART).
    simplify : bool, default=False
        PyPopART-specific extra (off by default): additionally smooth
        degree-2 median vertices into a single combined edge. PopART
        only ever removes medians of degree < 2.
    **kwargs : dict
        Additional parameters passed to the base class.
    """

    def __init__(
        self,
        distance_method: str = 'hamming',
        epsilon: float = 0.0,
        max_median_vectors: Optional[int] = None,
        simplify: bool = False,
        **kwargs,
    ):
        """
        Initialize MJN algorithm.

        Parameters
        ----------
        distance_method : str, default='hamming'
            Method for calculating distances.
        epsilon : float, default=0.0
            Weight parameter controlling network complexity.
        max_median_vectors : int, optional
            Optional cap on inferred medians (None = unlimited).
        simplify : bool, default=False
            Opt-in smoothing of degree-2 medians (not part of PopART).
        **kwargs : dict
            Additional parameters passed to the base class.
        """
        super().__init__(distance_method, epsilon=epsilon, **kwargs)
        self.max_median_vectors = max_median_vectors
        self.simplify = simplify

    def construct_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
        Construct MJN from sequence alignment.

        Parameters
        ----------
        alignment : Alignment
            Multiple sequence alignment.
        distance_matrix : DistanceMatrix, optional
            Ignored; MJN computes its own weighted distances in
            condensed site-pattern space.

        Returns
        -------
        HaplotypeNetwork
            Median-joining network with inferred median vertices. Median
            node sequences are expressed over condensed site patterns
            (as in PopART); sampled haplotypes keep their original
            sequences.
        """
        haplotypes = identify_haplotypes(alignment)

        if len(haplotypes) == 0:
            return HaplotypeNetwork()
        if len(haplotypes) == 1:
            network = HaplotypeNetwork()
            network.add_haplotype(haplotypes[0])
            return network

        n_samples = len(haplotypes)
        condensed, weights = condense_site_patterns([h.data for h in haplotypes])

        seqs: List[str] = list(condensed)
        labels: List[str] = [h.id for h in haplotypes]

        graph = self._compute_mjn(seqs, labels, n_samples, weights)

        # Build the HaplotypeNetwork: samples keep original sequences,
        # medians carry their condensed-space sequence.
        network = HaplotypeNetwork()
        for haplotype in haplotypes:
            if haplotype.id in graph:
                network.add_haplotype(haplotype)
        for idx, label in enumerate(labels):
            if idx >= n_samples and label in graph:
                median_seq = Sequence(
                    id=label,
                    data=seqs[idx],
                    description='Inferred median vector (condensed sites)',
                )
                network.add_haplotype(
                    Haplotype(sequence=median_seq, sample_ids=[]),
                    median_vector=True,
                )
                network.graph.nodes[label]['is_median'] = True
        for u, v, attrs in graph.edges(data=True):
            network.add_edge(u, v, distance=attrs['distance'])

        if self.simplify:
            self._smooth_degree2_medians(network)

        return network

    # ------------------------------------------------------------------
    # Core algorithm (indices into seqs/labels; medians appended at end)
    # ------------------------------------------------------------------

    def _compute_mjn(
        self,
        seqs: List[str],
        labels: List[str],
        n_samples: int,
        weights: List[int],
    ) -> nx.Graph:
        """
        Run the iterative median-joining refinement.

        Parameters
        ----------
        seqs : list of str
            Condensed sequences; medians are appended in place.
        labels : list of str
            Node labels parallel to seqs; medians are appended in place.
        n_samples : int
            Number of sampled haplotypes (prefix of seqs).
        weights : list of int
            Site weights for the condensed columns.

        Returns
        -------
        nx.Graph
            Final pruned graph over node labels with 'distance' edges.

        Raises
        ------
        RuntimeError
            If refinement fails to converge within the safety cap.
        """
        all_seq_set: Set[str] = set(seqs)
        removed: Set[int] = set()
        median_counter = 0

        for _ in range(_MAX_ITERATIONS):
            changed = False

            active = [i for i in range(len(seqs)) if i not in removed]
            matrix = self._weighted_matrix(seqs, weights, active)
            edges, feasible = self._msn_with_feasible_links(active, matrix)

            graph = nx.Graph()
            graph.add_nodes_from(active)
            for i, j, dist in edges:
                if (i, j) in feasible or (j, i) in feasible:
                    graph.add_edge(i, j, distance=dist)

            if self._remove_obsolete(graph, n_samples, seqs, all_seq_set, removed):
                changed = True

            # Pass 1: global minimum cost over quasi-medians of every
            # path v-u-w through a shared vertex u
            min_cost = float('inf')
            for u in graph.nodes:
                neighbours = list(graph.neighbors(u))
                for a in range(len(neighbours)):
                    for b in range(a):
                        v, w = neighbours[a], neighbours[b]
                        for median in self._quasi_medians(seqs[u], seqs[v], seqs[w]):
                            if median not in all_seq_set:
                                cost = self._median_cost(
                                    seqs[u], seqs[v], seqs[w], median, weights
                                )
                                min_cost = min(min_cost, cost)

            # Pass 2: add medians for feasible-link pairs sharing a vertex
            feasible_list = [(i, j) for i, j in feasible if graph.has_edge(i, j)]
            at_cap = (
                self.max_median_vectors is not None
                and sum(1 for i in range(n_samples, len(seqs)) if i not in removed)
                >= self.max_median_vectors
            )
            for e1 in range(len(feasible_list)):
                u, v = feasible_list[e1]
                for e2 in range(e1):
                    a, b = feasible_list[e2]
                    if a in (u, v):
                        w = b
                    elif b in (u, v):
                        w = a
                    else:
                        continue

                    for median in self._quasi_medians(seqs[u], seqs[v], seqs[w]):
                        if median in all_seq_set or at_cap:
                            continue
                        cost = self._median_cost(
                            seqs[u], seqs[v], seqs[w], median, weights
                        )
                        if cost <= min_cost + self.epsilon:
                            labels.append(f'Median_{median_counter}')
                            median_counter += 1
                            seqs.append(median)
                            all_seq_set.add(median)
                            changed = True
                            at_cap = (
                                self.max_median_vectors is not None
                                and sum(
                                    1
                                    for i in range(n_samples, len(seqs))
                                    if i not in removed
                                )
                                >= self.max_median_vectors
                            )

            if not changed:
                break
        else:
            raise RuntimeError('MJN refinement failed to converge')

        # Final phase (computeGraph): rebuild and prune until stable
        while True:
            active = [i for i in range(len(seqs)) if i not in removed]
            matrix = self._weighted_matrix(seqs, weights, active)
            edges, feasible = self._msn_with_feasible_links(active, matrix)

            graph = nx.Graph()
            graph.add_nodes_from(active)
            for i, j, dist in edges:
                if (i, j) in feasible or (j, i) in feasible:
                    graph.add_edge(i, j, distance=dist)

            if not self._remove_obsolete(graph, n_samples, seqs, all_seq_set, removed):
                break

        return nx.relabel_nodes(graph, {i: labels[i] for i in graph.nodes})

    def _weighted_matrix(
        self, seqs: List[str], weights: List[int], active: List[int]
    ) -> Dict[Tuple[int, int], float]:
        """
        Compute weighted pairwise distances between active sequences.

        Parameters
        ----------
        seqs : list of str
            All condensed sequences.
        weights : list of int
            Site weights.
        active : list of int
            Indices of live sequences.

        Returns
        -------
        dict
            Mapping (i, j) with i > j to weighted distance.
        """
        import numpy as np

        from ..core.site_patterns import AMBIGUOUS_CHARS

        n = len(active)
        if n == 0:
            return {}
        length = len(seqs[active[0]])
        encoded = np.zeros((n, length), dtype=np.uint8)
        for pos, idx in enumerate(active):
            encoded[pos] = np.frombuffer(seqs[idx].encode('ascii'), dtype=np.uint8)
        invalid = np.isin(
            encoded, np.frombuffer(''.join(AMBIGUOUS_CHARS).encode(), np.uint8)
        )
        valid = ~invalid
        weight_arr = np.asarray(weights, dtype=float)

        matrix: Dict[Tuple[int, int], float] = {}
        for a_pos, i in enumerate(active):
            comparable = valid[a_pos] & valid
            diffs = (encoded[a_pos] != encoded) & comparable
            row = diffs @ weight_arr
            for b_pos in range(a_pos):
                j = active[b_pos]
                matrix[(i, j)] = matrix[(j, i)] = float(row[b_pos])
        return matrix

    @staticmethod
    def _weighted_distance(s1: str, s2: str, weights: List[int]) -> float:
        """
        Site-weighted Hamming distance over condensed sequences.

        Matches HapNet::pairwiseDistance: ambiguous positions are
        skipped; a mismatch at column i contributes weights[i].

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
            Weighted distance.
        """
        total = 0
        for c1, c2, w in zip(s1, s2, weights):
            if is_ambiguous(c1) or is_ambiguous(c2):
                continue
            if c1 != c2:
                total += w
        return float(total)

    def _msn_with_feasible_links(
        self, active: List[int], matrix: Dict[Tuple[int, int], float]
    ) -> Tuple[List[Tuple[int, int, float]], Set[Tuple[int, int]]]:
        """
        Build the relaxed MSN and tag feasible links.

        Port of MedJoinNet::computeMSN: every pair at each accepted
        distance level becomes an edge; an edge is a *feasible link*
        when its endpoints lie in different components of the threshold
        graph, whose components merge for pairs at distance strictly
        below (threshold - epsilon).

        Parameters
        ----------
        active : list of int
            Indices of live sequences.
        matrix : dict
            Pairwise distances keyed by index pairs.

        Returns
        -------
        tuple of (list, set)
            All edges as (i, j, distance), and the feasible subset as
            (i, j) pairs.

        Raises
        ------
        RuntimeError
            If the pair queue empties before the graph connects.
        """
        n = len(active)
        position = {node: pos for pos, node in enumerate(active)}

        levels = defaultdict(list)
        for a_pos, i in enumerate(active):
            for j in active[:a_pos]:
                levels[matrix[(i, j)]].append((i, j))

        msn_comp = list(range(n))
        threshold_comp = list(range(n))
        ncomps = n
        max_value = float('inf')
        edges: List[Tuple[int, int, float]] = []
        feasible: Set[Tuple[int, int]] = set()

        def merge(components: List[int], i: int, j: int) -> None:
            """
            Merge two component ids with C++-style renumbering.

            Parameters
            ----------
            components : list of int
                Component id per position (modified in place).
            i : int
                First position.
            j : int
                Second position.
            """
            high = max(components[i], components[j])
            low = min(components[i], components[j])
            for k in range(n):
                if components[k] == high:
                    components[k] = low
                elif components[k] > high:
                    components[k] -= 1

        for dist in sorted(levels):
            if dist > max_value:
                break

            # Update the threshold graph for this level
            for a_pos, i in enumerate(active):
                for j in active[:a_pos]:
                    pi, pj = position[i], position[j]
                    if threshold_comp[pi] != threshold_comp[pj] and matrix[(i, j)] < (
                        dist - self.epsilon
                    ):
                        merge(threshold_comp, pi, pj)

            pairs = levels[dist]
            for i, j in pairs:
                edges.append((i, j, dist))
                if threshold_comp[position[i]] != threshold_comp[position[j]]:
                    feasible.add((i, j))

            for i, j in pairs:
                pi, pj = position[i], position[j]
                if msn_comp[pi] != msn_comp[pj]:
                    merge(msn_comp, pi, pj)
                    ncomps -= 1
                if ncomps == 1 and max_value == float('inf'):
                    max_value = dist + self.epsilon

        if ncomps > 1:
            raise RuntimeError('Pair queue empty before the graph is connected')

        return edges, feasible

    @staticmethod
    def _remove_obsolete(
        graph: nx.Graph,
        n_samples: int,
        seqs: List[str],
        all_seq_set: Set[str],
        removed: Set[int],
    ) -> bool:
        """
        Remove median vertices of degree < 2, iterating to a fixpoint.

        Port of MedJoinNet::removeObsoleteVerts (sampled haplotypes are
        never removed).

        Parameters
        ----------
        graph : nx.Graph
            Current graph over sequence indices (modified in place).
        n_samples : int
            Number of sampled haplotypes.
        seqs : list of str
            All condensed sequences.
        all_seq_set : set of str
            Live sequence strings (updated in place).
        removed : set of int
            Indices removed so far (updated in place).

        Returns
        -------
        bool
            True when any vertex was removed.
        """
        removed_any = False
        while True:
            obsolete = [
                node
                for node in graph.nodes
                if node >= n_samples and graph.degree(node) < 2
            ]
            if not obsolete:
                return removed_any
            for node in obsolete:
                graph.remove_node(node)
                all_seq_set.discard(seqs[node])
                removed.add(node)
                removed_any = True

    @staticmethod
    def _quasi_medians(seq_a: str, seq_b: str, seq_c: str) -> Set[str]:
        """
        Compute the quasi-median sequences of a triplet.

        Port of MedJoinNet::computeQuasiMedianSeqs: per position, take
        the majority character; positions where all three differ become
        stars, which are resolved recursively three ways.

        Parameters
        ----------
        seq_a : str
            First sequence.
        seq_b : str
            Second sequence.
        seq_c : str
            Third sequence.

        Returns
        -------
        set of str
            All quasi-median sequences of the triplet.
        """
        qm = []
        has_star = False
        for a, b, c in zip(seq_a, seq_b, seq_c):
            if a == b or a == c:
                qm.append(a)
            elif b == c:
                qm.append(b)
            else:
                qm.append('*')
                has_star = True

        seed = ''.join(qm)
        if not has_star:
            return {seed}

        medians: Set[str] = set()
        stack = [seed]
        while stack:
            seq = stack.pop()
            star = seq.find('*')
            resolved = [
                seq[:star] + source[star] + seq[star + 1 :]
                for source in (seq_a, seq_b, seq_c)
            ]
            if '*' in resolved[0]:
                stack.extend(resolved)
            else:
                medians.update(resolved)
        return medians

    def _median_cost(
        self, seq_u: str, seq_v: str, seq_w: str, median: str, weights: List[int]
    ) -> float:
        """
        Cost of a candidate median (MedJoinNet::computeCost).

        Parameters
        ----------
        seq_u : str
            First triplet sequence.
        seq_v : str
            Second triplet sequence.
        seq_w : str
            Third triplet sequence.
        median : str
            Candidate median sequence.
        weights : list of int
            Site weights.

        Returns
        -------
        float
            Sum of weighted distances from the median to the triplet.
        """
        return (
            self._weighted_distance(seq_u, median, weights)
            + self._weighted_distance(seq_v, median, weights)
            + self._weighted_distance(seq_w, median, weights)
        )

    def _smooth_degree2_medians(self, network: HaplotypeNetwork) -> None:
        """
        Opt-in smoothing: collapse degree-2 medians into single edges.

        Not part of PopART's MJN (removeObsoleteVerts only drops medians
        of degree < 2); enabled via simplify=True.

        Parameters
        ----------
        network : HaplotypeNetwork
            Network to smooth in place.
        """
        graph = network.graph
        changed = True
        while changed:
            changed = False
            for node in list(graph.nodes):
                if not graph.nodes[node].get('is_median'):
                    continue
                if graph.degree(node) != 2:
                    continue
                n1, n2 = graph.neighbors(node)
                combined = graph[node][n1].get('distance', 1) + graph[node][n2].get(
                    'distance', 1
                )
                network.remove_haplotype(node)
                if not graph.has_edge(n1, n2):
                    network.add_edge(n1, n2, distance=combined)
                changed = True
                break

    def get_parameters(self) -> dict:
        """
        Get algorithm parameters.

        Returns
        -------
        dict
            Parameters including epsilon, max_median_vectors, simplify.
        """
        params = super().get_parameters()
        params['max_median_vectors'] = self.max_median_vectors
        params['simplify'] = self.simplify
        return params
