"""
Parsimony trees and ancestral-state sampling for the Parsimony Network.

Ports the relevant behaviour of PopART's tree/ParsimonyTree.cpp with
uniform substitution costs: a Fitch up-pass computes the parsimony
score and per-site state sets, and a randomised down-pass samples one
concrete ancestral sequence per internal node (PopART's
computeAncestors resolves ties randomly on every call, which is what
makes repeated edge sampling explore the ancestor space).

State sets are numpy uint8 bitmasks (A=1, C=2, G=4, T=8; gaps and
ambiguity codes are the full mask, mirroring SankoffUp's zero-cost '-'
handling), so the up-pass is vectorised over sites.

PopART receives its parsimony trees from a Nexus TREES block; PyPopART
has no tree input path yet, so trees are generated here by random-order
stepwise addition, keeping each taxon's best-scoring insertion point.
"""

import random
from typing import Dict, List, Optional, Tuple
from typing import Sequence as SequenceType

import networkx as nx
import numpy as np

_BIT_OF = {'A': 1, 'C': 2, 'G': 4, 'T': 8, 'U': 8}
_ALL_BITS = 15

#: Byte-indexed lookup table: character -> state bitmask.
_ENCODE_TABLE = np.full(256, _ALL_BITS, dtype=np.uint8)
for _char, _bit in _BIT_OF.items():
    _ENCODE_TABLE[ord(_char)] = _bit
    _ENCODE_TABLE[ord(_char.lower())] = _bit

_CHAR_OF_BIT = {1: 'A', 2: 'C', 4: 'G', 8: 'T'}


def encode_sequence(sequence: str) -> np.ndarray:
    """
    Encode a sequence as per-site Fitch state bitmasks.

    Parameters
    ----------
    sequence : str
        Sequence string.

    Returns
    -------
    np.ndarray
        Uint8 array of state bitmasks (gaps/ambiguity = all states).
    """
    raw = np.frombuffer(sequence.encode('ascii'), dtype=np.uint8)
    return _ENCODE_TABLE[raw]


class ParsimonyTree:
    """
    An unrooted tree over haplotype indices with Fitch machinery.

    Parameters
    ----------
    topology : nx.Graph
        Tree with leaf nodes 0..n-1 (haplotype indices) and internal
        nodes labelled with negative integers.
    sequences : list of str
        Condensed sequence per leaf index.
    weights : list of int
        Site weights for the condensed columns.
    encoded : list of np.ndarray, optional
        Precomputed encode_sequence arrays per leaf (shared across the
        candidate trees stepwise addition evaluates).
    """

    def __init__(
        self,
        topology: nx.Graph,
        sequences: SequenceType[str],
        weights: SequenceType[int],
        encoded: Optional[List[np.ndarray]] = None,
    ):
        """
        Initialize the tree.

        Parameters
        ----------
        topology : nx.Graph
            Tree topology (leaves are haplotype indices >= 0).
        sequences : list of str
            Condensed leaf sequences.
        weights : list of int
            Site weights.
        encoded : list of np.ndarray, optional
            Precomputed leaf encodings.
        """
        self.topology = topology
        self.sequences = list(sequences)
        self.weights = np.asarray(weights, dtype=float)
        self._encoded = (
            encoded
            if encoded is not None
            else [encode_sequence(s) for s in self.sequences]
        )
        self._state_sets: Dict[int, np.ndarray] = {}
        self._score: float = 0.0
        self._up_pass_done = False

    def compute_score(self) -> float:
        """
        Run the Fitch up-pass and return the weighted parsimony score.

        Returns
        -------
        float
            Weighted number of state changes implied by the tree.
        """
        root = self._pick_root()
        order = list(nx.dfs_postorder_nodes(self.topology, source=root))
        parent = nx.dfs_predecessors(self.topology, source=root)

        self._state_sets = {}
        score = 0.0
        for node in order:
            if node >= 0:  # leaf
                self._state_sets[node] = self._encoded[node]
                continue
            children = [
                neighbour
                for neighbour in self.topology.neighbors(node)
                if parent.get(neighbour) == node
            ]
            merged = self._state_sets[children[0]]
            for child in children[1:]:
                child_set = self._state_sets[child]
                intersection = merged & child_set
                empty = intersection == 0
                if empty.any():
                    score += float(self.weights[empty].sum())
                    merged = np.where(empty, merged | child_set, intersection)
                else:
                    merged = intersection
            self._state_sets[node] = merged

        self._root = root
        self._parent = parent
        self._score = score
        self._up_pass_done = True
        return score

    def sample_ancestors(self, rng: random.Random) -> Dict[int, str]:
        """
        Sample one concrete sequence per internal node (random Fitch).

        Ties are resolved randomly on every call, as in PopART's
        computeAncestors, so repeated calls explore alternative equally
        parsimonious ancestral reconstructions.

        Parameters
        ----------
        rng : random.Random
            Random source.

        Returns
        -------
        dict
            Mapping internal node id to its sampled sequence.
        """
        if not self._up_pass_done:
            self.compute_score()

        ancestors: Dict[int, str] = {}
        assigned_bits: Dict[int, np.ndarray] = {}
        for node in nx.dfs_preorder_nodes(self.topology, source=self._root):
            if node >= 0:
                continue
            states = self._state_sets[node]
            parent = self._parent.get(node)

            if parent is not None and parent < 0:
                parent_bits = assigned_bits[parent]
                chosen = np.where(parent_bits & states, parent_bits, 0).astype(np.uint8)
            else:
                chosen = np.zeros(len(states), dtype=np.uint8)

            # Random resolution wherever the parent state isn't allowed
            for site in np.flatnonzero(chosen == 0):
                bits = int(states[site])
                options = [b for b in (1, 2, 4, 8) if bits & b]
                chosen[site] = rng.choice(options)

            assigned_bits[node] = chosen
            ancestors[node] = ''.join(_CHAR_OF_BIT[int(b)] for b in chosen)
        return ancestors

    def edge_sequences(self, ancestors: Dict[int, str]) -> List[Tuple[str, str]]:
        """
        List tree edges as (sequence, sequence) pairs.

        Parameters
        ----------
        ancestors : dict
            Sampled internal-node sequences from sample_ancestors.

        Returns
        -------
        list of tuple
            One (seq_from, seq_to) pair per tree edge.
        """

        def seq_of(node: int) -> str:
            """
            Return the sequence string for a leaf or sampled ancestor.

            Parameters
            ----------
            node : int
                Tree node id.

            Returns
            -------
            str
                The node's sequence.
            """
            return self.sequences[node] if node >= 0 else ancestors[node]

        return [(seq_of(u), seq_of(v)) for u, v in self.topology.edges()]

    def _pick_root(self) -> int:
        """
        Choose a root node for traversals (an internal node if any).

        Returns
        -------
        int
            Node id to root traversals at.
        """
        for node in self.topology.nodes:
            if node < 0:
                return node
        return next(iter(self.topology.nodes))


def stepwise_addition_tree(
    sequences: SequenceType[str],
    weights: SequenceType[int],
    rng: random.Random,
    encoded: Optional[List[np.ndarray]] = None,
) -> ParsimonyTree:
    """
    Build a parsimony tree by random-order stepwise addition.

    Taxa are inserted in random order; each is grafted onto the edge
    that minimises the Fitch score, breaking ties randomly. This is the
    classic stepwise-addition heuristic and serves as PyPopART's tree
    source in place of PopART's user-supplied Nexus trees.

    Parameters
    ----------
    sequences : list of str
        Condensed sequence per haplotype index.
    weights : list of int
        Site weights.
    rng : random.Random
        Random source.
    encoded : list of np.ndarray, optional
        Precomputed leaf encodings (shared across candidate trees).

    Returns
    -------
    ParsimonyTree
        The constructed tree.
    """
    n = len(sequences)
    if encoded is None:
        encoded = [encode_sequence(s) for s in sequences]
    order = list(range(n))
    rng.shuffle(order)

    topology = nx.Graph()
    internal_counter = [0]

    def new_internal() -> int:
        """
        Allocate the next internal-node id.

        Returns
        -------
        int
            A fresh negative node id.
        """
        internal_counter[0] -= 1
        return internal_counter[0]

    if n == 1:
        topology.add_node(order[0])
        return ParsimonyTree(topology, sequences, weights, encoded)
    if n == 2:
        topology.add_edge(order[0], order[1])
        return ParsimonyTree(topology, sequences, weights, encoded)

    # Seed with the first three taxa around one internal node
    hub = new_internal()
    for leaf in order[:3]:
        topology.add_edge(hub, leaf)

    for leaf in order[3:]:
        best_score = None
        best_edges: List[Tuple[int, int]] = []
        for u, v in list(topology.edges()):
            attach = new_internal()
            topology.remove_edge(u, v)
            topology.add_edge(u, attach)
            topology.add_edge(attach, v)
            topology.add_edge(attach, leaf)

            score = ParsimonyTree(topology, sequences, weights, encoded).compute_score()
            if best_score is None or score < best_score:
                best_score = score
                best_edges = [(u, v)]
            elif score == best_score:
                best_edges.append((u, v))

            topology.remove_node(attach)
            topology.add_edge(u, v)

        u, v = rng.choice(best_edges)
        attach = new_internal()
        topology.remove_edge(u, v)
        topology.add_edge(u, attach)
        topology.add_edge(attach, v)
        topology.add_edge(attach, leaf)

    return ParsimonyTree(topology, sequences, weights, encoded)


def sample_parsimony_trees(
    sequences: SequenceType[str],
    weights: SequenceType[int],
    n_trees: int,
    rng: random.Random,
) -> List[ParsimonyTree]:
    """
    Generate a set of stepwise-addition parsimony trees.

    Parameters
    ----------
    sequences : list of str
        Condensed sequence per haplotype index.
    weights : list of int
        Site weights.
    n_trees : int
        Number of trees to build.
    rng : random.Random
        Random source.

    Returns
    -------
    list of ParsimonyTree
        The sampled trees (scores already computed).
    """
    encoded = [encode_sequence(s) for s in sequences]
    trees = []
    for _ in range(n_trees):
        tree = stepwise_addition_tree(sequences, weights, rng, encoded)
        tree.compute_score()
        trees.append(tree)
    return trees


__all__ = [
    'ParsimonyTree',
    'encode_sequence',
    'sample_parsimony_trees',
    'stepwise_addition_tree',
]
