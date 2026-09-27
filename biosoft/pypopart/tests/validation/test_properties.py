"""
Property-based invariants the algorithms must satisfy.

Properties marked xfail encode C++ PopART semantics the current Python
implementations do not yet meet; the Phase 4 parity rewrites turn them
green (and the xfail markers must then be removed).
"""

import random

import pytest

from pypopart.algorithms import ALGORITHMS, build
from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence


def make_alignment(seed: int = 5, n_seqs: int = 10, length: int = 30) -> Alignment:
    """
    Build a reproducible alignment with clustered variation.

    Parameters
    ----------
    seed : int, default=5
        RNG seed.
    n_seqs : int, default=10
        Number of sequences.
    length : int, default=30
        Alignment length.

    Returns
    -------
    Alignment
        Synthetic alignment without ambiguity codes.
    """
    rng = random.Random(seed)
    ancestor = [rng.choice('ACGT') for _ in range(length)]
    sequences = []
    for i in range(n_seqs):
        seq = ancestor.copy()
        for _ in range(rng.randrange(1, 4)):
            pos = rng.randrange(length)
            seq[pos] = rng.choice('ACGT')
        sequences.append(Sequence(f's{i}', ''.join(seq)))
    return Alignment(sequences)


ALIGNMENT = make_alignment()


@pytest.mark.parametrize('name', sorted(ALGORITHMS))
def test_sampled_haplotypes_present(name):
    """Every sampled haplotype appears as a node in the network."""
    from pypopart.core.haplotype import identify_haplotypes_from_alignment

    params = {'random_seed': 42} if name == 'pn' else {}
    network = build(name, **params).build_network(ALIGNMENT)

    haplotype_ids = {h.id for h in identify_haplotypes_from_alignment(ALIGNMENT)}
    node_ids = {str(n) for n in network.graph.nodes()}
    assert haplotype_ids <= node_ids


@pytest.mark.parametrize('name', ['mst', 'msn', 'mjn', 'tsw'])
def test_network_connected(name):
    """Spanning-style networks connect all haplotypes."""
    network = build(name).build_network(ALIGNMENT)
    assert network.is_connected()


def test_tcs_connected():
    """PopART's TCS always produces a fully connected network."""
    network = build('tcs').build_network(ALIGNMENT)
    assert network.is_connected()


def test_msn_contains_mst_edges():
    """The MSN edge set is a superset of some minimum spanning tree."""
    mst_net = build('mst').build_network(ALIGNMENT)
    msn_net = build('msn').build_network(ALIGNMENT)

    mst_edges = {tuple(sorted(map(str, e))) for e in mst_net.graph.edges()}
    msn_edges = {tuple(sorted(map(str, e))) for e in msn_net.graph.edges()}
    assert mst_edges <= msn_edges


def test_mst_total_length_minimal():
    """No other spanning tree has smaller total distance than the MST."""
    import networkx as nx

    from pypopart.core.distance import pairwise_distance_matrix
    from pypopart.core.haplotype import identify_haplotypes_from_alignment

    network = build('mst').build_network(ALIGNMENT)
    mst_total = sum(
        attrs.get('distance', 0) for _, _, attrs in network.graph.edges(data=True)
    )

    haplotypes = identify_haplotypes_from_alignment(ALIGNMENT)
    matrix = pairwise_distance_matrix(haplotypes)
    complete = nx.Graph()
    for i, h1 in enumerate(haplotypes):
        for h2 in haplotypes[i + 1 :]:
            complete.add_edge(h1.id, h2.id, weight=matrix.get_distance(h1.id, h2.id))
    reference = nx.minimum_spanning_tree(complete)
    reference_total = sum(attrs['weight'] for _, _, attrs in reference.edges(data=True))
    assert mst_total == reference_total
