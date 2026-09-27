"""Unit tests for the Parsimony Network (PN) algorithm."""

import pytest

from pypopart.algorithms.parsimony_net import ParsimonyNetwork
from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence


def small_alignment():
    """Build a five-sequence alignment with simple structure."""
    return Alignment(
        [
            Sequence('seq1', 'AAAAAA'),
            Sequence('seq2', 'AAAAAT'),
            Sequence('seq3', 'AAAATT'),
            Sequence('seq4', 'AAATTT'),
            Sequence('seq5', 'AATTTT'),
        ]
    )


class TestParsimonyNetwork:
    """Test cases for the PN algorithm."""

    def test_pn_initialization(self):
        """Defaults follow the AncestralSeqNet port."""
        pn = ParsimonyNetwork()
        assert pn.distance_method == 'hamming'
        assert pn.n_trees == 20
        assert pn.alpha == 0.95
        assert pn.n_iterations is None
        assert pn.random_seed is None

    def test_pn_custom_parameters(self):
        """Custom parameters are stored and reported."""
        pn = ParsimonyNetwork(n_trees=5, alpha=0.5, n_iterations=100, random_seed=42)
        params = pn.get_parameters()
        assert params['n_trees'] == 5
        assert params['alpha'] == 0.5
        assert params['n_iterations'] == 100
        assert params['random_seed'] == 42

    def test_pn_rejects_removed_parameter(self):
        """The old min_edge_frequency parameter fails loudly."""
        with pytest.raises(TypeError, match='min_edge_frequency'):
            ParsimonyNetwork(min_edge_frequency=0.05)

    def test_pn_empty_alignment(self):
        """Empty alignment gives an empty network."""
        network = ParsimonyNetwork().construct_network(Alignment())
        assert len(network) == 0

    def test_pn_single_sequence(self):
        """Single sequence gives a single-node network."""
        network = ParsimonyNetwork().construct_network(
            Alignment([Sequence('seq1', 'ATCG')])
        )
        assert len(network) == 1
        assert len(network.edges) == 0

    def test_pn_two_identical_sequences(self):
        """Identical sequences collapse to one haplotype."""
        network = ParsimonyNetwork(n_trees=3, random_seed=1).construct_network(
            Alignment([Sequence('seq1', 'ATCG'), Sequence('seq2', 'ATCG')])
        )
        assert len(network) == 1

    def test_pn_two_different_sequences(self):
        """Two haplotypes connect with one edge."""
        network = ParsimonyNetwork(n_trees=3, random_seed=1).construct_network(
            Alignment([Sequence('seq1', 'ATCG'), Sequence('seq2', 'ATCC')])
        )
        assert len(network) == 2
        assert network.is_connected()

    def test_pn_sampled_haplotypes_always_present(self):
        """Every sampled haplotype appears in the network."""
        network = ParsimonyNetwork(n_trees=5, random_seed=7).construct_network(
            small_alignment()
        )
        for hap_id in ('H1', 'H2', 'H3', 'H4', 'H5'):
            assert network.has_node(hap_id)

    def test_pn_samples_connected(self):
        """Sampled haplotypes end up mutually reachable."""
        network = ParsimonyNetwork(n_trees=5, random_seed=7).construct_network(
            small_alignment()
        )
        assert network.is_connected()

    def test_pn_reproducibility_with_seed(self):
        """The same seed reproduces the same network."""
        net1 = ParsimonyNetwork(n_trees=5, random_seed=42).construct_network(
            small_alignment()
        )
        net2 = ParsimonyNetwork(n_trees=5, random_seed=42).construct_network(
            small_alignment()
        )
        assert set(net1.graph.nodes) == set(net2.graph.nodes)
        assert set(map(frozenset, net1.graph.edges)) == set(
            map(frozenset, net2.graph.edges)
        )

    def test_pn_ancestors_flagged_as_medians(self):
        """Inferred ancestral vertices carry median markers."""
        network = ParsimonyNetwork(n_trees=5, random_seed=3).construct_network(
            small_alignment()
        )
        for node, attrs in network.graph.nodes(data=True):
            if str(node).startswith('Median_'):
                assert attrs.get('median_vector') is True
                # Ancestors keep degree > 1 (isolated ones are dropped)
                assert network.get_degree(node) >= 2

    def test_pn_alpha_zero_keeps_all_edges(self):
        """alpha=0 disables frequency pruning entirely."""
        pruned = ParsimonyNetwork(
            n_trees=5, random_seed=5, alpha=0.95
        ).construct_network(small_alignment())
        unpruned = ParsimonyNetwork(
            n_trees=5, random_seed=5, alpha=0.0
        ).construct_network(small_alignment())
        assert len(unpruned.edges) >= len(pruned.edges)

    def test_pn_with_gaps(self):
        """Gapped sequences are handled."""
        network = ParsimonyNetwork(n_trees=3, random_seed=1).construct_network(
            Alignment([Sequence('seq1', 'AT-G'), Sequence('seq2', 'ATCG')])
        )
        assert len(network) >= 1


class TestPNEdgeDistances:
    """Edge distances come from the shared sequence_distance path."""

    def test_pn_edge_distance(self):
        """One mismatch gives distance 1."""
        from pypopart.core.distance import sequence_distance

        assert sequence_distance(Sequence('a', 'ATCG'), Sequence('b', 'ATCT')) == 1.0

    def test_pn_edge_distance_identical(self):
        """Identical sequences have zero distance."""
        from pypopart.core.distance import sequence_distance

        assert sequence_distance(Sequence('a', 'ATCG'), Sequence('b', 'ATCG')) == 0.0

    def test_pn_edge_distance_unequal_length_raises(self):
        """Unequal-length sequences raise, matching PopART's behaviour."""
        from pypopart.core.distance import sequence_distance

        with pytest.raises(ValueError, match='same length'):
            sequence_distance(Sequence('a', 'ATCG'), Sequence('b', 'ATC'))


class TestParsimonyTrees:
    """The stepwise-addition tree generator and Fitch machinery."""

    def test_stepwise_tree_contains_all_leaves(self):
        """Every haplotype index appears as a leaf."""
        import random

        from pypopart.algorithms.ancestral import stepwise_addition_tree

        seqs = ['AAA', 'AAT', 'ATT', 'TTT']
        tree = stepwise_addition_tree(seqs, [1, 1, 1], random.Random(1))
        leaves = [n for n in tree.topology.nodes if n >= 0]
        assert sorted(leaves) == [0, 1, 2, 3]

    def test_fitch_score_simple(self):
        """A three-leaf star scores weighted state-set unions."""
        import networkx as nx

        from pypopart.algorithms.ancestral import ParsimonyTree

        topology = nx.Graph()
        for leaf in (0, 1, 2):
            topology.add_edge(-1, leaf)
        # Column 0: all A (0 changes). Column 1 (weight 2): A/A/T -> one
        # union event costing 2. Column 2: T/A/A -> one event costing 1.
        tree = ParsimonyTree(topology, ['AAT', 'AAA', 'ATA'], [1, 2, 1])
        assert tree.compute_score() == 3

    def test_fitch_ancestor_sampling(self):
        """Internal nodes get concrete sequences from their state sets."""
        import random

        import networkx as nx

        from pypopart.algorithms.ancestral import ParsimonyTree

        topology = nx.Graph()
        for leaf in (0, 1, 2):
            topology.add_edge(-1, leaf)
        tree = ParsimonyTree(topology, ['AAT', 'AAA', 'ATA'], [1, 1, 1])
        tree.compute_score()
        ancestors = tree.sample_ancestors(random.Random(1))
        assert set(ancestors) == {-1}
        assert len(ancestors[-1]) == 3
        # First column is A in all leaves; the ancestor must agree
        assert ancestors[-1][0] == 'A'
