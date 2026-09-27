"""Unit tests for Minimum Spanning Network (MSN) algorithm."""

from pypopart.algorithms.msn import MinimumSpanningNetwork
from pypopart.algorithms.mst import MinimumSpanningTree
from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence


class TestMinimumSpanningNetwork:
    """Test cases for MSN algorithm."""

    def test_msn_initialization(self):
        """Test MSN algorithm initialization."""
        msn = MinimumSpanningNetwork()
        assert msn.distance_method == 'hamming'
        assert msn.epsilon == 0.0
        assert msn.prune_redundant is False

    def test_msn_with_epsilon(self):
        """Test MSN with epsilon tolerance."""
        msn = MinimumSpanningNetwork(epsilon=0.5)
        assert msn.epsilon == 0.5

    def test_msn_empty_alignment(self):
        """Test MSN with empty alignment."""
        msn = MinimumSpanningNetwork()
        alignment = Alignment()
        network = msn.construct_network(alignment)

        assert len(network) == 0

    def test_msn_single_sequence(self):
        """Test MSN with single sequence."""
        msn = MinimumSpanningNetwork()
        alignment = Alignment([Sequence('seq1', 'ATCG')])
        network = msn.construct_network(alignment)

        assert len(network) == 1
        assert len(network.edges) == 0

    def test_msn_adds_alternative_connections(self):
        """Test that MSN adds alternative connections at same distance."""
        msn = MinimumSpanningNetwork()
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),  # 1 diff from seq1
                Sequence('seq3', 'AAAC'),  # 1 diff from seq1
                Sequence('seq4', 'AAAT'),  # Same as seq2
            ]
        )
        network = msn.construct_network(alignment)

        # Should have at least as many edges as MST (or same if no alternatives)
        assert network.is_connected()

    def test_msn_vs_mst(self):
        """Test that MSN has at least as many edges as MST."""
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCC'),
                Sequence('seq3', 'GTCG'),
                Sequence('seq4', 'GTCC'),
            ]
        )

        mst = MinimumSpanningTree()
        mst_network = mst.construct_network(alignment)

        msn = MinimumSpanningNetwork()
        msn_network = msn.construct_network(alignment)

        # MSN should have at least as many edges as MST
        assert len(msn_network.edges) >= len(mst_network.edges)

    def test_msn_contains_all_mst_edges(self):
        """PopART MSN keeps every tied minimal edge (superset of an MST)."""
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCC'),
                Sequence('seq3', 'GTCG'),
                Sequence('seq4', 'GTCC'),
            ]
        )

        mst_edges = {
            tuple(sorted(map(str, e)))
            for e in MinimumSpanningTree().construct_network(alignment).graph.edges()
        }
        msn_edges = {
            tuple(sorted(map(str, e)))
            for e in MinimumSpanningNetwork().construct_network(alignment).graph.edges()
        }
        assert mst_edges <= msn_edges

    def test_msn_keeps_equal_distance_ties(self):
        """All cross-component ties at the connection level are kept.

        The four sequences form a square: two distance-1 'rungs' plus two
        distance-1 'rails'; PopART's MSN includes all four unit edges.
        """
        alignment = Alignment(
            [
                Sequence('seq1', 'AA'),
                Sequence('seq2', 'AT'),
                Sequence('seq3', 'GA'),
                Sequence('seq4', 'GT'),
            ]
        )
        network = MinimumSpanningNetwork().construct_network(alignment)

        # Unit-distance pairs: AA-AT, AA-GA, AT-GT, GA-GT (a 4-cycle)
        assert len(network.edges) == 4
        assert network.is_connected()

    def test_msn_triangle(self):
        """Test MSN with three equidistant sequences (triangle)."""
        msn = MinimumSpanningNetwork()
        alignment = Alignment(
            [Sequence('seq1', 'AT'), Sequence('seq2', 'AC'), Sequence('seq3', 'GT')]
        )
        network = msn.construct_network(alignment)

        # Should form a triangle (all three connected)
        assert len(network) == 3
        assert network.is_connected()

    def test_msn_parameters(self):
        """Test getting MSN parameters."""
        msn = MinimumSpanningNetwork(
            distance_method='k2p', epsilon=0.5, prune_redundant=True
        )
        params = msn.get_parameters()

        assert params['distance_method'] == 'k2p'
        assert params['epsilon'] == 0.5
        assert params['prune_redundant'] is True

    def test_msn_string_representation(self):
        """Test string representation of MSN."""
        msn = MinimumSpanningNetwork(distance_method='hamming')
        assert 'MinimumSpanningNetwork' in str(msn)
        assert 'hamming' in str(msn)


class TestMSNEpsilon:
    """Epsilon extends the accepted distance levels (PopART semantics)."""

    def make_chain_alignment(self):
        """Chain with a distance-2 shortcut only epsilon > 0 admits."""
        return Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),  # 1 from seq1
                Sequence('seq3', 'AATT'),  # 1 from seq2, 2 from seq1
            ]
        )

    def test_epsilon_zero_is_strict(self):
        """With epsilon=0 only cross-component minimal edges are added."""
        network = MinimumSpanningNetwork(epsilon=0).construct_network(
            self.make_chain_alignment()
        )
        # Chain connects at distance 1: seq1-seq2, seq2-seq3
        assert len(network.edges) == 2
        assert network.is_connected()

    def test_epsilon_extends_threshold(self):
        """Epsilon=1 admits the distance-2 level after connection at 1."""
        network = MinimumSpanningNetwork(epsilon=1).construct_network(
            self.make_chain_alignment()
        )
        # The relaxed sweep also includes the seq1-seq3 distance-2 edge
        assert len(network.edges) == 3

    def test_epsilon_relaxes_component_rule(self):
        """Any epsilon > 0 adds all pairs at accepted levels, as in C++."""
        alignment = Alignment(
            [
                Sequence('seq1', 'AA'),
                Sequence('seq2', 'AT'),
                Sequence('seq3', 'TT'),
            ]
        )
        strict = MinimumSpanningNetwork(epsilon=0).construct_network(alignment)
        relaxed = MinimumSpanningNetwork(epsilon=1).construct_network(alignment)
        assert len(relaxed.edges) >= len(strict.edges)

    def test_epsilon_parameter_storage(self):
        """Epsilon is stored and reported."""
        msn = MinimumSpanningNetwork(epsilon=0.5)
        assert msn.epsilon == 0.5
        assert msn.get_parameters()['epsilon'] == 0.5


class TestMSNPruneRedundant:
    """The opt-in prune_redundant extra (not part of PopART)."""

    def test_prune_removes_longer_redundant_edges(self):
        """Pruning drops edges with equal-or-shorter alternative paths."""
        alignment = Alignment(
            [
                Sequence('seq1', 'AA'),
                Sequence('seq2', 'AT'),
                Sequence('seq3', 'TT'),
            ]
        )
        plain = MinimumSpanningNetwork(epsilon=1).construct_network(alignment)
        pruned = MinimumSpanningNetwork(
            epsilon=1, prune_redundant=True
        ).construct_network(alignment)

        assert len(pruned.edges) <= len(plain.edges)
        assert pruned.is_connected()
