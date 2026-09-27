"""Unit tests for the Median-Joining Network (MJN) algorithm.

Consolidates the former test_mjn_epsilon_fix.py and
test_mjn_edge_not_found_fix.py regression files.
"""

from pypopart.algorithms.mjn import MedianJoiningNetwork
from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence


def median_nodes(network):
    """Return the ids of median-vector nodes in a network."""
    return [
        node
        for node, attrs in network.graph.nodes(data=True)
        if attrs.get('is_median') or attrs.get('median_vector')
    ]


class TestMedianJoiningNetwork:
    """Test cases for MJN algorithm."""

    def test_mjn_initialization(self):
        """Test MJN algorithm initialization and PopART-parity defaults."""
        mjn = MedianJoiningNetwork()
        assert mjn.distance_method == 'hamming'
        assert mjn.epsilon == 0.0
        # simplify (degree-2 smoothing) is a PyPopART extra, off by default
        assert mjn.simplify is False
        assert mjn.max_median_vectors is None

    def test_mjn_simplify_opt_in(self):
        """Test MJN with opt-in simplification."""
        mjn = MedianJoiningNetwork(simplify=True)
        assert mjn.simplify is True

    def test_mjn_empty_alignment(self):
        """Test MJN with empty alignment."""
        network = MedianJoiningNetwork().construct_network(Alignment())
        assert len(network) == 0

    def test_mjn_single_sequence(self):
        """Test MJN with single sequence."""
        network = MedianJoiningNetwork().construct_network(
            Alignment([Sequence('seq1', 'ATCG')])
        )
        assert len(network) == 1
        assert len(network.edges) == 0

    def test_mjn_two_sequences(self):
        """Test MJN with two sequences (no median vectors)."""
        network = MedianJoiningNetwork().construct_network(
            Alignment([Sequence('seq1', 'ATCG'), Sequence('seq2', 'ATCC')])
        )
        assert len(network) == 2
        assert network.is_connected()

    def test_mjn_connected(self):
        """MJN always yields a connected network."""
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AATT'),
                Sequence('seq3', 'TTAA'),
                Sequence('seq4', 'TTTT'),
            ]
        )
        network = MedianJoiningNetwork().construct_network(alignment)
        assert network.is_connected()
        # All four sampled haplotypes present
        for hap_id in ('H1', 'H2', 'H3', 'H4'):
            assert network.has_node(hap_id)

    def test_mjn_median_nodes_flagged(self):
        """Inferred medians carry the median markers."""
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAAAA'),
                Sequence('seq2', 'AATTTT'),
                Sequence('seq3', 'TTAATT'),
                Sequence('seq4', 'TTTTAA'),
            ]
        )
        network = MedianJoiningNetwork(epsilon=2).construct_network(alignment)
        for node in median_nodes(network):
            assert node.startswith('Median_')
            assert network.graph.nodes[node]['median_vector'] is True

    def test_mjn_max_median_vectors_cap(self):
        """The opt-in cap limits the number of inferred medians."""
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAAAA'),
                Sequence('seq2', 'AATTTT'),
                Sequence('seq3', 'TTAATT'),
                Sequence('seq4', 'TTTTAA'),
            ]
        )
        network = MedianJoiningNetwork(
            epsilon=2, max_median_vectors=1
        ).construct_network(alignment)
        assert len(median_nodes(network)) <= 1


class TestMJNQuasiMedians:
    """The quasi-median machinery ported from computeQuasiMedianSeqs."""

    def test_majority_positions(self):
        """Majority positions resolve without stars."""
        medians = MedianJoiningNetwork._quasi_medians('AAT', 'AAA', 'ATA')
        assert medians == {'AAA'}

    def test_star_positions_branch_three_ways(self):
        """All-different positions yield three resolutions."""
        medians = MedianJoiningNetwork._quasi_medians('A', 'C', 'G')
        assert medians == {'A', 'C', 'G'}

    def test_multiple_stars_expand_recursively(self):
        """Two star positions give up to nine resolutions."""
        medians = MedianJoiningNetwork._quasi_medians('AA', 'CC', 'GG')
        assert medians == {a + b for a in 'ACG' for b in 'ACG'}

    def test_weighted_cost(self):
        """Median cost is the weighted distance sum to the triplet."""
        mjn = MedianJoiningNetwork()
        cost = mjn._median_cost('AA', 'AT', 'TA', 'AA', [1, 2])
        # d(AA,AA)=0, d(AT,AA)=2 (weight of col 1), d(TA,AA)=1
        assert cost == 3.0


class TestMJNEpsilon:
    """Epsilon drives both the threshold graph and cost acceptance."""

    def make_alignment(self):
        """Four-haplotype alignment that can host medians."""
        return Alignment(
            [
                Sequence('seq1', 'AAAAAA'),
                Sequence('seq2', 'AATTTT'),
                Sequence('seq3', 'TTAATT'),
                Sequence('seq4', 'TTTTAA'),
            ]
        )

    def test_epsilon_zero_baseline(self):
        """Epsilon 0 builds a connected network."""
        network = MedianJoiningNetwork(epsilon=0.0).construct_network(
            self.make_alignment()
        )
        assert network.is_connected()

    def test_larger_epsilon_never_smaller_network(self):
        """Growing epsilon admits at least as many medians."""
        small = MedianJoiningNetwork(epsilon=0.0).construct_network(
            self.make_alignment()
        )
        large = MedianJoiningNetwork(epsilon=2.0).construct_network(
            self.make_alignment()
        )
        assert len(large) >= len(small)

    def test_median_ids_distinct_from_haplotype_ids(self):
        """Median ids never collide with H-number haplotype ids."""
        network = MedianJoiningNetwork(epsilon=2.0).construct_network(
            self.make_alignment()
        )
        sampled = {n for n in network.nodes if str(n).startswith('H')}
        medians = set(median_nodes(network))
        assert sampled.isdisjoint(medians)


class TestMJNSimplify:
    """The opt-in degree-2 smoothing extra."""

    def make_alignment(self):
        """Alignment prone to producing chained medians."""
        return Alignment(
            [
                Sequence('seq1', 'AAAAAAAA'),
                Sequence('seq2', 'AAAATTTT'),
                Sequence('seq3', 'TTTTAAAA'),
                Sequence('seq4', 'TTTTTTTT'),
            ]
        )

    def test_simplify_never_larger(self):
        """Smoothing cannot increase the node count."""
        plain = MedianJoiningNetwork(epsilon=2.0).construct_network(
            self.make_alignment()
        )
        smoothed = MedianJoiningNetwork(epsilon=2.0, simplify=True).construct_network(
            self.make_alignment()
        )
        assert len(smoothed) <= len(plain)
        assert smoothed.is_connected()

    def test_simplify_preserves_observed_haplotypes(self):
        """Smoothing never removes sampled haplotypes."""
        network = MedianJoiningNetwork(epsilon=3.0, simplify=True).construct_network(
            self.make_alignment()
        )
        for hap_id in ('H1', 'H2', 'H3', 'H4'):
            assert network.has_node(hap_id)

    def test_no_low_degree_medians_remain(self):
        """PopART invariant: every median has degree >= 2."""
        network = MedianJoiningNetwork(epsilon=2.0).construct_network(
            self.make_alignment()
        )
        for node in median_nodes(network):
            assert network.get_degree(node) >= 2


class TestMJNParameters:
    """Parameter reporting."""

    def test_get_parameters(self):
        """All MJN parameters are reported."""
        mjn = MedianJoiningNetwork(epsilon=1.5, max_median_vectors=7, simplify=True)
        params = mjn.get_parameters()
        assert params['epsilon'] == 1.5
        assert params['max_median_vectors'] == 7
        assert params['simplify'] is True

    def test_string_representation(self):
        """String representation includes the class name."""
        assert 'MedianJoiningNetwork' in str(MedianJoiningNetwork())
