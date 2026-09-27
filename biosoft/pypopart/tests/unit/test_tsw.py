"""Unit tests for Tight Span Walker (TSW) algorithm."""

import numpy as np

from pypopart.algorithms.tsw import TightSpanWalker
from pypopart.core.alignment import Alignment
from pypopart.core.distance import DistanceMatrix
from pypopart.core.sequence import Sequence


class TestTightSpanWalker:
    """Test cases for TSW algorithm."""

    def test_tsw_initialization(self):
        """Test TSW algorithm initialization."""
        tsw = TightSpanWalker()
        assert tsw.distance_method == 'hamming'
        assert tsw.tolerance == 1e-6

        tsw_custom = TightSpanWalker(tolerance=1e-5)
        assert tsw_custom.tolerance == 1e-5

    def test_tsw_empty_alignment(self):
        """Test TSW with empty alignment."""
        tsw = TightSpanWalker()
        alignment = Alignment()
        network = tsw.construct_network(alignment)

        assert len(network) == 0

    def test_tsw_single_sequence(self):
        """Test TSW with single sequence."""
        tsw = TightSpanWalker()
        alignment = Alignment([Sequence('seq1', 'ATCG')])
        network = tsw.construct_network(alignment)

        assert len(network) == 1
        assert len(network.edges) == 0

    def test_tsw_two_sequences(self):
        """Test TSW with two sequences."""
        tsw = TightSpanWalker()
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCC'),  # 1 difference
            ]
        )
        network = tsw.construct_network(alignment)

        # Should have at least 2 nodes (may add median vertices)
        assert len(network) >= 2
        # Should have edges connecting them
        assert len(network.edges) >= 1
        assert network.is_connected()

    def test_tsw_three_sequences_linear(self):
        """Test TSW with three sequences in linear arrangement."""
        tsw = TightSpanWalker()
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),  # 1 diff from seq1
                Sequence('seq3', 'AATT'),  # 1 diff from seq2, 2 from seq1
            ]
        )
        network = tsw.construct_network(alignment)

        assert len(network) >= 3
        assert network.is_connected()

    def test_tsw_three_sequences_star(self):
        """Test TSW with three sequences in star topology."""
        tsw = TightSpanWalker()
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),  # Central
                Sequence('seq2', 'AAAT'),  # 1 diff
                Sequence('seq3', 'AAAC'),  # 1 diff
            ]
        )
        network = tsw.construct_network(alignment)

        assert len(network) >= 3
        assert network.is_connected()

    def test_tsw_compute_dt(self):
        """DT includes k == i/j, so dT(i, j) >= d(i, j) always."""
        tsw = TightSpanWalker()
        tsw._n_samples = 3
        import numpy as np

        tsw._d = np.array([[0.0, 1.0, 2.0], [1.0, 0.0, 2.0], [2.0, 2.0, 0.0]])
        tsw._compute_dt()

        # Symmetric and bounded below by d
        for i in range(3):
            for j in range(3):
                assert tsw._dt_get(i, j) == tsw._dt_get(j, i)
                assert tsw._dt_get(i, j) >= tsw._d[i, j]

        # dT(0,1) = max(|0-1|, |1-0|, |2-2|) = 1 (k==i term dominates)
        assert tsw._dt_get(0, 1) == 1.0

    def test_tsw_internal_vertices_sequenceless(self):
        """Internal tight-span vertices carry no sequence (C++ parity)."""
        tsw = TightSpanWalker()
        alignment = Alignment(
            [
                Sequence('s1', 'AAT'),
                Sequence('s2', 'ATA'),
                Sequence('s3', 'TAA'),
            ]
        )
        network = tsw.construct_network(alignment)
        medians = [
            node
            for node, attrs in network.graph.nodes(data=True)
            if attrs.get('median_vector')
        ]
        for node in medians:
            assert network.get_haplotype(node).data == ''
            assert network.get_haplotype(node).frequency == 0

    def test_tsw_with_identical_sequences(self):
        """Test TSW with identical sequences (should be single haplotype)."""
        tsw = TightSpanWalker()
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCG'),
                Sequence('seq3', 'ATCG'),
            ]
        )
        network = tsw.construct_network(alignment)

        # Should collapse to single haplotype with frequency 3
        assert len(network) == 1
        assert len(network.edges) == 0
        hap = network.haplotypes[0]
        assert hap.frequency == 3

    def test_tsw_with_provided_distance_matrix(self):
        """Test TSW with pre-computed distance matrix."""
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCC'),
                Sequence('seq3', 'GTCG'),
            ]
        )

        # Create custom distance matrix with haplotype IDs
        # Note: TSW will create haplotypes H1, H2, H3
        labels = ['H1', 'H2', 'H3']
        matrix = np.array([[0.0, 1.0, 1.0], [1.0, 0.0, 2.0], [1.0, 2.0, 0.0]])
        dist_matrix = DistanceMatrix(labels, matrix)

        tsw = TightSpanWalker()
        network = tsw.construct_network(alignment, distance_matrix=dist_matrix)

        assert len(network) >= 3
        assert network.is_connected()

    def test_tsw_four_sequences_complex(self):
        """Test TSW with four sequences requiring complex network."""
        tsw = TightSpanWalker()
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),
                Sequence('seq3', 'ATAA'),
                Sequence('seq4', 'TAAA'),
            ]
        )
        network = tsw.construct_network(alignment)

        # Should have 4 original haplotypes
        original_haps = [
            h for h in network.haplotypes if not network.is_median_vector(h.id)
        ]
        assert len(original_haps) == 4

        # Network should be connected
        assert network.is_connected()

    def test_tsw_inferred_vs_observed_haplotypes(self):
        """Test that TSW correctly marks inferred vs observed haplotypes."""
        tsw = TightSpanWalker()
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AACC'),  # 2 differences
            ]
        )
        network = tsw.construct_network(alignment)

        # Count observed vs inferred (using median_vector flag)
        observed = [h for h in network.haplotypes if not network.is_median_vector(h.id)]
        inferred = [h.id for h in network.haplotypes if network.is_median_vector(h.id)]

        # Should have 2 observed
        assert len(observed) == 2

        # All inferred should have frequency 0
        for hap_id in inferred:
            hap = next(h for h in network.haplotypes if h.id == hap_id)
            assert hap.frequency == 0
            assert len(hap.sample_ids) == 0

    def test_tsw_parameters(self):
        """Test getting TSW parameters."""
        tsw = TightSpanWalker(distance_method='k2p', tolerance=1e-5)
        params = tsw.get_parameters()

        assert params['distance_method'] == 'k2p'
        assert params['tolerance'] == 1e-5

    def test_tsw_string_representation(self):
        """Test string representation of TSW."""
        tsw = TightSpanWalker(distance_method='hamming', tolerance=1e-6)
        str_repr = str(tsw)

        assert 'TightSpanWalker' in str_repr
        assert 'hamming' in str_repr
        assert '1e-06' in str_repr or '1e-6' in str_repr

    def test_tsw_handles_gaps(self):
        """Test TSW handles sequences with gaps."""
        tsw = TightSpanWalker()
        alignment = Alignment(
            [
                Sequence('seq1', 'AT-G'),
                Sequence('seq2', 'ATCG'),
            ]
        )
        network = tsw.construct_network(alignment)

        # Should handle gaps in distance calculation
        # May add median vertices
        assert len(network) >= 2
        assert network.is_connected()

    def test_tsw_longer_sequences(self):
        """Test TSW with longer sequences."""
        tsw = TightSpanWalker()
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCGATCGATCG'),
                Sequence('seq2', 'ATCGATCGATCC'),  # 1 diff at end
                Sequence('seq3', 'ATCGATCCATCG'),  # 2 diffs in middle
            ]
        )
        network = tsw.construct_network(alignment)

        assert len(network) >= 3
        assert network.is_connected()

    def test_tsw_different_tolerance_values(self):
        """Test TSW with different tolerance values."""
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),
                Sequence('seq2', 'AAAT'),
                Sequence('seq3', 'AATT'),
            ]
        )

        # Test with strict tolerance
        tsw1 = TightSpanWalker(tolerance=1e-10)
        network1 = tsw1.construct_network(alignment)

        # Test with relaxed tolerance
        tsw2 = TightSpanWalker(tolerance=1.0)
        network2 = tsw2.construct_network(alignment)

        # Both should create valid networks
        assert network1.is_connected()
        assert network2.is_connected()

        # Networks may have different structures
        # but both should have at least the original haplotypes
        assert len(network1) >= 3
        assert len(network2) >= 3
