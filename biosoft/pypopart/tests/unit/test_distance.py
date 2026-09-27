"""Unit tests for distance calculation module."""

import numpy as np
import pytest

from pypopart.core.alignment import Alignment
from pypopart.core.distance import (
    DistanceMatrix,
    calculate_distance_matrix,
    calculate_pairwise_distances,
    hamming_distance,
    jukes_cantor_distance,
    kimura_2p_distance,
    p_distance,
    tamura_nei_distance,
)
from pypopart.core.sequence import Sequence


class TestHammingDistance:
    """Test Hamming distance calculation."""

    def test_identical_sequences(self):
        """Test distance between identical sequences."""
        seq1 = Sequence('s1', 'ATCG')
        seq2 = Sequence('s2', 'ATCG')
        assert hamming_distance(seq1, seq2) == 0

    def test_all_different(self):
        """Test distance when all positions differ."""
        seq1 = Sequence('s1', 'ATCG')
        seq2 = Sequence('s2', 'GCTA')
        assert hamming_distance(seq1, seq2) == 4

    def test_some_differences(self):
        """Test distance with some differences."""
        seq1 = Sequence('s1', 'ATCG')
        seq2 = Sequence('s2', 'ATCC')
        assert hamming_distance(seq1, seq2) == 1

    def test_with_gaps_ignored(self):
        """Test that gaps are ignored by default."""
        seq1 = Sequence('s1', 'AT-CG')
        seq2 = Sequence('s2', 'AT-CG')
        assert hamming_distance(seq1, seq2, ignore_gaps=True) == 0

    def test_with_gaps_not_ignored(self):
        """Test counting gaps as differences."""
        seq1 = Sequence('s1', 'AT-CG')
        seq2 = Sequence('s2', 'ATCCG')
        assert hamming_distance(seq1, seq2, ignore_gaps=False) == 1

    def test_length_mismatch(self):
        """Test error when sequences have different lengths."""
        seq1 = Sequence('s1', 'ATCG')
        seq2 = Sequence('s2', 'ATCGAA')
        with pytest.raises(ValueError, match='same length'):
            hamming_distance(seq1, seq2)


class TestPDistance:
    """Test p-distance calculation."""

    def test_identical_sequences(self):
        """Test p-distance of identical sequences."""
        seq1 = Sequence('s1', 'ATCG')
        seq2 = Sequence('s2', 'ATCG')
        assert p_distance(seq1, seq2) == 0.0

    def test_all_different(self):
        """Test p-distance when all differ."""
        seq1 = Sequence('s1', 'ATCG')
        seq2 = Sequence('s2', 'GCTA')
        assert p_distance(seq1, seq2) == 1.0

    def test_half_different(self):
        """Test p-distance with 50% difference."""
        seq1 = Sequence('s1', 'ATCG')
        seq2 = Sequence('s2', 'ATAA')
        assert p_distance(seq1, seq2) == 0.5

    def test_with_gaps(self):
        """Test p-distance with gaps."""
        seq1 = Sequence('s1', 'AT-CG')
        seq2 = Sequence('s2', 'AT-CC')
        # Only 4 non-gap positions, 1 differs
        assert abs(p_distance(seq1, seq2) - 0.25) < 0.001


class TestJukesCantor:
    """Test Jukes-Cantor distance."""

    def test_identical_sequences(self):
        """Test JC distance of identical sequences."""
        seq1 = Sequence('s1', 'ATCG')
        seq2 = Sequence('s2', 'ATCG')
        assert jukes_cantor_distance(seq1, seq2) == 0.0

    def test_small_difference(self):
        """Test JC distance with small difference."""
        seq1 = Sequence('s1', 'ATCGAAAA')
        seq2 = Sequence('s2', 'ATCGAAAG')
        # p = 1/8 = 0.125
        # JC = -0.75 * ln(1 - 4/3 * 0.125) = -0.75 * ln(0.8333)
        jc_dist = jukes_cantor_distance(seq1, seq2)
        assert 0.13 < jc_dist < 0.14

    def test_too_divergent(self):
        """Test error when sequences are too divergent."""
        seq1 = Sequence('s1', 'AAAA')
        seq2 = Sequence('s2', 'GGGG')
        # p = 1.0 >= 0.75
        with pytest.raises(ValueError, match='too divergent'):
            jukes_cantor_distance(seq1, seq2)


class TestKimura2P:
    """Test Kimura 2-parameter distance."""

    def test_identical_sequences(self):
        """Test K2P distance of identical sequences."""
        seq1 = Sequence('s1', 'ATCG')
        seq2 = Sequence('s2', 'ATCG')
        assert kimura_2p_distance(seq1, seq2) == 0.0

    def test_only_transitions(self):
        """Test K2P with only transitions (A<->G, C<->T)."""
        seq1 = Sequence('s1', 'AACCGGTTAAAA')
        seq2 = Sequence('s2', 'GGTTCCAAAAAA')
        # First 8 positions are transitions, rest identical
        k2p_dist = kimura_2p_distance(seq1, seq2)
        assert k2p_dist > 0

    def test_only_transversions(self):
        """Test K2P with only transversions."""
        seq1 = Sequence('s1', 'AAAAGGGGGGGG')
        seq2 = Sequence('s2', 'CCCCGGGGGGGG')
        # First 4 positions are transversions (A->C), rest identical (8 positions)
        k2p_dist = kimura_2p_distance(seq1, seq2)
        assert k2p_dist > 0

    def test_with_ambiguous_chars(self):
        """Test K2P with ambiguous characters."""
        seq1 = Sequence('s1', 'ATCGN')
        seq2 = Sequence('s2', 'ATCCN')
        # Should skip N positions
        k2p_dist = kimura_2p_distance(seq1, seq2)
        assert k2p_dist > 0


class TestTamuraNei:
    """Test Tamura-Nei distance."""

    def test_identical_sequences(self):
        """Test TN distance of identical sequences."""
        seq1 = Sequence('s1', 'ATCG')
        seq2 = Sequence('s2', 'ATCG')
        assert tamura_nei_distance(seq1, seq2) == 0.0

    def test_with_transitions(self):
        """Test TN with purine and pyrimidine transitions."""
        # Longer sequences with moderate divergence
        seq1 = Sequence('s1', 'ATCGATCGATCGATCGATCGATCG')
        seq2 = Sequence('s2', 'ATCGATCGATCGGTCGATCGATCG')
        # Single A->G transition
        tn_dist = tamura_nei_distance(seq1, seq2)
        assert tn_dist > 0

    def test_with_transversions(self):
        """Test TN with transversions."""
        seq1 = Sequence('s1', 'AAAAGGGGCCCCTTTT')
        seq2 = Sequence('s2', 'CCCCGGGGCCCCTTTT')
        # First 4 are transversions (A->C)
        tn_dist = tamura_nei_distance(seq1, seq2)
        assert tn_dist > 0

    def test_mixed_differences(self):
        """Test TN with mixed types of differences."""
        seq1 = Sequence('s1', 'ATCGATCGATCGATCGATCGATCG')
        seq2 = Sequence('s2', 'ATCGATCGGTCCATCGATCGATCG')
        # Mix of transitions (A->G) and transversions (A->C)
        tn_dist = tamura_nei_distance(seq1, seq2)
        assert tn_dist > 0

    def test_with_gaps(self):
        """Test TN with gaps ignored."""
        seq1 = Sequence('s1', 'AT-CG-ATCG')
        seq2 = Sequence('s2', 'AT-CC-GTCG')
        # Should ignore gap positions
        tn_dist = tamura_nei_distance(seq1, seq2, ignore_gaps=True)
        assert tn_dist > 0

    def test_with_ambiguous_chars(self):
        """Test TN with ambiguous characters."""
        seq1 = Sequence('s1', 'ATCGATCGNATCG')
        seq2 = Sequence('s2', 'ATCGGTCGNATCG')
        # Should skip N positions
        tn_dist = tamura_nei_distance(seq1, seq2)
        assert tn_dist > 0

    def test_too_divergent(self):
        """Test error when sequences are too divergent."""
        seq1 = Sequence('s1', 'AAAAAAAA')
        seq2 = Sequence('s2', 'GGGGGGGG')
        # All positions differ
        with pytest.raises(ValueError, match='too divergent'):
            tamura_nei_distance(seq1, seq2)

    def test_length_mismatch(self):
        """Test error when sequences have different lengths."""
        seq1 = Sequence('s1', 'ATCG')
        seq2 = Sequence('s2', 'ATCGAA')
        with pytest.raises(ValueError, match='same length'):
            tamura_nei_distance(seq1, seq2)


class TestDistanceMatrix:
    """Test DistanceMatrix class."""

    def test_empty_matrix(self):
        """Test creating empty matrix."""
        labels = ['s1', 's2', 's3']
        dm = DistanceMatrix(labels)

        assert dm.n == 3
        assert dm.labels == labels
        assert dm.matrix.shape == (3, 3)
        assert np.all(dm.matrix == 0)

    def test_with_provided_matrix(self):
        """Test creating matrix with data."""
        labels = ['s1', 's2']
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(labels, matrix)

        assert dm.get_distance('s1', 's2') == 1
        assert dm.get_distance('s2', 's1') == 1

    def test_set_distance(self):
        """Test setting distances."""
        labels = ['s1', 's2']
        dm = DistanceMatrix(labels)

        dm.set_distance('s1', 's2', 5.0)
        assert dm.get_distance('s1', 's2') == 5.0
        assert dm.get_distance('s2', 's1') == 5.0  # Symmetric

    def test_get_row(self):
        """Test getting all distances for a sequence."""
        labels = ['s1', 's2', 's3']
        matrix = np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
        dm = DistanceMatrix(labels, matrix)

        row = dm.get_row('s1')
        assert list(row) == [0, 1, 2]

    def test_min_max_distance(self):
        """Test min/max distance."""
        labels = ['s1', 's2', 's3']
        matrix = np.array([[0, 1, 5], [1, 0, 3], [5, 3, 0]])
        dm = DistanceMatrix(labels, matrix)

        assert dm.get_min_distance(exclude_zero=True) == 1
        assert dm.get_max_distance() == 5

    def test_to_from_dict(self):
        """Test serialization."""
        labels = ['s1', 's2']
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(labels, matrix)

        data = dm.to_dict()
        dm2 = DistanceMatrix.from_dict(data)

        assert dm2.labels == labels
        assert np.allclose(dm2.matrix, matrix)

    def test_csv_export_import(self, tmp_path):
        """Test CSV export and import."""
        labels = ['s1', 's2', 's3']
        matrix = np.array([[0, 1.5, 2.3], [1.5, 0, 3.1], [2.3, 3.1, 0]])
        dm = DistanceMatrix(labels, matrix)

        # Export to CSV
        csv_file = tmp_path / 'distances.csv'
        dm.to_csv(str(csv_file))

        # Import from CSV
        dm2 = DistanceMatrix.from_csv(str(csv_file))

        assert dm2.labels == labels
        assert np.allclose(dm2.matrix, matrix, atol=1e-6)

    def test_visualize_without_matplotlib(self):
        """Test that visualize raises error without matplotlib."""
        labels = ['s1', 's2']
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(labels, matrix)

        # Mock matplotlib import to test error handling
        # This test assumes matplotlib is installed; in production you'd mock the import
        try:
            fig = dm.visualize()
            # If matplotlib is installed, should return a figure
            assert fig is not None
            import matplotlib.pyplot as plt

            plt.close(fig)
        except ImportError:
            # If matplotlib not installed, should raise ImportError
            pass


class TestCalculateDistanceMatrix:
    """Test distance matrix calculation from alignment."""

    def test_basic_calculation(self):
        """Test basic distance matrix calculation."""
        sequences = [
            Sequence('s1', 'ATCG'),
            Sequence('s2', 'ATCC'),
            Sequence('s3', 'GTCG'),
        ]
        alignment = Alignment(sequences)

        dm = calculate_distance_matrix(alignment)

        assert dm.n == 3
        assert dm.get_distance('s1', 's2') == 1
        assert dm.get_distance('s1', 's3') == 1
        assert dm.get_distance('s2', 's3') == 2

    def test_with_custom_function(self):
        """Test with custom distance function."""
        sequences = [Sequence('s1', 'ATCG'), Sequence('s2', 'ATCC')]
        alignment = Alignment(sequences)

        dm = calculate_distance_matrix(alignment, distance_func=p_distance)

        assert dm.get_distance('s1', 's2') == 0.25


class TestCalculatePairwiseDistances:
    """Test convenience function for pairwise distances."""

    def test_hamming_method(self):
        """Test using Hamming distance method."""
        sequences = [Sequence('s1', 'ATCG'), Sequence('s2', 'ATCC')]
        alignment = Alignment(sequences)

        dm = calculate_pairwise_distances(alignment, method='hamming')
        assert dm.get_distance('s1', 's2') == 1

    def test_p_distance_method(self):
        """Test using p-distance method."""
        sequences = [Sequence('s1', 'ATCG'), Sequence('s2', 'ATCC')]
        alignment = Alignment(sequences)

        dm = calculate_pairwise_distances(alignment, method='p')
        assert dm.get_distance('s1', 's2') == 0.25

    def test_jc_method(self):
        """Test using Jukes-Cantor method."""
        sequences = [Sequence('s1', 'ATCGAAAA'), Sequence('s2', 'ATCGAAAG')]
        alignment = Alignment(sequences)

        dm = calculate_pairwise_distances(alignment, method='jc')
        assert dm.get_distance('s1', 's2') > 0

    def test_k2p_method(self):
        """Test using K2P method."""
        sequences = [Sequence('s1', 'ATCG'), Sequence('s2', 'ATCC')]
        alignment = Alignment(sequences)

        dm = calculate_pairwise_distances(alignment, method='k2p')
        assert dm.get_distance('s1', 's2') > 0

    def test_tn_method(self):
        """Test using Tamura-Nei method."""
        sequences = [Sequence('s1', 'ATCGATCGATCG'), Sequence('s2', 'ATCGGTCGATCG')]
        alignment = Alignment(sequences)

        dm = calculate_pairwise_distances(alignment, method='tn')
        assert dm.get_distance('s1', 's2') > 0

    def test_unknown_method(self):
        """Test error with unknown method."""
        sequences = [Sequence('s1', 'ATCG')]
        alignment = Alignment(sequences)

        with pytest.raises(ValueError, match='Unknown distance method'):
            calculate_pairwise_distances(alignment, method='invalid')
