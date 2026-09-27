"""Unit tests for N and ? character handling in distance calculations."""

import pytest

from pypopart.core.alignment import Alignment
from pypopart.core.distance import (
    hamming_distance,
    kimura_2p_distance,
    p_distance,
)
from pypopart.core.sequence import Sequence


class TestNCharacterHandling:
    """Test N and ? character handling in distance calculations."""

    def test_n_vs_base_not_counted(self):
        """Test that N vs any base does not count as difference."""
        seq1 = Sequence('s1', 'ATNCG')
        seq2 = Sequence('s2', 'ATCCG')
        assert hamming_distance(seq1, seq2) == 0

    def test_question_vs_base_not_counted(self):
        """Test that ? vs any base does not count as difference."""
        seq1 = Sequence('s1', 'AT?CG')
        seq2 = Sequence('s2', 'ATCCG')
        assert hamming_distance(seq1, seq2) == 0

    def test_multiple_n_characters(self):
        """Test multiple N characters are handled correctly."""
        seq1 = Sequence('s1', 'ANNCG')
        seq2 = Sequence('s2', 'ATCCG')
        assert hamming_distance(seq1, seq2) == 0

    def test_n_vs_n(self):
        """Test N vs N does not count as difference."""
        seq1 = Sequence('s1', 'ATNCG')
        seq2 = Sequence('s2', 'ATNCG')
        assert hamming_distance(seq1, seq2) == 0

    def test_question_vs_question(self):
        """Test ? vs ? does not count as difference."""
        seq1 = Sequence('s1', 'AT?CG')
        seq2 = Sequence('s2', 'AT?CG')
        assert hamming_distance(seq1, seq2) == 0

    def test_n_vs_question(self):
        """Test N vs ? does not count as difference."""
        seq1 = Sequence('s1', 'ATNCG')
        seq2 = Sequence('s2', 'AT?CG')
        assert hamming_distance(seq1, seq2) == 0

    def test_real_difference_with_n_present(self):
        """Test real differences are still counted with N present."""
        seq1 = Sequence('s1', 'ATNCG')
        seq2 = Sequence('s2', 'GTNCG')
        # Only position 0 differs (A vs G), N at position 2 is ignored
        assert hamming_distance(seq1, seq2) == 1

    def test_multiple_differences_with_n(self):
        """Test multiple differences with N in between."""
        seq1 = Sequence('s1', 'ATNCGT')
        seq2 = Sequence('s2', 'GTNCAT')
        # Differences: position 0 (A->G) and position 4 (G->A)
        # Position 2 has N so it's skipped
        assert hamming_distance(seq1, seq2) == 2

    def test_n_with_gaps_both_ignored(self):
        """Test N and gaps can both be ignored."""
        seq1 = Sequence('s1', 'AT-NCG')
        seq2 = Sequence('s2', 'AT-CCG')
        assert hamming_distance(seq1, seq2, ignore_gaps=True) == 0

    def test_n_counted_gaps_counted(self):
        """Test with both N and gaps when gaps are not ignored."""
        seq1 = Sequence('s1', 'AT-NCG')
        seq2 = Sequence('s2', 'ATCCCG')
        # N is still ignored, but gap vs C counts as difference
        assert hamming_distance(seq1, seq2, ignore_gaps=False) == 1

    def test_p_distance_with_n(self):
        """Test p-distance with N characters."""
        seq1 = Sequence('s1', 'ATNCG')
        seq2 = Sequence('s2', 'GTNCG')
        # 1 difference out of 4 comparable positions (N position excluded)
        assert p_distance(seq1, seq2) == 0.25

    def test_p_distance_all_n(self):
        """Test p-distance when all positions have N."""
        seq1 = Sequence('s1', 'NNNNN')
        seq2 = Sequence('s2', 'ATCGT')
        with pytest.raises(ValueError, match='No valid sites to compare'):
            p_distance(seq1, seq2)

    def test_kimura_2p_with_n(self):
        """Test Kimura 2P distance with N characters."""
        seq1 = Sequence('s1', 'ATNCG')
        seq2 = Sequence('s2', 'GTNCG')
        # Should handle N correctly (already did before)
        dist = kimura_2p_distance(seq1, seq2)
        assert dist >= 0  # Just check it computes without error

    def test_alignment_with_n_characters(self):
        """Test that alignments with N characters are valid."""
        seqs = [
            Sequence('seq1', 'ATNCG'),
            Sequence('seq2', 'GTNCG'),
            Sequence('seq3', 'ATCCG'),
        ]
        aln = Alignment(seqs)
        assert aln.is_valid()
        assert len(aln) == 3
        assert aln.length == 5

    def test_all_positions_with_n_or_gaps(self):
        """Test sequences where all positions have N or gaps."""
        seq1 = Sequence('s1', 'N-N-N')
        seq2 = Sequence('s2', '-N-N-')
        # All positions have either N or gap, so no comparable sites
        assert hamming_distance(seq1, seq2, ignore_gaps=True) == 0

    def test_lowercase_n_converted_to_uppercase(self):
        """Test that lowercase n is converted to uppercase N by Sequence class."""
        seq1 = Sequence('s1', 'ATnCG')
        seq2 = Sequence('s2', 'ATCCG')
        # Sequence class converts to uppercase, so 'n' becomes 'N'
        # N vs C should not count as difference
        assert hamming_distance(seq1, seq2) == 0

    def test_mixed_case_normalized(self):
        """Test mixed case sequences are normalized to uppercase."""
        seq1 = Sequence('s1', 'ATNCG')
        seq2 = Sequence('s2', 'atccg')
        # Sequence class normalizes to uppercase
        # After normalization: ATNCG vs ATCCG
        # N is skipped, so distance is 0
        assert hamming_distance(seq1, seq2) == 0

    def test_n_at_different_positions(self):
        """Test N at different positions in sequences."""
        seq1 = Sequence('s1', 'NATCG')
        seq2 = Sequence('s2', 'ATCGN')
        # Both positions with N are skipped
        # Positions 1-3 in seq1 vs positions 0-2 in seq2
        # Position 0 of seq1 (N) vs position 0 of seq2 (A): skipped
        # Position 1 of seq1 (A) vs position 1 of seq2 (T): different
        # Position 2 of seq1 (T) vs position 2 of seq2 (C): different
        # Position 3 of seq1 (C) vs position 3 of seq2 (G): different
        # Position 4 of seq1 (G) vs position 4 of seq2 (N): skipped
        assert hamming_distance(seq1, seq2) == 3


class TestGapsAsValidCharacters:
    """Test that gaps are treated as valid characters when not ignored."""

    def test_alignment_with_gaps_accepted(self):
        """Test that alignments with gaps are accepted."""
        seqs = [
            Sequence('seq1', 'AT-CG'),
            Sequence('seq2', 'ATCCG'),
            Sequence('seq3', 'AT-CT'),
        ]
        aln = Alignment(seqs)
        assert aln.is_valid()
        assert len(aln) == 3
        assert aln.length == 5

    def test_gap_vs_base_counts_as_difference(self):
        """Test that gap vs base counts as difference when gaps not ignored."""
        seq1 = Sequence('s1', 'AT-CG')
        seq2 = Sequence('s2', 'ATCCG')
        assert hamming_distance(seq1, seq2, ignore_gaps=False) == 1

    def test_gap_vs_gap_no_difference(self):
        """Test that gap vs gap is not a difference."""
        seq1 = Sequence('s1', 'AT-CG')
        seq2 = Sequence('s2', 'AT-CG')
        assert hamming_distance(seq1, seq2, ignore_gaps=False) == 0

    def test_multiple_gaps(self):
        """Test sequences with multiple gaps."""
        seq1 = Sequence('s1', 'A--CG')
        seq2 = Sequence('s2', 'ATCCG')
        assert hamming_distance(seq1, seq2, ignore_gaps=False) == 2

    def test_all_gaps(self):
        """Test sequences that are all gaps."""
        seq1 = Sequence('s1', '-----')
        seq2 = Sequence('s2', '-----')
        assert hamming_distance(seq1, seq2, ignore_gaps=False) == 0

    def test_gaps_ignored_by_default(self):
        """Test that gaps are ignored by default."""
        seq1 = Sequence('s1', 'AT-CG')
        seq2 = Sequence('s2', 'AT-CT')
        # With gaps ignored, only position 3 differs (G vs T)
        assert hamming_distance(seq1, seq2, ignore_gaps=True) == 1

    def test_p_distance_with_gaps_not_ignored(self):
        """Test p-distance when gaps are not ignored."""
        seq1 = Sequence('s1', 'AT-CG')
        seq2 = Sequence('s2', 'ATCCG')
        # When gaps not ignored, 1 difference out of 5 positions
        assert p_distance(seq1, seq2, ignore_gaps=False) == 0.2

    def test_length_validation_includes_gaps(self):
        """Test that length validation includes gap characters."""
        seq1 = Sequence('s1', 'AT-CG')
        seq2 = Sequence('s2', 'ATCG')  # Different length
        with pytest.raises(ValueError, match='same length'):
            hamming_distance(seq1, seq2)


class TestNumbaOptimizedVersions:
    """Test that Numba-optimized versions handle N correctly."""

    def test_numba_n_handling(self):
        """Test that Numba version handles N correctly."""
        seq1 = Sequence('s1', 'ATNCG')
        seq2 = Sequence('s2', 'ATCCG')
        # Force use of Numba version
        assert hamming_distance(seq1, seq2, use_numba=True) == 0

    def test_numba_with_real_differences_and_n(self):
        """Test Numba version with real differences and N."""
        seq1 = Sequence('s1', 'ATNCG')
        seq2 = Sequence('s2', 'GTNCG')
        assert hamming_distance(seq1, seq2, use_numba=True) == 1

    def test_numba_multiple_n(self):
        """Test Numba version with multiple N characters."""
        seq1 = Sequence('s1', 'ANNCG')
        seq2 = Sequence('s2', 'ATCCG')
        assert hamming_distance(seq1, seq2, use_numba=True) == 0
