"""Tests for site-pattern condensation (port of HapNet::condenseSitePats)."""

import numpy as np

from pypopart.core.distance import pairwise_distance_matrix
from pypopart.core.sequence import Sequence
from pypopart.core.site_patterns import condense_site_patterns, is_ambiguous


class TestCondenseSitePatterns:
    """Hand-worked condensation cases mirroring the C++ semantics."""

    def test_identical_columns_merge(self):
        """Byte-identical columns collapse with summed weight."""
        # Columns 0/1 identical (A,A,C); columns 2/3 identical (T,G,T)
        condensed, weights = condense_site_patterns(['AATT', 'AAGG', 'CCTT'])
        assert condensed == ['AT', 'AG', 'CT']
        assert weights == [2, 2]

    def test_bijective_relabelling_merges(self):
        """Columns equivalent up to a character bijection collapse."""
        # col0 = (A,A,C), col1 = (G,G,T): A->G, C->T is a bijection
        condensed, weights = condense_site_patterns(['AG', 'AG', 'CT'])
        assert condensed == ['A', 'A', 'C']
        assert weights == [2]

    def test_non_bijective_stays_separate(self):
        """A many-to-one mapping is not an equivalence."""
        # col0 = (A,A,C), col1 = (G,T,T): A maps to both G and T
        condensed, weights = condense_site_patterns(['AG', 'AT', 'CT'])
        assert condensed == ['AG', 'AT', 'CT']
        assert weights == [1, 1]

    def test_all_ambiguous_column_dropped(self):
        """Columns ambiguous in every sequence are removed."""
        # Columns 1 (all N) and 2 (all -) drop; columns 0 and 3 are
        # identical patterns and merge with weight 2.
        condensed, weights = condense_site_patterns(['AN-A', 'CN-C', 'GN-G'])
        assert condensed == ['A', 'C', 'G']
        assert weights == [2]

    def test_ambiguous_chars_must_match_exactly(self):
        """A column with N pairs only with a column carrying the same N."""
        # col0 = (A,N,C) vs col1 = (A,A,C): differs where col0 is ambiguous
        condensed, weights = condense_site_patterns(['AA', 'NA', 'CC'])
        assert condensed == ['AA', 'NA', 'CC']
        assert weights == [1, 1]

        # Identical ambiguous placement does merge
        condensed, weights = condense_site_patterns(['AA', 'NN', 'CC'])
        assert condensed == ['A', 'N', 'C']
        assert weights == [2]

    def test_constant_columns_merge_across_letters(self):
        """Constant columns of different letters are pattern-equivalent."""
        condensed, weights = condense_site_patterns(['ATG', 'ATG', 'ATG'])
        assert condensed == ['A', 'A', 'A']
        assert weights == [3]

    def test_empty_input(self):
        """Empty input yields empty output."""
        assert condense_site_patterns([]) == ([], [])

    def test_is_ambiguous(self):
        """The ambiguity set covers gaps, N, and IUPAC codes."""
        for char in '-NYRMSVWKDHB?X':
            assert is_ambiguous(char)
        for char in 'ACGTU':
            assert not is_ambiguous(char)


class TestWeightedDistanceInvariance:
    """Weighted distances over condensed columns equal plain distances."""

    def test_pairwise_distances_unchanged_by_condensation(self):
        """For sampled sequences, condensation must not change distances."""
        raw = ['AATTGC', 'AAGGGC', 'CCTTAT', 'CCGGAT']
        sequences = [Sequence(f's{i}', s) for i, s in enumerate(raw)]
        plain = pairwise_distance_matrix(sequences)

        condensed, weights = condense_site_patterns(raw)
        condensed_seqs = [Sequence(f's{i}', s) for i, s in enumerate(condensed)]
        weighted = pairwise_distance_matrix(
            condensed_seqs, site_weights=np.asarray(weights)
        )

        np.testing.assert_array_equal(plain.matrix, weighted.matrix)


class TestMaskSupport:
    """Character masking excludes columns from the calculation."""

    def test_mask_excludes_columns(self):
        """Masked-out columns contribute nothing to distances."""
        sequences = [Sequence('a', 'AAAA'), Sequence('b', 'TTAA')]
        full = pairwise_distance_matrix(sequences)
        assert full.matrix[0, 1] == 2

        masked = pairwise_distance_matrix(
            sequences, mask=np.array([False, True, True, True])
        )
        assert masked.matrix[0, 1] == 1
