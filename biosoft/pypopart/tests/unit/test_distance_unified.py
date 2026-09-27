"""Tests for the unified pairwise distance path."""

import random

import numpy as np
import pytest

from pypopart.core import distance as distance_mod
from pypopart.core.alignment import Alignment
from pypopart.core.distance import (
    DistanceMatrix,
    hamming_distance,
    normalize_distance_method,
    pairwise_distance_matrix,
    sequence_distance,
)
from pypopart.core.haplotype import identify_haplotypes_from_alignment
from pypopart.core.sequence import Sequence


def random_alignment(n_seqs: int, length: int, seed: int) -> Alignment:
    """Build a random DNA alignment with occasional gaps and Ns."""
    rng = random.Random(seed)
    alphabet = 'ACGT' * 5 + '-N'
    return Alignment(
        [
            Sequence(f's{i}', ''.join(rng.choice(alphabet) for _ in range(length)))
            for i in range(n_seqs)
        ]
    )


class TestHammingPathEquivalence:
    """The numba kernel, numpy fallback, and scalar loop must agree."""

    @pytest.mark.parametrize('seed', [1, 2, 3])
    @pytest.mark.parametrize('ignore_gaps', [True, False])
    def test_matrix_paths_match_scalar(self, monkeypatch, seed, ignore_gaps):
        """Fast matrix paths equal per-pair scalar hamming distances."""
        alignment = random_alignment(8, 40, seed)
        seqs = list(alignment)
        n = len(seqs)

        expected = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                d = hamming_distance(
                    seqs[i], seqs[j], ignore_gaps=ignore_gaps, use_numba=False
                )
                expected[i, j] = expected[j, i] = d

        # Default path (numba when available)
        result = pairwise_distance_matrix(alignment, ignore_gaps=ignore_gaps)
        np.testing.assert_array_equal(result.matrix, expected)

        # Forced numpy fallback
        monkeypatch.setattr(distance_mod, '_NUMBA_AVAILABLE', False)
        fallback = pairwise_distance_matrix(alignment, ignore_gaps=ignore_gaps)
        np.testing.assert_array_equal(fallback.matrix, expected)

    def test_accepts_haplotype_list(self):
        """Haplotype lists work and are labelled by haplotype id."""
        alignment = random_alignment(6, 20, seed=7)
        haplotypes = identify_haplotypes_from_alignment(alignment)

        result = pairwise_distance_matrix(haplotypes)
        assert isinstance(result, DistanceMatrix)
        assert result.labels == [h.id for h in haplotypes]
        assert result.matrix.shape == (len(haplotypes), len(haplotypes))
        np.testing.assert_array_equal(result.matrix, result.matrix.T)

    def test_ambiguity_and_gap_semantics(self):
        """N/? never count as differences; gaps only when not ignored."""
        seqs = [Sequence('a', 'ANC-'), Sequence('b', 'GNCT')]
        with_gaps_ignored = pairwise_distance_matrix(seqs, ignore_gaps=True)
        # Position 0 differs (A/G); N skipped; C matches; gap skipped
        assert with_gaps_ignored.matrix[0, 1] == 1

        counting_gaps = pairwise_distance_matrix(seqs, ignore_gaps=False)
        # Gap position now counts as a difference
        assert counting_gaps.matrix[0, 1] == 2


class TestCorrectedDistances:
    """Named substitution models run through the same entry point."""

    @pytest.mark.parametrize('method', ['jc', 'k2p', 'tn', 'p'])
    def test_corrected_methods(self, method):
        """Corrected models produce symmetric, zero-diagonal matrices."""
        alignment = random_alignment(5, 60, seed=11)
        result = pairwise_distance_matrix(alignment, method=method)
        np.testing.assert_array_equal(result.matrix, result.matrix.T)
        assert np.all(np.diag(result.matrix) == 0)

    def test_alias_normalization(self):
        """tamura_nei and tn resolve to the same computation."""
        alignment = random_alignment(4, 50, seed=13)
        via_alias = pairwise_distance_matrix(alignment, method='tamura_nei')
        via_canonical = pairwise_distance_matrix(alignment, method='tn')
        np.testing.assert_array_equal(via_alias.matrix, via_canonical.matrix)

    def test_unknown_method_raises(self):
        """Unknown method names fail loudly."""
        with pytest.raises(ValueError, match='Unknown distance method'):
            normalize_distance_method('bogus')


class TestAlgorithmsHonourDistanceMethod:
    """distance_method now reaches every algorithm's haplotype distances."""

    def test_k2p_differs_from_hamming(self):
        """A k2p-configured algorithm computes k2p haplotype distances."""
        from pypopart.algorithms import MinimumSpanningTree

        alignment = random_alignment(6, 80, seed=17)
        haplotypes = identify_haplotypes_from_alignment(alignment)

        hamming = MinimumSpanningTree(distance_method='hamming')
        k2p = MinimumSpanningTree(distance_method='k2p')

        m1 = hamming.calculate_haplotype_distances(haplotypes).matrix
        m2 = k2p.calculate_haplotype_distances(haplotypes).matrix
        assert not np.array_equal(m1, m2)

    def test_scalar_sequence_distance(self):
        """sequence_distance dispatches by method name."""
        s1, s2 = Sequence('a', 'AAAA'), Sequence('b', 'AAAT')
        assert sequence_distance(s1, s2, method='hamming') == 1.0
        assert sequence_distance(s1, s2, method='p') == pytest.approx(0.25)
