"""
Performance benchmarks for distance calculation and network algorithms.

Excluded from the default test run (pytest addopts select `-m 'not
benchmark'`). Run explicitly with:

    pytest tests/benchmarks -m benchmark --benchmark-json=bench.json
"""

import random

import pytest

from pypopart.core.alignment import Alignment
from pypopart.core.distance import pairwise_distance_matrix
from pypopart.core.sequence import Sequence

pytestmark = pytest.mark.benchmark


def synthetic_alignment(n_seqs: int, length: int, seed: int = 42) -> Alignment:
    """
    Generate a reproducible synthetic alignment with realistic variation.

    Mutates a random ancestral sequence at ~5% of sites per sequence so
    haplotype structure resembles real intraspecific data.

    Parameters
    ----------
    n_seqs : int
        Number of sequences.
    length : int
        Alignment length in bp.
    seed : int, default=42
        RNG seed.

    Returns
    -------
    Alignment
        Synthetic alignment.
    """
    rng = random.Random(seed)
    ancestor = [rng.choice('ACGT') for _ in range(length)]
    sequences = []
    for i in range(n_seqs):
        seq = ancestor.copy()
        for _ in range(max(1, length // 20)):
            pos = rng.randrange(length)
            seq[pos] = rng.choice('ACGT')
        sequences.append(Sequence(f's{i}', ''.join(seq)))
    return Alignment(sequences)


SMALL = synthetic_alignment(20, 500)
MEDIUM = synthetic_alignment(100, 2000)


@pytest.mark.parametrize('method', ['hamming', 'jc', 'k2p'])
def test_distance_matrix_small(benchmark, method):
    """Distance matrix on a small alignment (20 x 500)."""
    benchmark(pairwise_distance_matrix, SMALL, method=method)


def test_distance_matrix_medium_hamming(benchmark):
    """Hamming distance matrix on a medium alignment (100 x 2k)."""
    benchmark(pairwise_distance_matrix, MEDIUM, method='hamming')


@pytest.mark.parametrize('name', ['mst', 'msn', 'tcs', 'mjn', 'pn', 'tsw'])
def test_algorithm_small(benchmark, name):
    """Each network algorithm on the small alignment."""
    from pypopart.algorithms import build

    params = {'random_seed': 42} if name == 'pn' else {}
    algo = build(name, **params)
    benchmark(algo.build_network, SMALL)


@pytest.mark.parametrize('name', ['mst', 'msn'])
def test_algorithm_medium(benchmark, name):
    """Fast algorithms on the medium alignment.

    TCS (~3 min: many intermediate inferences), MJN, TSW and PN are
    excluded here - they are inherently heavy at this size (as in
    PopART) and would dominate the benchmark job; track them at the
    small size instead.
    """
    from pypopart.algorithms import build

    algo = build(name)
    benchmark(algo.build_network, MEDIUM)
