"""Tests for the NetworkAlgorithm base class parameter handling."""

import pytest

from pypopart.algorithms import (
    TCS,
    MedianJoiningNetwork,
    MinimumSpanningNetwork,
    MinimumSpanningTree,
    ParsimonyNetwork,
    TightSpanWalker,
)

ALL_ALGORITHMS = [
    MinimumSpanningTree,
    MinimumSpanningNetwork,
    TCS,
    MedianJoiningNetwork,
    ParsimonyNetwork,
    TightSpanWalker,
]


@pytest.mark.parametrize('algorithm_cls', ALL_ALGORITHMS)
def test_unknown_kwarg_raises(algorithm_cls):
    """A misspelled parameter must raise instead of being silently ignored."""
    with pytest.raises(TypeError, match='distance_metric'):
        algorithm_cls(distance_metric='hamming')


@pytest.mark.parametrize('algorithm_cls', ALL_ALGORITHMS)
def test_valid_construction(algorithm_cls):
    """Standard construction with distance_method and ignore_gaps works."""
    algo = algorithm_cls(distance_method='hamming', ignore_gaps=False)
    assert algo.distance_method == 'hamming'
    assert algo.params['ignore_gaps'] is False
