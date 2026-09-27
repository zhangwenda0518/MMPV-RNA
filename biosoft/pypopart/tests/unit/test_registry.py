"""Tests for the algorithm registry."""

import pytest

from pypopart.algorithms import ALGORITHMS, build, list_algorithms
from pypopart.algorithms.base import NetworkAlgorithm


class TestRegistry:
    """The registry is the single source of truth for algorithms."""

    def test_all_six_algorithms_registered(self):
        """All six network algorithms are present."""
        assert set(ALGORITHMS) == {'mst', 'msn', 'tcs', 'mjn', 'pn', 'tsw'}

    @pytest.mark.parametrize('name', ['mst', 'msn', 'tcs', 'mjn', 'pn', 'tsw'])
    def test_build_each(self, name):
        """Every registered algorithm builds with defaults."""
        algo = build(name)
        assert isinstance(algo, NetworkAlgorithm)
        assert algo.distance_method == 'hamming'

    def test_build_case_insensitive(self):
        """Algorithm names are case-insensitive."""
        assert type(build('MST')) is ALGORITHMS['mst']

    def test_build_with_params(self):
        """Algorithm-specific parameters pass through."""
        algo = build('msn', distance_method='k2p', epsilon=1.0)
        assert algo.distance_method == 'k2p'
        assert algo.epsilon == 1.0

    def test_build_unknown_name(self):
        """Unknown names raise with the available list."""
        with pytest.raises(ValueError, match='Unknown algorithm'):
            build('bogus')

    def test_build_unknown_param(self):
        """Typoed parameters raise TypeError from the constructor."""
        with pytest.raises(TypeError, match='not_a_param'):
            build('mst', not_a_param=1)

    def test_list_algorithms(self):
        """Listing returns name/description pairs for every algorithm."""
        listing = list_algorithms()
        assert [entry['name'] for entry in listing] == list(ALGORITHMS)
        assert all(entry['description'] for entry in listing)
