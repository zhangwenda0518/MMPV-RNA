"""Tests for figure labelling: node labels and mutation counts."""

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pytest

from pypopart.core.graph import HaplotypeNetwork
from pypopart.core.haplotype import Haplotype
from pypopart.core.sequence import Sequence
from pypopart.visualization.static_plot import (
    TICK_SPAN_FRACTION,
    StaticNetworkPlotter,
    tick_offsets,
)


@pytest.fixture
def network():
    """Build a two-haplotype network four mutations apart."""
    net = HaplotypeNetwork(name='LabelTest')
    net.add_haplotype(Haplotype(Sequence('Alpha_01', 'ATCG'), sample_ids=['a']))
    net.add_haplotype(Haplotype(Sequence('Beta_01', 'TTCC'), sample_ids=['b']))
    net.add_edge('Alpha_01', 'Beta_01', distance=4)
    return net


def _texts(ax):
    """Collect the rendered text of every label on the axes."""
    return {t.get_text() for t in ax.texts}


class TestNodeLabels:
    """Nodes carry the same labels the interactive view shows."""

    def test_defaults_to_node_ids(self, network):
        """Without a mapping the node's own ID is drawn."""
        _, ax = StaticNetworkPlotter(network).plot(show_mutations=False)
        try:
            assert {'Alpha_01', 'Beta_01'} <= _texts(ax)
        finally:
            plt.close('all')

    def test_uses_supplied_labels(self, network):
        """The GUI passes H numbers so the figure matches the screen."""
        _, ax = StaticNetworkPlotter(network).plot(
            node_labels={'Alpha_01': 'H1', 'Beta_01': 'H2'},
            show_mutations=False,
        )
        try:
            drawn = _texts(ax)
            assert {'H1', 'H2'} <= drawn
            assert 'Alpha_01' not in drawn
        finally:
            plt.close('all')


class TestEdgeMutationCounts:
    """Edges report the real number of mutations."""

    def test_numeral_shows_the_distance(self, network):
        """'weight' defaults to 1.0 for every edge; 'distance' is the count."""
        _, ax = StaticNetworkPlotter(network).plot(
            show_labels=False, show_edge_ticks=False
        )
        try:
            assert '4' in _texts(ax)
            # The old code read weight and labelled every edge '1'.
            assert '1' not in _texts(ax)
        finally:
            plt.close('all')

    def test_ticks_replace_the_numeral(self, network):
        """Below the threshold the count is drawn as strokes, not text."""
        _, ax = StaticNetworkPlotter(network).plot(
            show_labels=False, show_edge_ticks=True, edge_tick_threshold=10
        )
        try:
            assert '4' not in _texts(ax)
            # One line per mutation, on top of the edge itself.
            assert len(ax.lines) == 4
        finally:
            plt.close('all')

    def test_long_edges_fall_back_to_a_numeral(self, network):
        """A comb of 40 strokes is unreadable, so past the threshold: text."""
        network.graph['Alpha_01']['Beta_01']['distance'] = 40
        _, ax = StaticNetworkPlotter(network).plot(
            show_labels=False, show_edge_ticks=True, edge_tick_threshold=10
        )
        try:
            assert '40' in _texts(ax)
            assert not ax.lines
        finally:
            plt.close('all')

    def test_zero_distance_is_not_marked(self, network):
        """Identical haplotypes need no mutation marker."""
        network.graph['Alpha_01']['Beta_01']['distance'] = 0
        _, ax = StaticNetworkPlotter(network).plot(show_labels=False)
        try:
            assert not ax.lines
            assert not _texts(ax) - {'LabelTest'}
        finally:
            plt.close('all')


class TestTickOffsets:
    """Ticks sit in the gap between the discs and tighten on short edges."""

    def test_single_tick_centres_on_the_visible_span(self):
        """A bigger disc at one end pushes the tick towards the other."""
        assert tick_offsets(10.0, 4.0, 0.0, 1, 0.1) == [2.0]
        assert tick_offsets(10.0, 1.0, 1.0, 1, 0.1) == [0.0]

    def test_comb_spans_a_fraction_of_the_visible_edge(self):
        """Ticks spread over the configured share of the gap, symmetric."""
        offsets = tick_offsets(20.0, 2.0, 2.0, 5, 0.1)

        assert offsets is not None
        assert offsets[-1] - offsets[0] == pytest.approx(16.0 * TICK_SPAN_FRACTION)
        assert offsets[0] == pytest.approx(-offsets[-1])

    def test_spacing_shrinks_with_the_edge(self):
        """Closer nodes, tighter comb."""
        wide = tick_offsets(20.0, 2.0, 2.0, 4, 0.1)
        narrow = tick_offsets(10.0, 2.0, 2.0, 4, 0.1)

        assert wide[1] - wide[0] > narrow[1] - narrow[0]

    @pytest.mark.parametrize(
        ('length', 'ra', 'rb', 'count', 'min_spacing'),
        [
            (10.0, 2.0, 2.0, 4, 2.0),  # 6 visible, 3 gaps of 1.2 < 2
            (4.0, 2.0, 2.0, 1, 0.1),  # discs touch
            (3.0, 2.0, 2.0, 1, 0.1),  # discs overlap
            (10.0, 1.0, 1.0, 0, 0.1),  # nothing to draw
        ],
    )
    def test_unfit_returns_none(self, length, ra, rb, count, min_spacing):
        """None tells the caller to draw a numeral."""
        assert tick_offsets(length, ra, rb, count, min_spacing) is None


class TestTicksClearTheNodes:
    """On the page, no tick may be hidden under a node disc."""

    def test_ticks_lie_outside_both_discs(self, network):
        """Every stroke endpoint is further than the radius from each centre."""
        layout = {'Alpha_01': (0.0, 0.0), 'Beta_01': (1.0, 0.0)}
        plotter = StaticNetworkPlotter(network)
        _, ax = plotter.plot(
            layout=layout, show_labels=False, node_size_scale=2000.0, figsize=(4, 4)
        )
        try:
            assert len(ax.lines) == 4
            radii = {
                node: plotter._points_to_data_units((size**0.5) / 2)
                for node, size in plotter._node_sizes.items()
            }
            for line in ax.lines:
                for x, y in zip(*line.get_data()):
                    for node, (cx, cy) in layout.items():
                        assert ((x - cx) ** 2 + (y - cy) ** 2) ** 0.5 > radii[node]
        finally:
            plt.close('all')

    def test_too_short_an_edge_falls_back_to_a_numeral(self, network):
        """Ten ticks cannot fit between two discs that nearly touch."""
        # A far-off third node sets the scale, so the close pair really
        # is close on the page rather than autoscaled to fill it.
        network.add_haplotype(Haplotype(Sequence('Gamma_01', 'GGGG'), sample_ids=['c']))
        network.add_edge('Beta_01', 'Gamma_01', distance=2)
        network.graph['Alpha_01']['Beta_01']['distance'] = 10
        layout = {
            'Alpha_01': (0.0, 0.0),
            'Beta_01': (0.05, 0.0),
            'Gamma_01': (10.0, 0.0),
        }
        _, ax = StaticNetworkPlotter(network).plot(
            layout=layout, show_labels=False, node_size_scale=3000.0, figsize=(4, 4)
        )
        try:
            # Only the long Beta-Gamma edge draws its two ticks.
            assert len(ax.lines) == 2
            assert '10' in _texts(ax)
        finally:
            plt.close('all')
