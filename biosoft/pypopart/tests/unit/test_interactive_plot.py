"""Unit tests for interactive network visualization."""

import os
import tempfile

import plotly.graph_objects as go
import pytest

from pypopart.core.graph import HaplotypeNetwork
from pypopart.core.haplotype import Haplotype
from pypopart.core.sequence import Sequence
from pypopart.visualization.interactive_plot import (
    InteractiveNetworkPlotter,
    create_interactive_figure,
    plot_interactive_network,
)
from pypopart.visualization.style import (
    DEFAULT_MEDIAN_COLOR,
    DEFAULT_NODE_COLOR,
)


@pytest.fixture
def simple_network():
    """Create a simple haplotype network for testing."""
    network = HaplotypeNetwork('TestNetwork')

    # Create haplotypes with sample IDs
    seq1 = Sequence('hap1', 'ATCG')
    seq2 = Sequence('hap2', 'ATCC')
    seq3 = Sequence('hap3', 'GTCG')

    # Create populations dict
    pops1 = {f's{i}': 'PopA' for i in range(7)}
    pops1.update({f's{i}': 'PopB' for i in range(7, 10)})

    pops2 = {f's{i}': 'PopB' for i in range(10, 15)}

    pops3 = {f's{i}': 'PopA' for i in range(15, 18)}

    hap1 = Haplotype(seq1, sample_ids=list(pops1.keys()), populations=pops1)
    hap2 = Haplotype(seq2, sample_ids=list(pops2.keys()), populations=pops2)
    hap3 = Haplotype(seq3, sample_ids=list(pops3.keys()), populations=pops3)

    # Add to network
    network.add_haplotype(hap1)
    network.add_haplotype(hap2)
    network.add_haplotype(hap3)

    # Add edges
    network.add_edge('hap1', 'hap2', distance=1)
    network.add_edge('hap1', 'hap3', distance=1)

    return network


@pytest.fixture
def network_with_median():
    """Create a network with median vectors."""
    network = HaplotypeNetwork('MedianNetwork')

    # Create haplotypes
    seq1 = Sequence('hap1', 'ATCG')
    seq2 = Sequence('hap2', 'GTCC')

    hap1 = Haplotype(seq1, sample_ids=[f's{i}' for i in range(5)])
    hap2 = Haplotype(seq2, sample_ids=[f's{i}' for i in range(5, 8)])

    network.add_haplotype(hap1)
    network.add_haplotype(hap2)

    # Add median vector
    median_seq = Sequence('mv1', 'ATCC')
    median_hap = Haplotype(median_seq)
    network.add_haplotype(median_hap, median_vector=True)

    # Connect through median
    network.add_edge('hap1', 'mv1', distance=1)
    network.add_edge('mv1', 'hap2', distance=1)

    return network


class TestInteractiveNetworkPlotter:
    """Test InteractiveNetworkPlotter class."""

    def test_init(self, simple_network):
        """Test plotter initialization."""
        plotter = InteractiveNetworkPlotter(simple_network)
        assert plotter.network == simple_network
        assert plotter.figure is None

    def test_plot_basic(self, simple_network):
        """Test basic plotting."""
        plotter = InteractiveNetworkPlotter(simple_network)
        fig = plotter.plot()

        assert fig is not None
        assert isinstance(fig, go.Figure)
        assert plotter.figure == fig

        # Check that figure has traces
        assert len(fig.data) > 0

    def test_plot_with_layout_algorithms(self, simple_network):
        """Test different layout algorithms."""
        plotter = InteractiveNetworkPlotter(simple_network)

        algorithms = ['spring', 'circular', 'kamada_kawai', 'spectral', 'shell']
        for algo in algorithms:
            fig = plotter.plot(layout_algorithm=algo)
            assert fig is not None
            assert isinstance(fig, go.Figure)

    def test_plot_with_invalid_layout(self, simple_network):
        """Test plot with invalid layout algorithm."""
        plotter = InteractiveNetworkPlotter(simple_network)

        with pytest.raises(ValueError, match='Unknown layout algorithm'):
            plotter.plot(layout_algorithm='invalid')

    def test_plot_with_custom_layout(self, simple_network):
        """Test plot with pre-computed layout."""
        plotter = InteractiveNetworkPlotter(simple_network)

        # Custom layout positions
        layout = {'hap1': (0, 0), 'hap2': (1, 0), 'hap3': (0.5, 1)}

        fig = plotter.plot(layout=layout)
        assert fig is not None

    def test_plot_customization(self, simple_network):
        """Test plot with customization options."""
        plotter = InteractiveNetworkPlotter(simple_network)

        fig = plotter.plot(
            node_size_scale=30,
            edge_width_scale=3.0,
            show_labels=False,
            width=1200,
            height=900,
            title='Custom Title',
        )

        assert fig is not None
        assert 'Custom Title' in fig.layout.title.text
        assert fig.layout.width == 1200
        assert fig.layout.height == 900

    def test_plot_with_population_colors(self, simple_network):
        """Test plot with population coloring."""
        plotter = InteractiveNetworkPlotter(simple_network)

        population_colors = {'PopA': 'red', 'PopB': 'blue'}

        fig = plotter.plot(population_colors=population_colors)
        assert fig is not None

    def test_plot_with_custom_node_colors(self, simple_network):
        """Test plot with custom node colors."""
        plotter = InteractiveNetworkPlotter(simple_network)

        node_colors = {'hap1': 'red', 'hap2': 'green', 'hap3': 'blue'}

        fig = plotter.plot(node_color_map=node_colors)
        assert fig is not None

    def test_plot_median_vectors(self, network_with_median):
        """Test plotting network with median vectors."""
        plotter = InteractiveNetworkPlotter(network_with_median)

        fig = plotter.plot(median_vector_color='yellow')

        assert fig is not None

    def test_add_population_legend_without_plot(self, simple_network):
        """Test adding legend without plot raises error."""
        plotter = InteractiveNetworkPlotter(simple_network)

        with pytest.raises(ValueError, match='No plot exists'):
            plotter.add_population_legend({'PopA': 'red'})

    def test_add_population_legend(self, simple_network):
        """Test adding population legend."""
        plotter = InteractiveNetworkPlotter(simple_network)
        fig = plotter.plot()

        initial_traces = len(fig.data)

        pop_colors = {'PopA': 'red', 'PopB': 'blue'}
        plotter.add_population_legend(pop_colors)

        # Should have added legend traces
        assert len(fig.data) > initial_traces

    def test_save_html_without_plot(self, simple_network):
        """Test saving without plot raises error."""
        plotter = InteractiveNetworkPlotter(simple_network)

        with pytest.raises(ValueError, match='No plot exists'):
            plotter.save_html('test.html')

    def test_save_html(self, simple_network):
        """Test saving plot as HTML."""
        plotter = InteractiveNetworkPlotter(simple_network)
        plotter.plot()

        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(tmpdir, 'test.html')
            plotter.save_html(filename)

            assert os.path.exists(filename)
            assert os.path.getsize(filename) > 0

            # Check file contains expected content
            with open(filename, 'r') as f:
                content = f.read()
                assert 'plotly' in content.lower()

    def test_show_without_plot(self, simple_network):
        """Test show without plot raises error."""
        plotter = InteractiveNetworkPlotter(simple_network)

        with pytest.raises(ValueError, match='No plot exists'):
            plotter.show()

    def test_get_node_color_median(self, network_with_median):
        """Test node color for median vectors."""
        plotter = InteractiveNetworkPlotter(network_with_median)

        color = plotter._get_node_color('mv1', None, None, None, is_median=True)
        assert color == DEFAULT_MEDIAN_COLOR

    def test_get_node_color_custom(self, simple_network):
        """Test node color with custom mapping."""
        plotter = InteractiveNetworkPlotter(simple_network)

        node_colors = {'hap1': 'red'}
        hap = simple_network.get_haplotype('hap1')

        color = plotter._get_node_color('hap1', hap, node_colors, None, False)
        assert color == 'red'

    def test_get_node_color_population(self, simple_network):
        """Test node color with population mapping."""
        plotter = InteractiveNetworkPlotter(simple_network)

        pop_colors = {'PopA': 'red', 'PopB': 'blue'}
        hap = simple_network.get_haplotype('hap1')

        color = plotter._get_node_color('hap1', hap, None, pop_colors, False)
        assert color in ['red', 'blue']  # Should be one of the population colors

    def test_get_node_color_default(self, simple_network):
        """Test default node color."""
        plotter = InteractiveNetworkPlotter(simple_network)

        hap = simple_network.get_haplotype('hap1')
        color = plotter._get_node_color('hap1', hap, None, None, False)
        assert color == DEFAULT_NODE_COLOR

    def test_create_hover_text_haplotype(self, simple_network):
        """Test hover text creation for haplotype."""
        plotter = InteractiveNetworkPlotter(simple_network)

        hap = simple_network.get_haplotype('hap1')
        hover_text = plotter._create_hover_text('hap1', hap, False)

        assert 'hap1' in hover_text
        assert 'Frequency' in hover_text
        # Sequence should NOT be in hover text (removed per requirements)
        assert 'Sequence' not in hover_text
        assert 'ATCG' not in hover_text

    def test_create_hover_text_median(self, network_with_median):
        """Test hover text creation for median vector."""
        plotter = InteractiveNetworkPlotter(network_with_median)

        mv_hap = network_with_median.get_haplotype('mv1')
        hover_text = plotter._create_hover_text('mv1', mv_hap, True)

        assert 'mv1' in hover_text
        assert 'Median Vector' in hover_text

    def test_create_hover_text_with_populations(self, simple_network):
        """Test hover text with population information."""
        plotter = InteractiveNetworkPlotter(simple_network)

        hap = simple_network.get_haplotype('hap1')
        hover_text = plotter._create_hover_text('hap1', hap, False)

        assert 'Populations' in hover_text or 'Population' in hover_text

    def test_create_hover_text_many_samples(self):
        """Test hover text with many samples (>10)."""
        network = HaplotypeNetwork('ManyS samples')
        seq = Sequence('hap1', 'ATCG')

        # Create haplotype with many samples
        sample_ids = [f's{i}' for i in range(20)]
        hap = Haplotype(seq, sample_ids=sample_ids)
        network.add_haplotype(hap)

        plotter = InteractiveNetworkPlotter(network)
        hover_text = plotter._create_hover_text('hap1', hap, False)

        assert '20 total' in hover_text


class TestConvenienceFunctions:
    """Test convenience functions."""

    def test_plot_interactive_network(self, simple_network):
        """Test plot_interactive_network convenience function."""
        fig = plot_interactive_network(simple_network)

        assert fig is not None
        assert isinstance(fig, go.Figure)

    def test_plot_interactive_network_with_options(self, simple_network):
        """Test plot_interactive_network with options."""
        fig = plot_interactive_network(
            simple_network,
            layout_algorithm='circular',
            show_labels=False,
            title='Test Network',
        )

        assert fig is not None
        assert 'Test Network' in fig.layout.title.text

    def test_create_interactive_figure(self, simple_network):
        """Test create_interactive_figure."""
        pop_colors = {'PopA': 'red', 'PopB': 'blue'}

        fig = create_interactive_figure(simple_network, population_colors=pop_colors)

        assert fig is not None
        assert isinstance(fig, go.Figure)

    def test_create_interactive_figure_save(self, simple_network):
        """Test create_interactive_figure with save."""
        pop_colors = {'PopA': 'red', 'PopB': 'blue'}

        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(tmpdir, 'interactive.html')

            create_interactive_figure(
                simple_network,
                population_colors=pop_colors,
                filename=filename,
                auto_open=False,
            )

            assert os.path.exists(filename)
            assert os.path.getsize(filename) > 0

    def test_create_interactive_figure_no_populations(self, simple_network):
        """Test interactive figure without population colors."""
        fig = create_interactive_figure(simple_network)

        assert fig is not None


class TestEmptyNetwork:
    """Test plotting with empty or minimal networks."""

    def test_plot_empty_network(self):
        """Test plotting empty network."""
        network = HaplotypeNetwork('Empty')
        plotter = InteractiveNetworkPlotter(network)

        fig = plotter.plot()
        assert fig is not None

    def test_plot_single_node(self):
        """Test plotting network with single node."""
        network = HaplotypeNetwork('SingleNode')
        seq = Sequence('hap1', 'ATCG')
        hap = Haplotype(seq, sample_ids=[f's{i}' for i in range(5)])
        network.add_haplotype(hap)

        plotter = InteractiveNetworkPlotter(network)
        fig = plotter.plot()

        assert fig is not None

    def test_plot_disconnected_nodes(self):
        """Test plotting network with disconnected components."""
        network = HaplotypeNetwork('Disconnected')

        seq1 = Sequence('hap1', 'ATCG')
        seq2 = Sequence('hap2', 'GTCC')

        hap1 = Haplotype(seq1, sample_ids=[f's{i}' for i in range(3)])
        hap2 = Haplotype(seq2, sample_ids=[f's{i}' for i in range(3, 5)])

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        # No edges - disconnected

        plotter = InteractiveNetworkPlotter(network)
        fig = plotter.plot()

        assert fig is not None


class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_plot_with_zero_frequency(self):
        """Test plotting haplotype with zero frequency."""
        network = HaplotypeNetwork('ZeroFreq')
        seq = Sequence('hap1', 'ATCG')
        hap = Haplotype(seq)  # No sample IDs = frequency 0
        network.add_haplotype(hap)

        plotter = InteractiveNetworkPlotter(network)
        fig = plotter.plot()

        assert fig is not None

    def test_node_without_haplotype(self, simple_network):
        """Test handling of nodes without haplotype data."""
        plotter = InteractiveNetworkPlotter(simple_network)

        # This should not crash
        color = plotter._get_node_color('unknown', None, None, None, False)
        assert color == DEFAULT_NODE_COLOR

        hover_text = plotter._create_hover_text('unknown', None, False)
        assert 'unknown' in hover_text

    def test_multiple_plots_same_plotter(self, simple_network):
        """Test creating multiple plots with same plotter."""
        plotter = InteractiveNetworkPlotter(simple_network)

        fig1 = plotter.plot()
        fig2 = plotter.plot(layout_algorithm='circular')

        assert fig1 != fig2

    def test_plot_with_edge_without_weight(self):
        """Test plotting edges that have no weight attribute."""
        network = HaplotypeNetwork('NoWeight')

        seq1 = Sequence('hap1', 'ATCG')
        seq2 = Sequence('hap2', 'ATCC')

        hap1 = Haplotype(seq1, sample_ids=[f's{i}' for i in range(5)])
        hap2 = Haplotype(seq2, sample_ids=[f's{i}' for i in range(5, 8)])

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        # Add edge without explicit weight
        network._graph.add_edge('hap1', 'hap2')

        plotter = InteractiveNetworkPlotter(network)
        fig = plotter.plot()

        assert fig is not None

    def test_layout_options(self, simple_network):
        """Test passing layout options."""
        plotter = InteractiveNetworkPlotter(simple_network)

        fig = plotter.plot(
            margin={'l': 20, 'r': 20, 't': 40, 'b': 20}, paper_bgcolor='lightgray'
        )

        assert fig is not None
        assert fig.layout.paper_bgcolor == 'lightgray'
