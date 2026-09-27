"""Unit tests for static network visualization."""

import matplotlib
import pytest

matplotlib.use('Agg')  # Use non-interactive backend for testing
import os
import tempfile

import matplotlib.pyplot as plt

from pypopart.core.graph import HaplotypeNetwork
from pypopart.core.haplotype import Haplotype
from pypopart.core.sequence import Sequence
from pypopart.visualization.static_plot import (
    StaticNetworkPlotter,
    create_publication_figure,
    plot_network,
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


class TestStaticNetworkPlotter:
    """Test StaticNetworkPlotter class."""

    def test_init(self, simple_network):
        """Test plotter initialization."""
        plotter = StaticNetworkPlotter(simple_network)
        assert plotter.network == simple_network
        assert plotter.figure is None
        assert plotter.ax is None

    def test_plot_basic(self, simple_network):
        """Test basic plotting."""
        plotter = StaticNetworkPlotter(simple_network)
        fig, ax = plotter.plot()

        assert fig is not None
        assert ax is not None
        assert plotter.figure == fig
        assert plotter.ax == ax

        plt.close(fig)

    def test_plot_with_layout_algorithms(self, simple_network):
        """Test different layout algorithms."""
        plotter = StaticNetworkPlotter(simple_network)

        algorithms = ['spring', 'circular', 'kamada_kawai', 'spectral', 'shell']
        for algo in algorithms:
            fig, ax = plotter.plot(layout_algorithm=algo)
            assert fig is not None
            plt.close(fig)

    def test_plot_with_invalid_layout(self, simple_network):
        """Test plot with invalid layout algorithm."""
        plotter = StaticNetworkPlotter(simple_network)

        with pytest.raises(ValueError, match='Unknown layout algorithm'):
            plotter.plot(layout_algorithm='invalid')

    def test_plot_with_custom_layout(self, simple_network):
        """Test plot with pre-computed layout."""
        plotter = StaticNetworkPlotter(simple_network)

        # Custom layout positions
        layout = {'hap1': (0, 0), 'hap2': (1, 0), 'hap3': (0.5, 1)}

        fig, ax = plotter.plot(layout=layout)
        assert fig is not None
        plt.close(fig)

    def test_plot_customization(self, simple_network):
        """Test plot with customization options."""
        plotter = StaticNetworkPlotter(simple_network)

        fig, ax = plotter.plot(
            node_size_scale=500,
            edge_width_scale=2.0,
            show_labels=False,
            show_mutations=False,
            figsize=(15, 12),
            title='Custom Title',
        )

        assert fig is not None
        assert ax.get_title() == 'Custom Title'
        plt.close(fig)

    def test_plot_with_population_colors(self, simple_network):
        """Test plot with population coloring."""
        plotter = StaticNetworkPlotter(simple_network)

        population_colors = {'PopA': 'red', 'PopB': 'blue'}

        fig, ax = plotter.plot(population_colors=population_colors)
        assert fig is not None
        plt.close(fig)

    def test_plot_with_custom_node_colors(self, simple_network):
        """Test plot with custom node colors."""
        plotter = StaticNetworkPlotter(simple_network)

        node_colors = {'hap1': 'red', 'hap2': 'green', 'hap3': 'blue'}

        fig, ax = plotter.plot(node_color_map=node_colors)
        assert fig is not None
        plt.close(fig)

    def test_plot_median_vectors(self, network_with_median):
        """Test plotting network with median vectors."""
        plotter = StaticNetworkPlotter(network_with_median)

        fig, ax = plotter.plot(
            median_vector_color='yellow',
            median_vector_marker='D',  # Diamond
        )

        assert fig is not None
        plt.close(fig)

    def test_compute_node_sizes(self, simple_network):
        """Test node size computation."""
        plotter = StaticNetworkPlotter(simple_network)

        sizes = plotter._compute_node_sizes(scale=300)

        assert 'hap1' in sizes
        assert 'hap2' in sizes
        assert 'hap3' in sizes
        assert sizes['hap1'] > sizes['hap2']  # More frequent = larger
        assert sizes['hap2'] > sizes['hap3']

    def test_compute_node_colors(self, simple_network):
        """Test node color computation."""
        plotter = StaticNetworkPlotter(simple_network)

        # Test with population colors
        pop_colors = {'PopA': 'red', 'PopB': 'blue'}
        colors = plotter._compute_node_colors(None, pop_colors, 'gray')

        assert 'hap1' in colors
        assert colors['hap1'] == 'red'  # Dominant population is PopA
        assert colors['hap2'] == 'blue'  # Only PopB

    def test_compute_node_colors_median(self, network_with_median):
        """Test median vector coloring."""
        plotter = StaticNetworkPlotter(network_with_median)

        colors = plotter._compute_node_colors(None, None, 'lightgray')

        assert colors['mv1'] == 'lightgray'

    def test_compute_edge_widths(self, simple_network):
        """Test edge width computation."""
        plotter = StaticNetworkPlotter(simple_network)

        widths = plotter._compute_edge_widths(scale=1.0)

        assert len(widths) == 2  # Two edges
        assert all(w > 0 for w in widths)

    def test_add_legend_without_plot(self, simple_network):
        """Test adding legend without plot raises error."""
        plotter = StaticNetworkPlotter(simple_network)

        with pytest.raises(ValueError, match='No plot exists'):
            plotter.add_legend()

    def test_add_legend(self, simple_network):
        """Test adding legend to plot."""
        plotter = StaticNetworkPlotter(simple_network)
        fig, ax = plotter.plot()

        pop_colors = {'PopA': 'red', 'PopB': 'blue'}
        plotter.add_legend(
            population_colors=pop_colors, show_median_vectors=True, show_size_scale=True
        )

        # Check that legend was added
        legend = ax.get_legend()
        assert legend is not None

        plt.close(fig)

    def test_add_legend_minimal(self, simple_network):
        """Test adding minimal legend."""
        plotter = StaticNetworkPlotter(simple_network)
        fig, ax = plotter.plot()

        plotter.add_legend(show_median_vectors=False, show_size_scale=False)

        plt.close(fig)

    def test_add_scale_bar_without_plot(self, simple_network):
        """Test adding scale bar without plot raises error."""
        plotter = StaticNetworkPlotter(simple_network)

        with pytest.raises(ValueError, match='No plot exists'):
            plotter.add_scale_bar()

    def test_add_scale_bar(self, simple_network):
        """Test adding scale bar to plot."""
        plotter = StaticNetworkPlotter(simple_network)
        fig, ax = plotter.plot()

        plotter.add_scale_bar(num_mutations=1)
        plotter.add_scale_bar(num_mutations=5, position=(0.8, 0.1))

        plt.close(fig)

    def test_add_statistics_annotation_without_plot(self, simple_network):
        """Test adding statistics without plot raises error."""
        plotter = StaticNetworkPlotter(simple_network)

        with pytest.raises(ValueError, match='No plot exists'):
            plotter.add_statistics_annotation()

    def test_add_statistics_annotation(self, simple_network):
        """Test adding statistics annotation."""
        plotter = StaticNetworkPlotter(simple_network)
        fig, ax = plotter.plot()

        plotter.add_statistics_annotation()

        plt.close(fig)

    def test_add_statistics_annotation_custom(self, simple_network):
        """Test adding custom statistics."""
        plotter = StaticNetworkPlotter(simple_network)
        fig, ax = plotter.plot()

        custom_stats = {'Metric A': 42, 'Metric B': 3.14159, 'Metric C': 'value'}

        plotter.add_statistics_annotation(stats=custom_stats)

        plt.close(fig)

    def test_save_without_plot(self, simple_network):
        """Test saving without plot raises error."""
        plotter = StaticNetworkPlotter(simple_network)

        with pytest.raises(ValueError, match='No plot exists'):
            plotter.save('test.png')

    def test_save_png(self, simple_network):
        """Test saving plot as PNG."""
        plotter = StaticNetworkPlotter(simple_network)
        fig, ax = plotter.plot()

        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(tmpdir, 'test.png')
            plotter.save(filename, dpi=150)

            assert os.path.exists(filename)
            assert os.path.getsize(filename) > 0

        plt.close(fig)

    def test_save_pdf(self, simple_network):
        """Test saving plot as PDF."""
        plotter = StaticNetworkPlotter(simple_network)
        fig, ax = plotter.plot()

        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(tmpdir, 'test.pdf')
            plotter.save(filename)

            assert os.path.exists(filename)
            assert os.path.getsize(filename) > 0

        plt.close(fig)

    def test_save_svg(self, simple_network):
        """Test saving plot as SVG."""
        plotter = StaticNetworkPlotter(simple_network)
        fig, ax = plotter.plot()

        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(tmpdir, 'test.svg')
            plotter.save(filename)

            assert os.path.exists(filename)
            assert os.path.getsize(filename) > 0

        plt.close(fig)


class TestConvenienceFunctions:
    """Test convenience functions."""

    def test_plot_network(self, simple_network):
        """Test plot_network convenience function."""
        fig, ax = plot_network(simple_network)

        assert fig is not None
        assert ax is not None

        plt.close(fig)

    def test_plot_network_with_options(self, simple_network):
        """Test plot_network with options."""
        fig, ax = plot_network(
            simple_network,
            layout_algorithm='circular',
            show_labels=False,
            title='Test Network',
        )

        assert fig is not None
        assert ax.get_title() == 'Test Network'

        plt.close(fig)

    def test_create_publication_figure(self, simple_network):
        """Test create_publication_figure."""
        pop_colors = {'PopA': 'red', 'PopB': 'blue'}

        fig, ax = create_publication_figure(
            simple_network, population_colors=pop_colors
        )

        assert fig is not None
        assert ax is not None

        # Check that legend was added
        legend = ax.get_legend()
        assert legend is not None

        plt.close(fig)

    def test_create_publication_figure_save(self, simple_network):
        """Test create_publication_figure with save."""
        pop_colors = {'PopA': 'red', 'PopB': 'blue'}

        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(tmpdir, 'publication.pdf')

            fig, ax = create_publication_figure(
                simple_network, population_colors=pop_colors, filename=filename
            )

            assert os.path.exists(filename)
            assert os.path.getsize(filename) > 0

            plt.close(fig)

    def test_create_publication_figure_no_populations(self, simple_network):
        """Test publication figure without population colors."""
        fig, ax = create_publication_figure(simple_network)

        assert fig is not None

        plt.close(fig)


class TestEmptyNetwork:
    """Test plotting with empty or minimal networks."""

    def test_plot_empty_network(self):
        """Test plotting empty network."""
        network = HaplotypeNetwork('Empty')
        plotter = StaticNetworkPlotter(network)

        fig, ax = plotter.plot()
        assert fig is not None

        plt.close(fig)

    def test_plot_single_node(self):
        """Test plotting network with single node."""
        network = HaplotypeNetwork('SingleNode')
        seq = Sequence('hap1', 'ATCG')
        hap = Haplotype(seq, sample_ids=[f's{i}' for i in range(5)])
        network.add_haplotype(hap)

        plotter = StaticNetworkPlotter(network)
        fig, ax = plotter.plot()

        assert fig is not None
        plt.close(fig)

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

        plotter = StaticNetworkPlotter(network)
        fig, ax = plotter.plot()

        assert fig is not None
        plt.close(fig)


class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_plot_with_zero_frequency(self):
        """Test plotting haplotype with zero frequency."""
        network = HaplotypeNetwork('ZeroFreq')
        seq = Sequence('hap1', 'ATCG')
        hap = Haplotype(seq)  # No sample IDs = frequency 0
        network.add_haplotype(hap)

        plotter = StaticNetworkPlotter(network)
        fig, ax = plotter.plot()

        assert fig is not None
        plt.close(fig)

    def test_node_without_haplotype(self, simple_network):
        """Test handling of nodes without haplotype data."""
        plotter = StaticNetworkPlotter(simple_network)

        # This should not crash even if a node has no haplotype
        sizes = plotter._compute_node_sizes(300)
        colors = plotter._compute_node_colors(None, None, 'gray')

        assert len(sizes) > 0
        assert len(colors) > 0

    def test_multiple_plots_same_plotter(self, simple_network):
        """Test creating multiple plots with same plotter."""
        plotter = StaticNetworkPlotter(simple_network)

        fig1, ax1 = plotter.plot()
        plt.close(fig1)

        fig2, ax2 = plotter.plot(layout_algorithm='circular')
        plt.close(fig2)

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

        plotter = StaticNetworkPlotter(network)
        fig, ax = plotter.plot(show_mutations=True)

        assert fig is not None
        plt.close(fig)
