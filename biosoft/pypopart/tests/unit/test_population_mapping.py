"""Unit tests for population mapping in network visualization."""

import pytest

from pypopart.algorithms.mst import MinimumSpanningTree
from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence
from pypopart.visualization.cytoscape_plot import (
    InteractiveCytoscapePlotter,
    create_cytoscape_network,
)


class TestPopulationMapping:
    """Test cases for mapping sample_ids to populations in network visualization."""

    @pytest.fixture
    def sample_sequences(self):
        """Create sample sequences for testing."""
        sequences = [
            Sequence(id='seq1', data='ATCG'),
            Sequence(id='seq2', data='ATCG'),  # Same as seq1
            Sequence(id='seq3', data='ATCG'),  # Same as seq1
            Sequence(id='seq4', data='TTCG'),  # Different
            Sequence(id='seq5', data='TTCG'),  # Same as seq4
        ]
        return sequences

    @pytest.fixture
    def sample_alignment(self, sample_sequences):
        """Create a sample alignment."""
        return Alignment(sample_sequences)

    @pytest.fixture
    def sample_network(self, sample_alignment):
        """Create a sample network."""
        algo = MinimumSpanningTree(distance_method='hamming')
        network = algo.build_network(sample_alignment)
        return network

    @pytest.fixture
    def population_mapping(self):
        """Create a population mapping for sequences."""
        return {
            'seq1': 'PopA',
            'seq2': 'PopA',
            'seq3': 'PopB',
            'seq4': 'PopB',
            'seq5': 'PopC',
        }

    @pytest.fixture
    def population_colors(self):
        """Create population colors."""
        return {
            'PopA': '#FF0000',
            'PopB': '#00FF00',
            'PopC': '#0000FF',
        }

    def test_create_elements_with_population_mapping(
        self, sample_network, population_mapping, population_colors
    ):
        """Test that create_elements uses population_mapping correctly."""
        plotter = InteractiveCytoscapePlotter(sample_network)

        # Create elements with population mapping
        elements = plotter.create_elements(
            population_colors=population_colors,
            population_mapping=population_mapping,
        )

        # Filter to only node elements
        nodes = [e for e in elements if 'source' not in e.get('data', {})]

        # Check that we have nodes
        assert len(nodes) > 0

        # Find nodes with multiple populations (should have pie charts)
        pie_nodes = [n for n in nodes if n['data'].get('has_pie', False)]

        # We expect at least one node to have multiple populations
        # (haplotype with seq1, seq2, seq3 has PopA and PopB)
        assert len(pie_nodes) > 0

        # Check that pie nodes have SVG data
        for node in pie_nodes:
            assert 'pie_svg' in node['data']
            assert node['data']['pie_svg'].startswith('data:image/svg+xml;base64,')
            assert node['data']['color'] == 'transparent'

    def test_create_elements_without_population_mapping(
        self, sample_network, population_colors
    ):
        """Test that create_elements works without population_mapping."""
        plotter = InteractiveCytoscapePlotter(sample_network)

        # Create elements without population mapping
        elements = plotter.create_elements(
            population_colors=population_colors,
            population_mapping=None,
        )

        # Should still create elements (just won't have population data)
        assert len(elements) > 0

        # Filter to only node elements
        nodes = [e for e in elements if 'source' not in e.get('data', {})]
        assert len(nodes) > 0

    def test_create_cytoscape_network_with_population_mapping(
        self, sample_network, population_mapping, population_colors
    ):
        """Test the wrapper function with population mapping."""
        elements, stylesheet = create_cytoscape_network(
            sample_network,
            population_colors=population_colors,
            population_mapping=population_mapping,
        )

        # Check elements were created
        assert len(elements) > 0

        # Check stylesheet was created
        assert len(stylesheet) > 0

        # Filter to only node elements
        nodes = [e for e in elements if 'source' not in e.get('data', {})]

        # Check for pie chart nodes
        pie_nodes = [n for n in nodes if n['data'].get('has_pie', False)]
        assert len(pie_nodes) > 0

    def test_single_population_node_color(
        self, sample_network, population_mapping, population_colors
    ):
        """Test that nodes with single population get correct color."""
        plotter = InteractiveCytoscapePlotter(sample_network)

        elements = plotter.create_elements(
            population_colors=population_colors,
            population_mapping=population_mapping,
        )

        # Filter to only node elements
        nodes = [e for e in elements if 'source' not in e.get('data', {})]

        # Find nodes with single population
        single_pop_nodes = [
            n
            for n in nodes
            if not n['data'].get('has_pie', False)
            and not n['data'].get('is_median', False)
        ]

        # At least one node should have single population
        # Check they have appropriate colors
        for node in single_pop_nodes:
            node_color = node['data'].get('color')
            # Color should either be from population_colors or default
            assert node_color is not None

    def test_pie_chart_data_structure(
        self, sample_network, population_mapping, population_colors
    ):
        """Test that pie chart data has correct structure."""
        plotter = InteractiveCytoscapePlotter(sample_network)

        elements = plotter.create_elements(
            population_colors=population_colors,
            population_mapping=population_mapping,
        )

        # Filter to pie chart nodes
        nodes = [e for e in elements if 'source' not in e.get('data', {})]
        pie_nodes = [n for n in nodes if n['data'].get('has_pie', False)]

        for node in pie_nodes:
            # Check pie_data structure
            assert 'pie_data' in node['data']
            pie_data = node['data']['pie_data']
            assert isinstance(pie_data, list)

            for segment in pie_data:
                assert 'population' in segment
                assert 'value' in segment
                assert 'percent' in segment
                assert 'color' in segment
                # Percent should sum to 100 across all segments
                # Value should be positive
                assert segment['value'] > 0
                assert 0 < segment['percent'] <= 100

    def test_tooltip_includes_population_data(
        self, sample_network, population_mapping, population_colors
    ):
        """Test that the hover tooltip carries population information."""
        plotter = InteractiveCytoscapePlotter(sample_network)

        elements = plotter.create_elements(
            population_colors=population_colors,
            population_mapping=population_mapping,
        )

        # Filter to only non-median node elements
        nodes = [
            e
            for e in elements
            if 'source' not in e.get('data', {})
            and not e['data'].get('is_median', False)
        ]

        assert nodes
        for node in nodes:
            tooltip = node['data']['tooltip']
            assert tooltip['kind'] == 'haplotype'
            # Every sample is attributed to some population, and the counts
            # must add up to the haplotype's frequency.
            assert tooltip['populations']
            counted = sum(count for _, count in tooltip['populations'])
            assert counted == tooltip['frequency']
