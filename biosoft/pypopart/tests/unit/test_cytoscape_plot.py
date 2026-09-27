"""Tests for Cytoscape-based interactive network visualization."""

import json

import pytest

from pypopart.core.alignment import Sequence
from pypopart.core.graph import HaplotypeNetwork
from pypopart.core.haplotype import Haplotype
from pypopart.visualization.cytoscape_plot import (
    DEFAULT_TICK_THRESHOLD,
    MAX_TICK_MARKS,
    MAX_TOOLTIP_SAMPLES,
    MIN_TICK_FONT_SIZE,
    TICK_FONT_SIZE,
    InteractiveCytoscapePlotter,
    build_node_tooltip,
    create_cytoscape_network,
    create_edge_tick_stylesheet,
    fit_tick_font_size,
    format_edge_ticks,
    resolve_population_counts,
)


@pytest.fixture
def simple_network():
    """Create a simple test network with 3 haplotypes."""
    network = HaplotypeNetwork(name='TestNetwork')

    # Create sequences
    seq1 = Sequence('H1', 'ATCG')
    seq2 = Sequence('H2', 'ATGG')
    seq3 = Sequence('H3', 'GTCG')

    # Create haplotypes
    h1 = Haplotype(
        sequence=seq1,
        sample_ids=['S1', 'S2'],
        populations={'S1': 'PopA', 'S2': 'PopA'},
    )
    h2 = Haplotype(sequence=seq2, sample_ids=['S3'], populations={'S3': 'PopA'})
    h3 = Haplotype(
        sequence=seq3,
        sample_ids=['S4', 'S5', 'S6'],
        populations={'S4': 'PopB', 'S5': 'PopB', 'S6': 'PopB'},
    )

    # Add to network
    network.add_haplotype(h1)
    network.add_haplotype(h2)
    network.add_haplotype(h3)

    # Add edges
    network.add_edge('H1', 'H2', distance=1)
    network.add_edge('H2', 'H3', distance=2)

    return network


@pytest.fixture
def network_with_populations():
    """Create a network with population data."""
    network = HaplotypeNetwork(name='PopNetwork')

    # Create sequences
    seq1 = Sequence('H1', 'ATCG')
    seq2 = Sequence('H2', 'ATGG')
    seq3 = Sequence('H3', 'GTCG')

    # Create haplotypes with population data
    h1 = Haplotype(
        sequence=seq1,
        sample_ids=['S1', 'S2'],
        populations={'S1': 'PopA', 'S2': 'PopA'},
    )
    h2 = Haplotype(
        sequence=seq2,
        sample_ids=['S3', 'S4'],
        populations={'S3': 'PopA', 'S4': 'PopB'},
    )
    h3 = Haplotype(
        sequence=seq3,
        sample_ids=['S5', 'S6'],
        populations={'S5': 'PopB', 'S6': 'PopB'},
    )

    network.add_haplotype(h1)
    network.add_haplotype(h2)
    network.add_haplotype(h3)

    network.add_edge('H1', 'H2', distance=1)
    network.add_edge('H2', 'H3', distance=1)

    return network


class TestInteractiveCytoscapePlotter:
    """Tests for InteractiveCytoscapePlotter class."""

    def test_init(self, simple_network):
        """Test plotter initialization."""
        plotter = InteractiveCytoscapePlotter(simple_network)
        assert plotter.network == simple_network
        assert plotter.elements is None
        assert plotter.stylesheet is None

    def test_create_elements_basic(self, simple_network):
        """Test creating basic Cytoscape elements."""
        plotter = InteractiveCytoscapePlotter(simple_network)
        layout = {'H1': (0, 0), 'H2': (1, 0), 'H3': (2, 0)}

        elements = plotter.create_elements(layout=layout)

        # Should have 3 nodes + 2 edges = 5 elements
        assert len(elements) == 5

        # Count nodes and edges
        nodes = [e for e in elements if 'source' not in e.get('data', {})]
        edges = [e for e in elements if 'source' in e.get('data', {})]

        assert len(nodes) == 3
        assert len(edges) == 2

    def test_create_elements_with_positions(self, simple_network):
        """Test that positions are properly scaled."""
        plotter = InteractiveCytoscapePlotter(simple_network)
        layout = {'H1': (0.5, 0.5), 'H2': (1.0, 1.0), 'H3': (1.5, 1.5)}

        elements = plotter.create_elements(layout=layout)

        # Find a node and check position
        node = [e for e in elements if e.get('data', {}).get('id') == 'H1'][0]
        assert 'position' in node
        # Positions should be scaled by 100
        assert node['position']['x'] == 50.0
        assert node['position']['y'] == 50.0

    def test_create_elements_node_sizes(self, simple_network):
        """Test that node sizes scale with frequency."""
        plotter = InteractiveCytoscapePlotter(simple_network)
        layout = {'H1': (0, 0), 'H2': (1, 0), 'H3': (2, 0)}

        elements = plotter.create_elements(layout=layout)

        # Get nodes
        nodes = {
            e['data']['id']: e for e in elements if 'source' not in e.get('data', {})
        }

        # H3 has frequency 3, should be larger than H2 (frequency 1)
        assert nodes['H3']['data']['size'] > nodes['H2']['data']['size']

    def test_create_elements_with_labels(self, simple_network):
        """Test label display."""
        plotter = InteractiveCytoscapePlotter(simple_network)
        layout = {'H1': (0, 0), 'H2': (1, 0), 'H3': (2, 0)}

        elements_with = plotter.create_elements(layout=layout, show_labels=True)
        elements_without = plotter.create_elements(layout=layout, show_labels=False)

        # With labels
        node_with = [e for e in elements_with if e['data']['id'] == 'H1'][0]
        assert node_with['data']['label'] == 'H1'

        # Without labels
        node_without = [e for e in elements_without if e['data']['id'] == 'H1'][0]
        assert node_without['data']['label'] == ''

    def test_create_elements_edge_labels(self, simple_network):
        """Test edge label display."""
        plotter = InteractiveCytoscapePlotter(simple_network)
        layout = {'H1': (0, 0), 'H2': (1, 0), 'H3': (2, 0)}

        elements = plotter.create_elements(layout=layout, show_edge_labels=True)

        # Find edge
        edge = [e for e in elements if e.get('data', {}).get('source') == 'H1'][0]
        assert edge['data']['label'] == '1'  # Distance is 1

    def test_create_elements_with_populations(self, network_with_populations):
        """Test elements with population data."""
        plotter = InteractiveCytoscapePlotter(network_with_populations)
        layout = {'H1': (0, 0), 'H2': (1, 0), 'H3': (2, 0)}
        pop_colors = {'PopA': '#FF0000', 'PopB': '#0000FF'}

        elements = plotter.create_elements(layout=layout, population_colors=pop_colors)

        # Find H2 which has both populations
        node = [e for e in elements if e['data']['id'] == 'H2'][0]

        assert node['data']['has_pie'] is True
        assert 'pie_data' in node['data']
        assert len(node['data']['pie_data']) == 2  # Two populations

    def test_create_stylesheet(self, simple_network):
        """Test stylesheet creation."""
        plotter = InteractiveCytoscapePlotter(simple_network)
        stylesheet = plotter.create_stylesheet()

        assert len(stylesheet) > 0
        # Should have rules for nodes, edges, etc.
        selectors = [rule['selector'] for rule in stylesheet]
        assert 'node' in selectors
        assert 'edge' in selectors

    def test_generate_population_colors(self, simple_network):
        """Test color generation for populations."""
        plotter = InteractiveCytoscapePlotter(simple_network)
        populations = ['PopA', 'PopB', 'PopC']

        colors = plotter.generate_population_colors(populations)

        assert len(colors) == 3
        assert all(pop in colors for pop in populations)
        # Colors should be hex strings
        assert all(color.startswith('#') for color in colors.values())
        # Colors should be unique
        assert len(set(colors.values())) == 3

    def test_median_vector_styling(self, simple_network):
        """Test median vector nodes are styled differently."""
        # Add a median vector
        seq_mv = Sequence('MV1', 'ATCG')
        mv = Haplotype(sequence=seq_mv, sample_ids=[])
        simple_network.add_haplotype(mv, median_vector=True)
        simple_network.add_edge('H1', 'MV1', distance=1)

        plotter = InteractiveCytoscapePlotter(simple_network)
        layout = {'H1': (0, 0), 'H2': (1, 0), 'H3': (2, 0), 'MV1': (0.5, 0.5)}

        elements = plotter.create_elements(layout=layout)

        # Find median vector node
        mv_node = [e for e in elements if e['data']['id'] == 'MV1'][0]
        assert mv_node['data']['is_median'] is True


class TestCreateCytoscapeNetwork:
    """Tests for create_cytoscape_network convenience function."""

    def test_basic_creation(self, simple_network):
        """Test basic network creation."""
        layout = {'H1': (0, 0), 'H2': (1, 0), 'H3': (2, 0)}

        elements, stylesheet = create_cytoscape_network(simple_network, layout=layout)

        assert len(elements) == 5  # 3 nodes + 2 edges
        assert len(stylesheet) > 0

    def test_auto_generate_colors(self, network_with_populations):
        """Test automatic color generation for populations."""
        layout = {'H1': (0, 0), 'H2': (1, 0), 'H3': (2, 0)}

        elements, stylesheet = create_cytoscape_network(
            network_with_populations, layout=layout
        )

        # Should have generated colors for populations
        # Nodes with population data should have pie data
        nodes_with_pie = [
            e
            for e in elements
            if e.get('data', {}).get('has_pie') and 'source' not in e.get('data', {})
        ]
        assert len(nodes_with_pie) > 0

    def test_custom_colors(self, network_with_populations):
        """Test with custom population colors."""
        layout = {'H1': (0, 0), 'H2': (1, 0), 'H3': (2, 0)}
        custom_colors = {'PopA': '#AAAAAA', 'PopB': '#BBBBBB'}

        elements, stylesheet = create_cytoscape_network(
            network_with_populations, layout=layout, population_colors=custom_colors
        )

        # Check that custom colors are used
        node = [e for e in elements if e['data']['id'] == 'H1'][0]
        # H1 is all PopA
        assert node['data']['color'] == '#AAAAAA'

    def test_without_layout(self, simple_network):
        """Test network creation without pre-computed layout."""
        # Should still work, using default positions
        elements, stylesheet = create_cytoscape_network(simple_network)

        assert len(elements) == 5
        # All nodes should have positions (0, 0) by default
        nodes = [e for e in elements if 'position' in e]
        assert len(nodes) == 3


class TestEmptyAndEdgeCases:
    """Tests for edge cases and empty networks."""

    def test_empty_network(self):
        """Test with empty network."""
        network = HaplotypeNetwork(name='Empty')
        elements, stylesheet = create_cytoscape_network(network)

        assert len(elements) == 0
        assert len(stylesheet) > 0  # Should still have default styles

    def test_single_node(self):
        """Test with single node network."""
        network = HaplotypeNetwork(name='Single')
        seq1 = Sequence('H1', 'ATCG')
        h1 = Haplotype(sequence=seq1, sample_ids=['S1'])
        network.add_haplotype(h1)

        layout = {'H1': (0, 0)}
        elements, stylesheet = create_cytoscape_network(network, layout=layout)

        assert len(elements) == 1  # Just one node
        assert elements[0]['data']['id'] == 'H1'

    def test_node_without_frequency(self):
        """Test node with zero frequency."""
        network = HaplotypeNetwork(name='Test')
        seq1 = Sequence('H1', 'ATCG')
        h1 = Haplotype(sequence=seq1, sample_ids=[])  # No samples
        network.add_haplotype(h1)

        layout = {'H1': (0, 0)}
        elements = InteractiveCytoscapePlotter(network).create_elements(layout=layout)

        # Should still create node with minimum size (even if frequency is 0)
        node = elements[0]
        assert node['data']['size'] > 0  # Should use minimum size


class TestEdgeTickMarks:
    """Tests for the mutation tick marks drawn along edges."""

    @pytest.mark.parametrize(
        'distance,expected',
        [
            (0, ''),
            (-3, ''),
            (1, '|'),
            (2, '| |'),
            (3, '| | |'),
            (MAX_TICK_MARKS + 1, ''),
        ],
    )
    def test_format_edge_ticks(self, distance, expected):
        """Ticks are space-joined pipes, empty outside the valid range."""
        assert format_edge_ticks(distance) == expected

    def test_format_edge_ticks_respects_custom_cap(self):
        """A lower cap sends the caller to the numeral fallback sooner."""
        assert format_edge_ticks(5, max_ticks=4) == ''
        assert format_edge_ticks(4, max_ticks=4) == '| | | |'

    def test_create_elements_emits_ticks(self, simple_network):
        """Edges carry a ticks string alongside the numeral label."""
        plotter = InteractiveCytoscapePlotter(simple_network)
        elements = plotter.create_elements()

        edges = {
            el['data']['id']: el['data']
            for el in elements
            if 'source' in el.get('data', {})
        }
        assert edges['H1-H2']['ticks'] == '|'
        assert edges['H2-H3']['ticks'] == '| |'
        # Regression guard: the numeral label must survive untouched.
        assert edges['H1-H2']['label'] == '1'
        assert edges['H2-H3']['label'] == '2'

    def test_create_elements_fit_tick_font(self, simple_network):
        """Every edge carries a tick font sized to its visible length."""
        plotter = InteractiveCytoscapePlotter(simple_network)
        wide = {'H1': (0.0, 0.0), 'H2': (5.0, 0.0), 'H3': (10.0, 0.0)}
        elements = plotter.create_elements(layout=wide)
        edges = [el['data'] for el in elements if 'source' in el.get('data', {})]

        assert all(e['tick_font'] == TICK_FONT_SIZE for e in edges)

        # Nodes at the same spot leave no visible edge: numeral fallback.
        stacked = dict.fromkeys(wide, (0.0, 0.0))
        elements = plotter.create_elements(layout=stacked)
        edges = [el['data'] for el in elements if 'source' in el.get('data', {})]

        assert all(e['tick_font'] == 0 for e in edges)

    def test_create_elements_no_ticks_when_labels_hidden(self, simple_network):
        """Hiding edge labels hides the ticks too."""
        plotter = InteractiveCytoscapePlotter(simple_network)
        elements = plotter.create_elements(show_edge_labels=False)

        for el in elements:
            if 'source' in el.get('data', {}):
                assert el['data']['ticks'] == ''

    def test_tick_stylesheet_split_at_threshold(self):
        """Short edges get ticks, long edges get the numeral."""
        rules = create_edge_tick_stylesheet(threshold=7)
        selectors = [rule['selector'] for rule in rules]

        assert selectors == [
            'edge[distance <= 7][tick_font > 0]',
            'edge[distance <= 7][tick_font = 0]',
            'edge[distance > 7]',
        ]
        assert rules[0]['style']['label'] == 'data(ticks)'
        assert rules[0]['style']['text-rotation'] == 'autorotate'
        # The font shrinks per edge so the comb fits between the nodes.
        assert rules[0]['style']['font-size'] == 'data(tick_font)'
        # No white pill behind the ticks (overrides the edge[label] rule).
        assert rules[0]['style']['text-background-opacity'] == 0
        # An edge too short for even the smallest comb shows the numeral.
        assert rules[1]['style']['label'] == 'data(label)'
        assert rules[2]['style']['label'] == 'data(label)'

    def test_tick_stylesheet_disabled(self):
        """With ticks off every edge falls back to a single numeral rule."""
        rules = create_edge_tick_stylesheet(show_ticks=False)

        assert len(rules) == 1
        assert rules[0]['selector'] == 'edge[distance > 0]'
        assert rules[0]['style']['label'] == 'data(label)'

    def test_tick_selectors_share_a_prefix(self):
        """The GUI strips rules by prefix, so every selector must match it."""
        for show_ticks in (True, False):
            rules = create_edge_tick_stylesheet(show_ticks=show_ticks)
            assert all(r['selector'].startswith('edge[distance') for r in rules)

    def test_tick_rules_come_after_edge_label(self, simple_network):
        """Cytoscape resolves conflicts by order, so ticks must win."""
        _, stylesheet = create_cytoscape_network(simple_network)
        selectors = [str(rule.get('selector', '')) for rule in stylesheet]

        edge_label = selectors.index('edge[label]')
        last_tick = max(
            i for i, sel in enumerate(selectors) if sel.startswith('edge[distance')
        )
        assert last_tick > edge_label

    def test_tick_threshold_passed_through(self, simple_network):
        """The threshold kwarg reaches the generated selectors."""
        _, stylesheet = create_cytoscape_network(simple_network, edge_tick_threshold=3)
        selectors = [str(rule.get('selector', '')) for rule in stylesheet]

        assert 'edge[distance <= 3][tick_font > 0]' in selectors
        assert not any(f'distance <= {DEFAULT_TICK_THRESHOLD}]' in s for s in selectors)

    def test_no_tick_rules_when_edge_labels_hidden(self, simple_network):
        """Hiding edge labels drops the mutation-count rules entirely."""
        _, stylesheet = create_cytoscape_network(simple_network, show_edge_labels=False)
        selectors = [str(rule.get('selector', '')) for rule in stylesheet]

        assert not any(sel.startswith('edge[distance') for sel in selectors)


class TestFitTickFontSize:
    """The tick comb shrinks with the visible edge and has a floor."""

    def test_long_edge_uses_the_full_size(self):
        """Plenty of room: no shrinking."""
        assert fit_tick_font_size(400.0, 10.0, 10.0, 5) == TICK_FONT_SIZE

    def test_shrinks_as_the_edge_shortens(self):
        """Closer nodes, smaller font, never below the floor."""
        sizes = [fit_tick_font_size(length, 10.0, 10.0, 5) for length in (120, 100, 80)]

        assert sizes[0] > sizes[1] > sizes[2] >= MIN_TICK_FONT_SIZE
        assert sizes[0] < TICK_FONT_SIZE

    def test_more_ticks_need_more_room(self):
        """A longer comb on the same edge gets a smaller font."""
        assert fit_tick_font_size(120.0, 10.0, 10.0, 8) < fit_tick_font_size(
            120.0, 10.0, 10.0, 2
        )

    @pytest.mark.parametrize(
        ('length', 'ra', 'rb', 'count'),
        [
            (30.0, 10.0, 10.0, 5),  # visible 10px: too short for five ticks
            (20.0, 10.0, 10.0, 1),  # discs touch
            (10.0, 10.0, 10.0, 1),  # discs overlap
            (100.0, 5.0, 5.0, 0),  # nothing to draw
        ],
    )
    def test_unfit_returns_zero(self, length, ra, rb, count):
        """Zero means 'show the numeral'."""
        assert fit_tick_font_size(length, ra, rb, count) == 0.0

    def test_floor_is_honoured_exactly(self):
        """Just enough room for the minimum font: keep the ticks."""
        # glyphs = 3, so visible * 0.7 / (3 * 0.6) == MIN when visible == MIN * 1.8 / 0.7
        visible = MIN_TICK_FONT_SIZE * 1.8 / 0.7
        assert fit_tick_font_size(visible + 20.0, 10.0, 10.0, 2) == pytest.approx(
            MIN_TICK_FONT_SIZE, abs=0.01
        )


class TestMedianVectorSelector:
    """The median vector rule must not match ordinary haplotypes."""

    def test_uses_truthy_selector_syntax(self, simple_network):
        """'[field = true]' is invalid and matched every node."""
        _, stylesheet = create_cytoscape_network(simple_network)
        selectors = [str(rule.get('selector', '')) for rule in stylesheet]

        assert 'node[?is_median]' in selectors
        assert 'node[is_median = true]' not in selectors

    def test_haplotype_nodes_keep_their_population_colour(self, simple_network):
        """The greying bug: only the base rule may colour a haplotype."""
        _, stylesheet = create_cytoscape_network(simple_network)
        median_rules = [
            rule for rule in stylesheet if 'is_median' in str(rule.get('selector', ''))
        ]

        assert median_rules
        for rule in median_rules:
            assert rule['selector'].startswith('node[?')


class TestNodeTooltip:
    """The hover tooltip payload carried in each node's data."""

    def test_haplotype_tooltip(self, simple_network):
        """A haplotype node reports its label, count and populations."""
        plotter = InteractiveCytoscapePlotter(simple_network)
        nodes = {
            el['data']['id']: el['data']['tooltip']
            for el in plotter.create_elements()
            if 'source' not in el.get('data', {})
        }

        assert nodes['H1']['kind'] == 'haplotype'
        assert nodes['H1']['label'] == 'H1'
        assert nodes['H1']['frequency'] == 2
        assert nodes['H1']['populations'] == [['PopA', 2]]
        assert nodes['H1']['samples'] == ['S1', 'S2']
        assert nodes['H1']['extra'] == 0

    def test_tooltip_uses_custom_label(self, simple_network):
        """The tooltip title follows the H-number mapping, not the node id."""
        plotter = InteractiveCytoscapePlotter(simple_network)
        elements = plotter.create_elements(node_labels={'H1': 'Hap-One'})

        tooltips = {
            el['data']['id']: el['data']['tooltip']
            for el in elements
            if 'source' not in el.get('data', {})
        }
        assert tooltips['H1']['label'] == 'Hap-One'

    def test_tooltip_falls_back_to_node_id(self, simple_network):
        """With labels hidden the tooltip still names the node."""
        plotter = InteractiveCytoscapePlotter(simple_network)
        elements = plotter.create_elements(show_labels=False)

        tooltips = [
            el['data']['tooltip']
            for el in elements
            if 'source' not in el.get('data', {})
        ]
        assert all(t['label'] for t in tooltips)

    def test_median_vector_tooltip(self):
        """A median vector says so instead of listing samples."""
        tooltip = build_node_tooltip('H9', None, True)

        assert tooltip['kind'] == 'median'
        assert tooltip['samples'] == []
        assert tooltip['frequency'] == 0

    def test_sample_list_is_truncated(self):
        """Long sample lists are capped and the remainder is counted."""

        class FakeHap:
            sample_ids = [f'S{i}' for i in range(25)]
            frequency = 25

            @staticmethod
            def get_frequency_by_population():
                return {'PopA': 25}

        tooltip = build_node_tooltip('H1', FakeHap(), False)

        assert len(tooltip['samples']) == MAX_TOOLTIP_SAMPLES
        assert tooltip['extra'] == 25 - MAX_TOOLTIP_SAMPLES

    def test_tooltip_is_json_safe(self, simple_network):
        """The payload rides to the browser inside the elements prop."""
        plotter = InteractiveCytoscapePlotter(simple_network)
        elements = plotter.create_elements()

        json.dumps(elements)

    def test_hover_string_is_gone(self, simple_network):
        """The old unused hover string must not come back."""
        plotter = InteractiveCytoscapePlotter(simple_network)

        for el in plotter.create_elements():
            assert 'hover' not in el.get('data', {})


class TestResolvePopulationCounts:
    """Shared population counting for pie charts and tooltips."""

    def test_prefers_haplotype_counts(self):
        """A haplotype that knows its populations is trusted."""

        class FakeHap:
            sample_ids = ['S1', 'S2']

            @staticmethod
            def get_frequency_by_population():
                return {'PopA': 2}

        counts = resolve_population_counts(FakeHap(), {'S1': 'PopZ', 'S2': 'PopZ'})

        assert counts == {'PopA': 2}

    def test_falls_back_to_mapping(self):
        """Unassigned-only counts mean the metadata arrived later."""

        class FakeHap:
            sample_ids = ['S1', 'S2']

            @staticmethod
            def get_frequency_by_population():
                return {'Unassigned': 2}

        counts = resolve_population_counts(FakeHap(), {'S1': 'PopA', 'S2': 'PopB'})

        assert counts == {'PopA': 1, 'PopB': 1}

    def test_samples_missing_from_mapping_are_unassigned(self):
        """A sample with no metadata row still gets counted."""

        class FakeHap:
            sample_ids = ['S1', 'S2']

            @staticmethod
            def get_frequency_by_population():
                return {}

        counts = resolve_population_counts(FakeHap(), {'S1': 'PopA'})

        assert counts == {'PopA': 1, 'Unassigned': 1}
