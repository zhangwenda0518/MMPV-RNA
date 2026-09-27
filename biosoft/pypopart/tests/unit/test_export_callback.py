"""Tests for the GUI export callback.

The underlying writers are covered in test_network_export.py; this
exercises the Dash wiring around them, which had no coverage at all and
where a missing branch silently returned None.
"""

import json
import logging

import dash
from dash.exceptions import PreventUpdate
import pytest

from pypopart.core.graph import HaplotypeNetwork
from pypopart.core.haplotype import Haplotype
from pypopart.core.sequence import Sequence
from pypopart.gui.callbacks import export as export_callbacks
from pypopart.gui.serialization import network_to_store


class CallbackCollector:
    """Stand-in for a Dash app that just records callback functions."""

    def __init__(self):
        self.callbacks = {}

    def callback(self, *args, **kwargs):
        """
        Register a callback by function name.

        Parameters
        ----------
        *args : object
            Ignored Dash dependency arguments.
        **kwargs : object
            Ignored Dash options.

        Returns
        -------
        callable
            Decorator storing the function under its name.
        """

        def decorator(func):
            self.callbacks[func.__name__] = func
            return func

        return decorator


@pytest.fixture
def export_network_callback():
    """Register the export callbacks and hand back the main one."""
    app = CallbackCollector()
    export_callbacks.register(app, logging.getLogger('test'))
    return app.callbacks['export_network']


@pytest.fixture
def network_store():
    """Build a serialized two-haplotype network, as the store holds it."""
    network = HaplotypeNetwork(name='ExportTest')
    network.add_haplotype(Haplotype(Sequence('H1', 'ATCG'), sample_ids=['S1']))
    network.add_haplotype(Haplotype(Sequence('H2', 'ATGG'), sample_ids=['S2']))
    network.add_edge('H1', 'H2', distance=1)
    return network_to_store(network)


class TestTextExports:
    """GraphML, GML and JSON are written server-side."""

    @pytest.mark.parametrize(
        'fmt,filename,mimetype',
        [
            ('graphml', 'network.graphml', 'text/xml'),
            ('gml', 'network.gml', 'text/plain'),
            ('json', 'network.json', 'application/json'),
        ],
    )
    def test_payload_shape(
        self, export_network_callback, network_store, fmt, filename, mimetype
    ):
        """dcc.Download needs content, filename and type."""
        payload, image, _, _, is_open = export_network_callback(1, network_store, fmt)

        assert payload['filename'] == filename
        assert payload['type'] == mimetype
        assert payload['content']
        assert image is dash.no_update
        # The user gets told it worked.
        assert is_open is True

    def test_json_content_is_parseable(self, export_network_callback, network_store):
        """The bytes handed to the browser must be a usable file."""
        payload, _, _, _, _ = export_network_callback(1, network_store, 'json')
        parsed = json.loads(payload['content'])

        assert {'nodes', 'edges', 'metadata'} <= set(parsed)
        assert len(parsed['nodes']) == 2

    def test_graphml_content_is_xml(self, export_network_callback, network_store):
        """A GraphML file Cytoscape or Gephi could actually open."""
        payload, _, _, _, _ = export_network_callback(1, network_store, 'graphml')

        assert payload['content'].lstrip().startswith('<?xml')
        assert 'graphml' in payload['content']


class TestImageExports:
    """PNG comes off the canvas; SVG is rendered server-side."""

    def test_png_asks_the_canvas(self, export_network_callback, network_store):
        """The component renders exactly what the user is looking at."""
        payload, image, _, _, is_open = export_network_callback(1, network_store, 'png')

        assert payload is dash.no_update
        assert image['type'] == 'png'
        assert image['action'] == 'download'
        # dash-cytoscape appends the extension, so passing one gave
        # the user a file called network.png.png.
        assert image['filename'] == 'network'
        assert is_open is True

    def test_svg_is_rendered_server_side(self, export_network_callback, network_store):
        """Render SVG server-side.

        Cytoscape.js has no SVG output without an extension dash-cytoscape
        does not bundle, so asking the canvas produced nothing at all.
        """
        payload, image, _, _, is_open = export_network_callback(1, network_store, 'svg')

        assert image is dash.no_update
        assert payload['filename'] == 'network.svg'
        assert payload['type'] == 'image/svg+xml'
        assert '<svg' in payload['content']
        assert is_open is True

    def test_svg_uses_the_current_positions(
        self, export_network_callback, network_store
    ):
        """The exported figure should match the layout on screen."""
        payload, _, _, _, _ = export_network_callback(
            1,
            network_store,
            'svg',
            {'H1': [0.0, 0.0], 'H2': [1.0, 1.0]},
            {'H2': [5.0, 5.0]},
            None,
        )

        assert '<svg' in payload['content']


class TestExportGuards:
    """The failure paths that used to be silent."""

    def test_no_network_prevents_update(self, export_network_callback):
        """Nothing computed yet, nothing to export."""
        with pytest.raises(PreventUpdate):
            export_network_callback(1, None, 'json')

    def test_unknown_format_reports_instead_of_returning_none(
        self, export_network_callback, network_store
    ):
        """Report the problem instead of falling off the end returning None."""
        payload, image, children, header, is_open = export_network_callback(
            1, network_store, 'not-a-format'
        )

        assert payload is dash.no_update
        assert image is dash.no_update
        assert is_open is True
        assert header == 'Export failed'

    def test_broken_network_data_reports(self, export_network_callback):
        """A failure surfaces as a toast, not as a dead button."""
        _, _, _, header, is_open = export_network_callback(
            1, {'nodes': 'not-a-list'}, 'json'
        )

        assert header == 'Export failed'
        assert is_open is True

    def test_every_dropdown_option_is_handled(self, export_network_callback):
        """Guard against the UI offering a format the callback lacks."""
        from pypopart.gui.callbacks.export import IMAGE_FORMATS, TEXT_EXPORTERS
        from pypopart.gui.layout.cards import create_export_card

        def find_options(component):
            if getattr(component, 'id', None) == 'export-format':
                return component.options
            children = getattr(component, 'children', None)
            if isinstance(children, (list, tuple)):
                for child in children:
                    hit = find_options(child)
                    if hit is not None:
                        return hit
            elif children is not None:
                return find_options(children)
            return None

        offered = {opt['value'] for opt in find_options(create_export_card())}
        handled = set(TEXT_EXPORTERS) | set(IMAGE_FORMATS)

        assert offered == handled


@pytest.fixture
def coloured_store():
    """Build a network plus metadata whose populations drive the colours."""
    network = HaplotypeNetwork(name='ColourTest')
    network.add_haplotype(Haplotype(Sequence('A', 'ATCG'), sample_ids=['s1', 's2']))
    network.add_haplotype(Haplotype(Sequence('B', 'ATGG'), sample_ids=['s3']))
    network.add_edge('A', 'B', distance=4)
    metadata = {
        'populations': {'s1': 'Alpha', 's2': 'Alpha', 's3': 'Beta'},
        'population_colors': {'Alpha': '#123456', 'Beta': '#abcdef'},
    }
    return network_to_store(network), metadata


class TestFigureFidelity:
    """The exported figure must match what the canvas shows."""

    def test_h_numbers_are_supplied(self, export_network_callback, network_store):
        """Raw node IDs are meaningless to a reader of the figure.

        Asserted on the labels handed to the plotter rather than on the
        SVG text, because matplotlib writes glyph outlines, not
        characters -- checking the output for 'H1' passes by accident.
        """
        from pypopart.gui.callbacks.export import _default_h_numbers

        # Real node IDs are sequence names, not H numbers.
        network = HaplotypeNetwork()
        network.add_haplotype(Haplotype(Sequence('Zeta_02', 'ATCG'), sample_ids=['a']))
        network.add_haplotype(Haplotype(Sequence('Alpha_01', 'ATGG'), sample_ids=['b']))

        assert _default_h_numbers(network) == {'Alpha_01': 'H1', 'Zeta_02': 'H2'}

    def test_population_colours_reach_the_figure(
        self, export_network_callback, coloured_store
    ):
        """Colours come from the metadata store, not the haplotypes.

        A network rebuilt from the store usually has no population data of
        its own, so resolving from the haplotypes alone left every node
        the default colour.
        """
        store, metadata = coloured_store
        payload, _, _, _, _ = export_network_callback(
            1, store, 'svg', None, None, metadata
        )

        content = payload['content'].lower()
        assert '#123456' in content or '123456' in content

    def test_legend_is_optional(self, export_network_callback, coloured_store):
        """The legend switch controls whether a key is drawn."""
        store, metadata = coloured_store
        with_legend, _, _, _, _ = export_network_callback(
            1, store, 'svg', None, None, metadata, None, True, 10, True
        )
        without, _, _, _, _ = export_network_callback(
            1, store, 'svg', None, None, metadata, None, True, 10, False
        )

        assert len(with_legend['content']) > len(without['content'])

    def test_legend_names_the_populations(
        self, export_network_callback, coloured_store
    ):
        """A key is only useful if it carries the population labels."""
        store, metadata = coloured_store
        payload, _, _, _, _ = export_network_callback(
            1, store, 'svg', None, None, metadata, None, True, 10, True
        )

        assert 'Alpha' in payload['content']
        assert 'Beta' in payload['content']


class TestResolveNodeColors:
    """Mapping nodes to population colours for the figure."""

    def test_dominant_population_wins(self):
        """Matplotlib cannot draw a pie, so a mixed node takes the mode."""
        from pypopart.gui.callbacks.export import resolve_node_colors

        network = HaplotypeNetwork()
        network.add_haplotype(
            Haplotype(Sequence('A', 'ATCG'), sample_ids=['s1', 's2', 's3'])
        )
        colors = resolve_node_colors(
            network,
            {'Alpha': '#111111', 'Beta': '#222222'},
            {'s1': 'Alpha', 's2': 'Alpha', 's3': 'Beta'},
        )

        assert colors['A'] == '#111111'

    def test_unassigned_is_not_a_population(self):
        """Samples with no metadata must not out-vote a real population."""
        from pypopart.gui.callbacks.export import resolve_node_colors

        network = HaplotypeNetwork()
        network.add_haplotype(
            Haplotype(Sequence('A', 'ATCG'), sample_ids=['s1', 's2', 's3'])
        )
        colors = resolve_node_colors(network, {'Alpha': '#111111'}, {'s1': 'Alpha'})

        assert colors['A'] == '#111111'

    def test_no_colours_yields_nothing(self):
        """Without a palette there is nothing to resolve."""
        from pypopart.gui.callbacks.export import resolve_node_colors

        network = HaplotypeNetwork()
        network.add_haplotype(Haplotype(Sequence('A', 'ATCG'), sample_ids=['s1']))

        assert resolve_node_colors(network, None, {'s1': 'Alpha'}) == {}
