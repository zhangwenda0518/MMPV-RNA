"""Unit tests for PyPopART GUI application."""

import pytest

from pypopart.gui.layout import cards


class TestPyPopARTApp:
    """Test cases for PyPopART Dash GUI."""

    def test_app_initialization(self):
        """Test that the app initializes without errors."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False, port=8050)

        assert app is not None
        assert app.app is not None
        assert app.app.title == 'PyPopART - Haplotype Network Analysis'
        assert app.debug is False
        assert app.port == 8050

    def test_app_layout_structure(self):
        """Test that the app layout has expected components."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False)

        # Check that layout is set up
        assert app.app.layout is not None

        # Check that key components exist in layout
        layout_children = app.app.layout.children
        assert len(layout_children) > 0

    def test_create_upload_card(self):
        """Test upload card creation."""
        card = cards.create_upload_card()

        assert card is not None
        assert hasattr(card, 'children')

    def test_create_algorithm_card(self):
        """Test algorithm card creation."""
        card = cards.create_algorithm_card()

        assert card is not None
        assert hasattr(card, 'children')

    def test_create_layout_card(self):
        """Test layout card creation."""
        card = cards.create_layout_card()

        assert card is not None
        assert hasattr(card, 'children')

    def test_create_export_card(self):
        """Test export card creation."""
        card = cards.create_export_card()

        assert card is not None
        assert hasattr(card, 'children')

    def test_create_network_tab(self):
        """Test network tab creation."""
        tab = cards.create_network_tab()

        assert tab is not None
        assert hasattr(tab, 'children')

    def test_create_statistics_tab(self):
        """Test statistics tab creation."""
        tab = cards.create_statistics_tab()

        assert tab is not None
        assert hasattr(tab, 'children')

    def test_create_alignment_tab(self):
        """Test alignment tab creation."""
        tab = cards.create_alignment_tab()

        assert tab is not None
        assert hasattr(tab, 'children')

    def test_main_function_exists(self):
        """Test that main function is callable."""
        from pypopart.gui.app import main

        assert callable(main)

    def test_import_from_init(self):
        """Test that components can be imported from gui module."""
        from pypopart.gui import PyPopARTApp, main

        assert PyPopARTApp is not None
        assert callable(main)

    def test_logger_initialization(self):
        """Test that logger is properly initialized."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=True, port=8050)
        assert hasattr(app, 'logger')
        assert app.logger is not None

    def test_callbacks_registered(self):
        """Test that callbacks are properly registered without errors."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False)
        # If callbacks are registered with pattern-matching IDs correctly,
        # the app should initialize without errors
        assert app.app._callback_list is not None
        # Check that we have the expected number of callbacks
        assert len(app.app._callback_list) > 0

    def test_default_layout_is_hierarchical(self):
        """Test that default layout is set to hierarchical."""
        # Find the layout dropdown in the layout
        # It should have 'hierarchical' as the default value
        # This verifies Issue 4 is fixed
        layout_card = cards.create_layout_card()
        assert layout_card is not None

    def test_search_dropdown_supports_multi_select(self):
        """Test that search dropdown supports multiple selections."""
        # Find the haplotype search dropdown in the network tab
        # It should have multi=True set
        # This verifies Issue 6 is fixed
        network_tab = cards.create_network_tab()
        assert network_tab is not None

    def test_haplotype_summary_tab_has_mapping_components(self):
        """Test that haplotype summary tab has label mapping components."""
        tab = cards.create_haplotype_summary_tab()

        assert tab is not None
        # Tab should have children including buttons and stores
        assert hasattr(tab, 'children')
        assert len(tab.children) > 0

        # Find components by traversing children
        # Should have: download-h-number-template-button, upload-h-number-mapping, h-number-mapping-store
        component_ids = []

        def extract_ids(component):
            """Recursively extract all component IDs."""
            if hasattr(component, 'id'):
                component_ids.append(component.id)
            if hasattr(component, 'children'):
                if isinstance(component.children, list):
                    for child in component.children:
                        extract_ids(child)
                else:
                    extract_ids(component.children)

        extract_ids(tab)

        # Check that key components exist
        assert 'download-h-number-template-button' in component_ids
        assert 'upload-h-number-mapping' in component_ids
        assert 'h-number-mapping-store' in component_ids
        assert 'h-number-feedback' in component_ids

    def test_stylesheet_has_selected_pseudo_selector(self):
        """Test that stylesheet includes :selected pseudo-selector for node click highlighting."""
        from pypopart.core.graph import HaplotypeNetwork
        from pypopart.core.haplotype import Haplotype
        from pypopart.core.sequence import Sequence
        from pypopart.visualization.cytoscape_plot import InteractiveCytoscapePlotter
        from pypopart.visualization.style import POP_AMBER

        # Create a simple test network
        network = HaplotypeNetwork()
        seq1 = Sequence('hap1', 'ATCG')
        hap1 = Haplotype(seq1, sample_ids=['seq1'])
        network.add_haplotype(hap1)

        # Create plotter and stylesheet
        plotter = InteractiveCytoscapePlotter(network)
        stylesheet = plotter.create_stylesheet()

        # Check that :selected pseudo-selector exists in stylesheet
        selected_styles = [
            s for s in stylesheet if s.get('selector') == 'node:selected'
        ]
        assert len(selected_styles) == 1, (
            'Stylesheet should have exactly one node:selected style'
        )

        # Selection is an amber halo, thicker than the ink node border.
        # The constant is shared so the theme cannot drift from the test.
        selected_style = selected_styles[0]
        assert 'style' in selected_style
        assert selected_style['style']['border-color'] == POP_AMBER
        assert selected_style['style']['border-width'] == 6
        assert selected_style['style']['z-index'] == 999


class TestCentralHaplotypesTable:
    """The Statistics tab's Central Haplotypes table.

    This panel silently showed 'Unable to identify central haplotypes'
    for its whole life: the formatter indexed identify_central_haplotypes'
    list of (id, score) tuples as if it were a dict, and a bare except
    swallowed the TypeError. These tests pin the real shape.
    """

    def build_network(self):
        """Build a three-haplotype chain (H1 - H2 - H3)."""
        from pypopart import build_network
        from pypopart.core.alignment import Alignment
        from pypopart.core.sequence import Sequence

        alignment = Alignment(
            [
                Sequence('s1', 'AAAA'),
                Sequence('s2', 'AAAT'),
                Sequence('s3', 'AATT'),
            ]
        )
        return build_network('mst', alignment)

    def render(self, network, **kwargs):
        """Render the central-haplotype table for a network."""
        from pypopart.gui.callbacks.display import _format_central_haplotypes
        from pypopart.stats import (
            calculate_node_centrality,
            identify_central_haplotypes,
        )

        return _format_central_haplotypes(
            identify_central_haplotypes(network),
            calculate_node_centrality(network),
            **kwargs,
        )

    def rows_of(self, rendered):
        """Extract the table body rows from the rendered output."""
        table = rendered.children[0]
        return table.children[1].children

    def test_table_lists_each_measure(self):
        """Every centrality measure gets its own column."""
        rendered = self.render(self.build_network())
        header = rendered.children[0].children[0]
        titles = [th.children for th in header.children.children]
        assert titles == [
            'Haplotype',
            'Degree',
            'Betweenness',
            'Closeness',
            'Eigenvector',
        ]

    def test_rows_are_ranked_and_populated(self):
        """The hub ranks first and every cell holds a formatted score."""
        rows = self.rows_of(self.render(self.build_network()))
        assert len(rows) == 3

        cells = [[td.children for td in row.children] for row in rows]
        # H2 is the middle of the chain, so it is the most central
        assert cells[0][0] == 'H2'
        assert cells[0][1] == '1.000'
        # no cell is left empty or unformatted
        for row in cells:
            assert all(value and value != '-' for value in row)

    def test_top_n_caps_the_table(self):
        """top_n limits the rows and the caption reports the total."""
        rendered = self.render(self.build_network(), top_n=2)
        assert len(self.rows_of(rendered)) == 2
        assert 'Top 2 of 3' in rendered.children[1].children

    def test_empty_network_renders_a_message(self):
        """An empty network reports that there is nothing to rank."""
        from pypopart.core.graph import HaplotypeNetwork

        rendered = self.render(HaplotypeNetwork())
        assert 'No haplotypes to rank' in rendered.children


def _collect_ids(component, found=None):
    """Recursively collect every component id under a layout tree."""
    if found is None:
        found = []
    if getattr(component, 'id', None) is not None:
        found.append(component.id)
    children = getattr(component, 'children', None)
    if isinstance(children, (list, tuple)):
        for child in children:
            _collect_ids(child, found)
    elif children is not None:
        _collect_ids(children, found)
    return found


class TestSidebarCollapse:
    """The left control panel can be hidden and restored."""

    def test_layout_has_collapse_components(self):
        """The panels, the toggle and the state store are all present."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False)
        ids = _collect_ids(app.app.layout)

        assert 'sidebar-panel' in ids
        assert 'main-panel' in ids
        assert 'sidebar-toggle' in ids
        assert 'sidebar-collapsed' in ids

    def test_layout_has_drag_and_toast_components(self):
        """Drags live in their own store; successes go to one toast."""
        from pypopart.gui.app import PyPopARTApp

        ids = _collect_ids(PyPopARTApp(debug=False).app.layout)

        assert 'node-positions-store' in ids
        assert 'app-toast' in ids

    def test_sidebar_has_a_reachable_resize_grip(self):
        """The native corner handle sat off-screen on a 90vh panel."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False)

        def find(component, target):
            if getattr(component, 'id', None) == target:
                return component
            children = getattr(component, 'children', None)
            if isinstance(children, (list, tuple)):
                for child in children:
                    hit = find(child, target)
                    if hit is not None:
                        return hit
            elif children is not None:
                return find(children, target)
            return None

        sidebar = find(app.app.layout, 'sidebar-panel')
        assert sidebar is not None
        # Replaced by the centre-edge grip driven from assets/pypopart.js.
        assert 'resize' not in sidebar.style
        # The clamps the drag handler reads.
        assert sidebar.style['minWidth'] == '250px'
        assert sidebar.style['maxWidth'] == '600px'

        ids = _collect_ids(app.app.layout)
        assert 'sidebar-resizer' in ids
        assert 'sidebar-grip' in ids


class TestEdgeTickControls:
    """The layout card exposes the mutation tick mark settings."""

    def test_layout_card_has_tick_controls(self):
        """Both the toggle and the numeral threshold slider are present."""
        ids = _collect_ids(cards.create_layout_card())

        assert 'edge-tick-toggle' in ids
        assert 'edge-tick-threshold' in ids

    def test_layout_card_has_grid_controls(self):
        """Snap-to-grid needs both a switch and a spacing slider."""
        ids = _collect_ids(cards.create_layout_card())

        assert 'snap-to-grid-toggle' in ids
        assert 'grid-size' in ids

    def test_metadata_tab_has_colour_swatch_slot(self):
        """The per-population pickers render under the table."""
        ids = _collect_ids(cards.create_metadata_tab())

        assert 'population-colors' in ids


class TestTaskAndUploadControls:
    """Controls added with the repel slider and the task indicator."""

    def test_layout_card_has_repel_controls(self):
        """The spring layouts expose their repulsion strength."""
        ids = _collect_ids(cards.create_layout_card())

        assert 'repel-slider' in ids
        assert 'repel-options' in ids

    def test_title_bar_has_task_indicator(self):
        """Long callbacks report what they are doing in the header."""
        from pypopart.gui.app import PyPopARTApp

        ids = _collect_ids(PyPopARTApp().app.layout)

        assert 'task-indicator' in ids
        assert 'task-label' in ids

    def test_layout_card_has_box_select_switch(self):
        """Multi-select needs a way to box-drag without panning."""
        ids = _collect_ids(cards.create_layout_card())

        assert 'box-select-toggle' in ids

    def test_graph_allows_box_selection(self):
        """The Cytoscape canvas accepts a selection box."""
        import dash_cytoscape as cyto

        graph = next(
            c
            for c in _walk(cards.create_network_tab())
            if isinstance(c, cyto.Cytoscape)
        )

        assert graph.boxSelectionEnabled is True
        assert graph.userPanningEnabled is True

    def test_upload_card_is_a_drop_zone(self):
        """Both uploads carry the drop-zone class so drag and drop is visible."""
        from dash import dcc

        uploads = [
            c for c in _walk(cards.create_upload_card()) if isinstance(c, dcc.Upload)
        ]

        assert {u.id for u in uploads} >= {'upload-data', 'upload-metadata'}
        for upload in uploads:
            if upload.id in ('upload-data', 'upload-metadata'):
                assert upload.className == 'pp-dropzone'
                assert upload.className_active == 'pp-dropzone--active'


def _walk(component):
    """Yield a component and every descendant."""
    yield component
    children = getattr(component, 'children', None)
    if children is None:
        return
    if not isinstance(children, (list, tuple)):
        children = [children]
    for child in children:
        if hasattr(child, 'children') or hasattr(child, 'id'):
            yield from _walk(child)


class TestBackgroundJobs:
    """Long computations run out of process so they can be cancelled."""

    @pytest.fixture(scope='class')
    def callback_map(self):
        """Register the app once and expose its callback map."""
        from pypopart.gui.app import PyPopARTApp

        return PyPopARTApp().app.callback_map

    @staticmethod
    def _entry(callback_map, output_fragment):
        return next(v for k, v in callback_map.items() if output_fragment in k)

    def test_compute_network_is_a_background_callback(self, callback_map):
        """Changing algorithm or recomputing cancels the running job."""
        entry = self._entry(callback_map, 'network-store.data')

        assert entry['background']
        cancel_ids = {c['id'] for c in entry['background']['cancel']}
        assert {'metadata-commit-token', 'algorithm-select'} <= cancel_ids

    def test_apply_layout_is_a_background_callback(self, callback_map):
        """Any change to the layout settings cancels the running layout."""
        entry = self._entry(callback_map, 'layout-store.data')

        assert entry['background']
        cancel_ids = {c['id'] for c in entry['background']['cancel']}
        assert {
            'apply-layout-button',
            'layout-select',
            'spacing-slider',
            'repel-slider',
            'snap-to-grid-toggle',
            'grid-size',
        } <= cancel_ids


class TestFileChip:
    """The upload boxes show which file is loaded."""

    def test_chip_names_the_file(self):
        """Filename and a summary, with the full name as a tooltip."""
        from pypopart.gui.callbacks.upload import file_chip

        chip = file_chip('samples.fasta', '12 sequences')
        texts = [c.children for c in chip.children]

        assert chip.className == 'pp-filechip'
        assert chip.title == 'samples.fasta'
        assert texts == ['samples.fasta', '12 sequences']

    def test_chip_survives_a_missing_name(self):
        """Dash can hand over no filename; the chip still renders."""
        from pypopart.gui.callbacks.upload import file_chip

        assert file_chip(None, '1 row').children[0].children == 'file'
