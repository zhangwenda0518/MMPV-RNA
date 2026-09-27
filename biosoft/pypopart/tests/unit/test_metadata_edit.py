"""Tests for the editable metadata table helpers."""

import pytest

from pypopart.gui.metadata_edit import (
    apply_population_color,
    build_metadata_datatable,
    build_metadata_rows,
    build_population_color_controls,
    diff_metadata_rows,
    population_colors_from_rows,
    rows_to_metadata_store,
)


@pytest.fixture
def metadata_store():
    """Build a metadata store in the shape the CSV upload writes."""
    return {
        'raw': {
            'S1': {'population': 'PopA', 'color': '#ff0000'},
            'S2': {'population': 'PopB', 'color': '#00ff00'},
        },
        'coordinates': {'S1': {'lat': -35.3, 'lon': 149.1}},
        'populations': {'S1': 'PopA', 'S2': 'PopB'},
        'population_colors': {'PopA': '#ff0000', 'PopB': '#00ff00'},
        'sequence_ids': ['S1', 'S2'],
    }


class TestBuildMetadataRows:
    """Rows shown in the metadata table."""

    def test_rows_without_metadata(self):
        """Alignment-only sequences get blank metadata and a cross."""
        rows = build_metadata_rows(['S2', 'S1'], None)

        assert [row['id'] for row in rows] == ['S1', 'S2']
        assert all(row['in_alignment'] == '✓' for row in rows)
        assert all(row['in_metadata'] == '✗' for row in rows)
        assert all(row['population'] == '' for row in rows)

    def test_rows_merge_metadata(self, metadata_store):
        """Population, colour and coordinates are pulled through."""
        rows = {
            row['id']: row for row in build_metadata_rows(['S1', 'S2'], metadata_store)
        }

        assert rows['S1']['population'] == 'PopA'
        assert rows['S1']['color'] == '#ff0000'
        assert rows['S1']['latitude'] == '-35.3'
        assert rows['S1']['longitude'] == '149.1'
        # S2 has no coordinates.
        assert rows['S2']['latitude'] == ''

    def test_rows_include_metadata_only_ids(self, metadata_store):
        """An ID only in the CSV still gets a row, flagged as missing."""
        rows = {row['id']: row for row in build_metadata_rows(['S1'], metadata_store)}

        assert set(rows) == {'S1', 'S2'}
        assert rows['S2']['in_alignment'] == '✗'
        assert rows['S2']['in_metadata'] == '✓'

    def test_values_are_all_strings(self, metadata_store):
        """Rows go through a dcc.Store, so nothing may be a float or None."""
        for row in build_metadata_rows(['S1', 'S2'], metadata_store):
            assert all(isinstance(value, str) for value in row.values())


class TestDiffMetadataRows:
    """Detecting which cells the user changed."""

    def test_identical_rows_have_no_changes(self, metadata_store):
        """A freshly rendered table is not a draft."""
        rows = build_metadata_rows(['S1', 'S2'], metadata_store)

        assert diff_metadata_rows(rows, rows) == []

    def test_shuffled_rows_have_no_changes(self, metadata_store):
        """Native sorting reorders rows; that is not an edit."""
        rows = build_metadata_rows(['S1', 'S2'], metadata_store)

        assert diff_metadata_rows(rows, list(reversed(rows))) == []

    def test_edited_cell_is_reported(self, metadata_store):
        """Changing a population is reported against its sequence ID."""
        rows = build_metadata_rows(['S1', 'S2'], metadata_store)
        edited = [dict(row) for row in rows]
        edited[0]['population'] = 'PopC'

        assert diff_metadata_rows(rows, edited) == [('S1', 'population')]

    def test_non_editable_columns_are_ignored(self, metadata_store):
        """Only the editable columns can produce a draft."""
        rows = build_metadata_rows(['S1', 'S2'], metadata_store)
        edited = [dict(row) for row in rows]
        edited[0]['in_alignment'] = '✗'

        assert diff_metadata_rows(rows, edited) == []


class TestRowsToMetadataStore:
    """Turning edited rows back into the metadata store shape."""

    def test_round_trip_is_stable(self, metadata_store):
        """Committing an unedited table reproduces the same store."""
        rows = build_metadata_rows(['S1', 'S2'], metadata_store)
        store, warnings = rows_to_metadata_store(rows, metadata_store)

        assert warnings == []
        assert store['populations'] == metadata_store['populations']
        assert store['population_colors'] == metadata_store['population_colors']
        assert store['coordinates'] == metadata_store['coordinates']
        assert set(store['sequence_ids']) == {'S1', 'S2'}

    def test_edited_population_reaches_the_store(self, metadata_store):
        """The whole point: a typed population is applied."""
        rows = build_metadata_rows(['S1', 'S2'], metadata_store)
        rows[0]['population'] = 'PopC'
        store, _ = rows_to_metadata_store(rows, metadata_store)

        assert store['populations']['S1'] == 'PopC'

    def test_colour_first_wins_per_population(self):
        """Matching the CSV path: the first colour seen for a population."""
        rows = [
            {'id': 'S1', 'population': 'PopA', 'color': '#111111'},
            {'id': 'S2', 'population': 'PopA', 'color': '#222222'},
        ]
        store, _ = rows_to_metadata_store(rows)

        assert store['population_colors']['PopA'] == '#111111'

    def test_uncoloured_population_is_filled_in(self):
        """A population typed without a colour still gets one."""
        rows = [{'id': 'S1', 'population': 'PopA', 'color': ''}]
        store, _ = rows_to_metadata_store(rows)

        assert store['population_colors']['PopA'].startswith('#')

    def test_no_populations_gives_none_colors(self):
        """Downstream code checks for None, not an empty dict."""
        rows = [{'id': 'S1', 'population': '', 'color': ''}]
        store, _ = rows_to_metadata_store(rows)

        assert store['population_colors'] is None

    def test_coordinates_parse(self):
        """Numeric strings become floats in the store."""
        rows = [
            {'id': 'S1', 'population': 'A', 'latitude': '-35.3', 'longitude': '149.1'}
        ]
        store, warnings = rows_to_metadata_store(rows)

        assert warnings == []
        assert store['coordinates']['S1'] == {'lat': -35.3, 'lon': 149.1}

    @pytest.mark.parametrize(
        'lat,lon',
        [('abc', '149.1'), ('-35.3', 'xyz'), ('120.0', '149.1'), ('-35.3', '999')],
    )
    def test_bad_coordinates_warn_and_drop(self, lat, lon):
        """Unparseable or out-of-range values warn instead of crashing."""
        rows = [{'id': 'S1', 'population': 'A', 'latitude': lat, 'longitude': lon}]
        store, warnings = rows_to_metadata_store(rows)

        assert warnings
        assert 'S1' not in store['coordinates']

    def test_half_a_coordinate_is_dropped(self):
        """A latitude with no longitude cannot place anything on a map."""
        rows = [{'id': 'S1', 'population': 'A', 'latitude': '-35.3', 'longitude': ''}]
        store, warnings = rows_to_metadata_store(rows)

        assert warnings == []
        assert store['coordinates'] == {}


class TestMetadataDataTable:
    """The rendered table component."""

    def test_only_metadata_columns_are_editable(self, metadata_store):
        """Sequence IDs and the presence flags are derived, not typed."""
        rows = build_metadata_rows(['S1', 'S2'], metadata_store)
        table = build_metadata_datatable(rows, metadata_store['population_colors'])

        editable = {col['id']: col['editable'] for col in table.columns}
        assert editable == {
            'id': False,
            'in_alignment': False,
            'in_metadata': False,
            'population': True,
            'color': True,
            'latitude': True,
            'longitude': True,
        }

    def test_rows_cannot_be_added_or_deleted(self, metadata_store):
        """Row membership follows the alignment, not the user."""
        table = build_metadata_datatable(
            build_metadata_rows(['S1'], metadata_store), {}
        )

        assert getattr(table, 'row_deletable', False) is False
        assert table.editable is True

    def test_one_swatch_rule_per_colour(self, metadata_store):
        """Colours are tinted cells, since DataTable cannot render HTML."""
        table = build_metadata_datatable([], metadata_store['population_colors'])

        backgrounds = {rule['backgroundColor'] for rule in table.style_data_conditional}
        assert backgrounds == {'#ff0000', '#00ff00'}

    def test_no_colours_means_no_conditional_styles(self):
        """An uncoloured table needs no swatch rules at all."""
        table = build_metadata_datatable([], None)

        assert table.style_data_conditional == []


def _collect_swatches(component, found=None):
    """Recursively collect the pattern-matching colour input ids."""
    if found is None:
        found = []
    cid = getattr(component, 'id', None)
    if isinstance(cid, dict) and cid.get('type') == 'pop-color':
        found.append((cid['pop'], component.value))
    children = getattr(component, 'children', None)
    if isinstance(children, (list, tuple)):
        for child in children:
            _collect_swatches(child, found)
    elif children is not None:
        _collect_swatches(children, found)
    return found


class TestPopulationColorPicker:
    """Per-population colour swatches under the metadata table."""

    def test_one_swatch_per_population(self, metadata_store):
        """Populations, not rows, get a control."""
        rows = build_metadata_rows(['S1', 'S2'], metadata_store)
        swatches = _collect_swatches(build_population_color_controls(rows))

        assert dict(swatches) == {'PopA': '#ff0000', 'PopB': '#00ff00'}

    def test_repeated_population_gets_one_swatch(self):
        """Two rows of one population share a single control."""
        rows = [
            {'id': 'S1', 'population': 'PopA', 'color': '#111111'},
            {'id': 'S2', 'population': 'PopA', 'color': '#111111'},
        ]
        swatches = _collect_swatches(build_population_color_controls(rows))

        assert swatches == [('PopA', '#111111')]

    def test_unnamed_populations_are_skipped(self):
        """A row with no population has nothing to colour."""
        rows = [{'id': 'S1', 'population': '', 'color': ''}]

        assert _collect_swatches(build_population_color_controls(rows)) == []

    def test_uncoloured_population_gets_a_generated_swatch(self):
        """A picker must always open on a real colour, never blank."""
        rows = [{'id': 'S1', 'population': 'PopA', 'color': ''}]
        swatches = _collect_swatches(build_population_color_controls(rows))

        assert len(swatches) == 1
        assert swatches[0][0] == 'PopA'
        assert swatches[0][1].startswith('#')

    def test_colors_from_rows_is_first_wins(self):
        """Matches how rows_to_metadata_store resolves a conflict."""
        rows = [
            {'id': 'S1', 'population': 'PopA', 'color': '#111111'},
            {'id': 'S2', 'population': 'PopA', 'color': '#222222'},
        ]

        assert population_colors_from_rows(rows)['PopA'] == '#111111'


class TestApplyPopulationColor:
    """Recolouring every row of one population."""

    def test_only_matching_rows_change(self, metadata_store):
        """Other populations keep their colour."""
        rows = build_metadata_rows(['S1', 'S2'], metadata_store)
        updated = {r['id']: r for r in apply_population_color(rows, 'PopA', '#123456')}

        assert updated['S1']['color'] == '#123456'
        assert updated['S2']['color'] == '#00ff00'

    def test_all_rows_of_the_population_change(self):
        """The whole point: one pick, every row."""
        rows = [
            {'id': 'S1', 'population': 'PopA', 'color': '#111111'},
            {'id': 'S2', 'population': 'PopA', 'color': '#222222'},
        ]
        updated = apply_population_color(rows, 'PopA', '#abcdef')

        assert [r['color'] for r in updated] == ['#abcdef', '#abcdef']

    def test_input_is_not_mutated(self):
        """The draft diff compares against the originals."""
        rows = [{'id': 'S1', 'population': 'PopA', 'color': '#111111'}]
        apply_population_color(rows, 'PopA', '#abcdef')

        assert rows[0]['color'] == '#111111'

    def test_unknown_population_is_a_no_op(self):
        """A stale swatch must not blank out the table."""
        rows = [{'id': 'S1', 'population': 'PopA', 'color': '#111111'}]

        assert apply_population_color(rows, 'PopZ', '#abcdef') == rows
