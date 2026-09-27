"""
Pure helpers backing the editable metadata table in the GUI.

Kept free of callbacks so the row/store transformations can be tested
without a running Dash app.
"""

from typing import Dict, Iterable, List, Optional, Tuple

from dash import dash_table, html
import dash_bootstrap_components as dbc

from ..visualization.style import (
    POP_BONE,
    POP_INK,
    POP_PAPER,
    generate_population_colors,
    luma,
)

#: Columns the user may type into. Everything else is derived.
EDITABLE_COLUMNS = ('population', 'color', 'latitude', 'longitude')

#: Column order of the metadata table, as (field, display label).
COLUMNS = (
    ('id', 'Sequence ID'),
    ('in_alignment', 'In Alignment'),
    ('in_metadata', 'In Metadata'),
    ('population', 'Population'),
    ('color', 'Color'),
    ('latitude', 'Latitude'),
    ('longitude', 'Longitude'),
)


def build_metadata_rows(
    alignment_ids: Iterable[str], metadata_data: Optional[Dict]
) -> List[Dict[str, str]]:
    """
    Build the metadata table rows for every known sequence.

    Values are all strings so the rows survive a round trip through a
    ``dcc.Store`` unchanged.

    Parameters
    ----------
    alignment_ids : iterable of str
        Sequence IDs present in the uploaded alignment.
    metadata_data : dict, optional
        Contents of the metadata store, if a CSV has been uploaded.

    Returns
    -------
    list of dict
        One row per sequence ID, sorted by ID.
    """
    alignment_ids = set(alignment_ids)
    metadata_data = metadata_data or {}

    metadata_ids = list(metadata_data.get('sequence_ids', []))
    populations = metadata_data.get('populations', {}) or {}
    coordinates = metadata_data.get('coordinates', {}) or {}
    population_colors = metadata_data.get('population_colors', {}) or {}

    rows = []
    for sid in sorted(alignment_ids.union(metadata_ids)):
        coords = coordinates.get(sid, {}) or {}
        pop = populations.get(sid, '') or ''
        rows.append(
            {
                'id': sid,
                'in_alignment': '✓' if sid in alignment_ids else '✗',
                'in_metadata': '✓' if sid in metadata_ids else '✗',
                'population': str(pop),
                'color': str(population_colors.get(pop, '') or ''),
                'latitude': _as_text(coords.get('lat')),
                'longitude': _as_text(coords.get('lon')),
            }
        )
    return rows


def _as_text(value) -> str:
    """
    Render a store value as table text, mapping missing values to ''.

    Parameters
    ----------
    value : object
        Value from the metadata store.

    Returns
    -------
    str
        The value as a string, or '' when absent.
    """
    return '' if value is None or value == '' else str(value)


def build_metadata_datatable(
    rows: List[Dict[str, str]], population_colors: Optional[Dict[str, str]] = None
) -> dash_table.DataTable:
    """
    Build the editable metadata table.

    Only the metadata columns are editable; rows cannot be added or
    removed, since the sequence IDs come from the alignment.

    Parameters
    ----------
    rows : list of dict
        Rows from :func:`build_metadata_rows`.
    population_colors : dict, optional
        Population to hex colour mapping, used to tint the colour cell.
        ``DataTable`` cannot render HTML, so the swatch is a cell
        background rather than a coloured dot.

    Returns
    -------
    dash_table.DataTable
        The configured table component.
    """
    columns = [
        {
            'name': label,
            'id': field,
            'editable': field in EDITABLE_COLUMNS,
        }
        for field, label in COLUMNS
    ]

    style_data_conditional = []
    for color_hex in sorted({c for c in (population_colors or {}).values() if c}):
        style_data_conditional.append(
            {
                'if': {
                    'column_id': 'color',
                    'filter_query': f'{{color}} = "{color_hex}"',
                },
                'backgroundColor': color_hex,
                'color': _readable_text_color(color_hex),
            }
        )

    return dash_table.DataTable(
        id='metadata-table',
        data=rows,
        columns=columns,
        editable=True,
        sort_action='native',
        filter_action='native',
        page_action='native',
        page_size=50,
        style_table={'overflowX': 'auto'},
        style_cell={
            'fontSize': '13px',
            'padding': '6px 8px',
            'textAlign': 'left',
            'fontFamily': 'inherit',
            'border': f'1px solid {POP_BONE}',
        },
        style_header={
            'fontWeight': 'bold',
            'backgroundColor': POP_INK,
            'color': POP_PAPER,
            'border': f'1px solid {POP_INK}',
        },
        style_data_conditional=style_data_conditional,
    )


def _readable_text_color(hex_color: str) -> str:
    """
    Pick black or white text for a swatch background.

    Parameters
    ----------
    hex_color : str
        Background colour as ``#rrggbb``.

    Returns
    -------
    str
        ``'#000000'`` or ``'#ffffff'``.
    """
    # Shared with the node-brightness check in style.py.
    return '#000000' if luma(hex_color) > 140 else '#ffffff'


def diff_metadata_rows(
    committed: List[Dict[str, str]], edited: List[Dict[str, str]]
) -> List[Tuple[str, str]]:
    """
    List the cells that differ between two sets of rows.

    Rows are matched on their ``id``, never on position: the table has
    native sort and filter, so the order the user sees is not the order
    the rows were built in.

    Parameters
    ----------
    committed : list of dict
        Rows as last applied to the figure.
    edited : list of dict
        Rows currently in the table.

    Returns
    -------
    list of tuple of str
        ``(sequence_id, column)`` pairs that changed.
    """
    baseline = {row['id']: row for row in committed}

    changes = []
    for row in edited:
        original = baseline.get(row['id'])
        if original is None:
            continue
        for column in EDITABLE_COLUMNS:
            if _as_text(row.get(column)) != _as_text(original.get(column)):
                changes.append((row['id'], column))
    return changes


def rows_to_metadata_store(
    rows: List[Dict[str, str]], previous: Optional[Dict] = None
) -> Tuple[Dict, List[str]]:
    """
    Convert edited table rows back into the metadata store shape.

    The output matches exactly what the CSV upload callback writes, so
    every downstream consumer keeps working unchanged.

    Parameters
    ----------
    rows : list of dict
        Rows from the metadata table.
    previous : dict, optional
        The metadata store being replaced. Unused keys are carried over.

    Returns
    -------
    tuple of (dict, list of str)
        The new store contents and any warnings raised while parsing.
    """
    warnings: List[str] = []
    raw: Dict[str, Dict[str, str]] = {}
    populations: Dict[str, str] = {}
    coordinates: Dict[str, Dict[str, float]] = {}
    population_colors: Dict[str, str] = {}

    for row in rows:
        sid = row.get('id')
        if not sid:
            continue

        pop = _as_text(row.get('population')).strip()
        color = _as_text(row.get('color')).strip()

        entry: Dict[str, str] = {}
        if pop:
            populations[sid] = pop
            entry['population'] = pop
            # First value wins, matching the CSV upload path.
            if color and pop not in population_colors:
                population_colors[pop] = color
        if color:
            entry['color'] = color

        lat = _parse_coordinate(row.get('latitude'), -90.0, 90.0, sid, warnings, 'lat')
        lon = _parse_coordinate(
            row.get('longitude'), -180.0, 180.0, sid, warnings, 'lon'
        )
        if lat is not None and lon is not None:
            coordinates[sid] = {'lat': lat, 'lon': lon}
            entry['latitude'] = str(lat)
            entry['longitude'] = str(lon)

        if entry:
            raw[sid] = entry

    # Fill in any population left without a colour.
    uncoloured = sorted(set(populations.values()) - set(population_colors))
    if uncoloured:
        generated = generate_population_colors(sorted(set(populations.values())))
        for pop in uncoloured:
            population_colors[pop] = generated[pop]

    store = dict(previous or {})
    store.update(
        {
            'raw': raw,
            'coordinates': coordinates,
            'populations': populations,
            'population_colors': population_colors if population_colors else None,
            'sequence_ids': list(raw.keys()),
        }
    )
    return store, warnings


def _parse_coordinate(
    value,
    low: float,
    high: float,
    sid: str,
    warnings: List[str],
    kind: str,
) -> Optional[float]:
    """
    Parse one coordinate cell, recording a warning when it is unusable.

    Parameters
    ----------
    value : object
        Raw cell contents.
    low : float
        Lowest valid value.
    high : float
        Highest valid value.
    sid : str
        Sequence ID, used in the warning text.
    warnings : list of str
        Accumulator for warning messages, appended to in place.
    kind : str
        ``'lat'`` or ``'lon'``, used in the warning text.

    Returns
    -------
    float or None
        The parsed value, or ``None`` if blank or invalid.
    """
    text = _as_text(value).strip()
    if not text:
        return None
    try:
        number = float(text)
    except ValueError:
        warnings.append(f'{sid}: {kind} "{text}" is not a number')
        return None
    if not low <= number <= high:
        warnings.append(f'{sid}: {kind} {number} is outside [{low}, {high}]')
        return None
    return number


def population_colors_from_rows(
    rows: List[Dict[str, str]],
) -> Dict[str, str]:
    """
    Read the colour currently assigned to each population.

    Colours are stored per population, but the table holds one per row, so
    the first non-empty value for a population wins -- the same rule
    :func:`rows_to_metadata_store` applies on commit.

    Parameters
    ----------
    rows : list of dict
        Rows from the metadata table.

    Returns
    -------
    dict
        Population name to hex colour, for every named population.
    """
    colors: Dict[str, str] = {}
    for row in rows or []:
        pop = _as_text(row.get('population')).strip()
        if not pop:
            continue
        color = _as_text(row.get('color')).strip()
        if color and pop not in colors:
            colors[pop] = color

    missing = sorted(
        {
            _as_text(row.get('population')).strip()
            for row in rows or []
            if _as_text(row.get('population')).strip()
        }
        - set(colors)
    )
    if missing:
        generated = generate_population_colors(sorted(set(missing) | set(colors)))
        for pop in missing:
            colors[pop] = generated[pop]
    return colors


def build_population_color_controls(rows: List[Dict[str, str]]) -> html.Div:
    """
    Build one colour swatch per population.

    A native colour input is used rather than a per-cell editor because
    colours belong to the population, not the row: editing them per row
    lets two rows of one population disagree, and the disagreement is then
    silently resolved first-wins on commit.

    Parameters
    ----------
    rows : list of dict
        Rows from the metadata table.

    Returns
    -------
    html.Div
        A swatch per population, or an empty div when none are named.
    """
    colors = population_colors_from_rows(rows)
    if not colors:
        return html.Div()

    swatches = [
        html.Div(
            [
                dbc.Input(
                    type='color',
                    id={'type': 'pop-color', 'pop': pop},
                    value=color,
                    className='pp-color-swatch',
                ),
                html.Span(pop, className='fw-bold ms-2'),
                html.Small(color, className='text-muted ms-2'),
            ],
            className='d-flex align-items-center me-4 mb-2',
        )
        for pop, color in sorted(colors.items())
    ]

    return html.Div(
        [
            html.Strong('Population colours'),
            html.Small(
                ' — a swatch recolours every row of that population.',
                className='text-muted',
            ),
            html.Div(swatches, className='d-flex flex-wrap mt-2'),
        ],
        className='pp-population-colors',
    )


def apply_population_color(
    rows: List[Dict[str, str]], population: str, hex_color: str
) -> List[Dict[str, str]]:
    """
    Recolour every row belonging to one population.

    Parameters
    ----------
    rows : list of dict
        Rows from the metadata table.
    population : str
        Population whose rows are recoloured.
    hex_color : str
        New colour as ``#rrggbb``.

    Returns
    -------
    list of dict
        New rows; the originals are left untouched.
    """
    updated = []
    for row in rows or []:
        new_row = dict(row)
        if _as_text(new_row.get('population')).strip() == population:
            new_row['color'] = hex_color
        updated.append(new_row)
    return updated
