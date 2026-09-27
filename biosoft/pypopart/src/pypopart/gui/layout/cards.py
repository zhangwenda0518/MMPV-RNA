"""Layout builders for the PyPopART GUI: page skeleton, cards and tabs."""

from dash import dcc, html
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto

from ...visualization.cytoscape_plot import DEFAULT_TICK_THRESHOLD, MAX_TICK_MARKS
from ..callbacks.feedback import TOAST_DURATION_MS


def build_layout(app) -> None:
    """
    Set up the application layout with all components.

    Parameters
    ----------
    app : dash.Dash
        The Dash application to attach the layout to.
    """
    app.layout = html.Div(
        [
            # Header
            html.Div(
                [
                    html.H1(
                        [
                            html.Span('Py', className='pp-t1'),
                            html.Span('Pop', className='pp-t2'),
                            html.Span('ART', className='pp-t3'),
                        ],
                        className='pp-title',
                    ),
                    # Shown by the long callbacks' `running=` spec while
                    # they work, so the user can see what the app is doing.
                    html.Div(
                        [
                            dbc.Spinner(
                                size='sm',
                                color='dark',
                                spinner_class_name='pp-task-spinner',
                            ),
                            html.Span(id='task-label', className='pp-task-label'),
                        ],
                        id='task-indicator',
                        className='pp-task',
                        style={'display': 'none'},
                    ),
                    html.Span(
                        'Haplotype Network Analysis',
                        className='pp-subtitle',
                    ),
                ],
                className='pp-titlebar',
            ),
            # Main resizable container
            html.Div(
                [
                    # Left panel - Controls (resizable sidebar)
                    html.Div(
                        [
                            create_upload_card(),
                            html.Br(),
                            create_algorithm_card(),
                            html.Br(),
                            create_layout_card(),
                            html.Br(),
                            create_export_card(),
                        ],
                        id='sidebar-panel',
                        className='pp-sidebar',
                        # The width is dragged by assets/pypopart.js, which
                        # writes it straight onto this node's inline style.
                        # 'minWidth'/'maxWidth' are the drag clamps.
                        style={
                            'minWidth': '250px',
                            'width': '300px',
                            'maxWidth': '600px',
                            'height': '90vh',
                            'overflowY': 'auto',
                            'padding': '20px',
                        },
                    ),
                    # Collapse button at the top, drag grip at the centre of
                    # the edge. The browser's native 'resize' corner sat at
                    # the bottom of a 90vh scrolling panel and was almost
                    # always off-screen.
                    html.Div(
                        [
                            dbc.Button(
                                '«',
                                id='sidebar-toggle',
                                color='light',
                                size='sm',
                                title='Hide the control panel',
                                className='border pp-sidebar-toggle',
                            ),
                            html.Div(
                                '⋮',
                                id='sidebar-grip',
                                className='pp-sidebar-grip',
                                title='Drag to resize the control panel',
                            ),
                        ],
                        id='sidebar-resizer',
                        className='pp-sidebar-resizer',
                    ),
                    # Right panel - Visualization
                    html.Div(
                        [
                            dbc.Tabs(
                                [
                                    dbc.Tab(
                                        create_network_tab(),
                                        label='Network',
                                    ),
                                    dbc.Tab(
                                        create_statistics_tab(),
                                        label='Statistics',
                                    ),
                                    dbc.Tab(
                                        create_haplotype_summary_tab(),
                                        label='Haplotype Summary',
                                    ),
                                    dbc.Tab(
                                        create_metadata_tab(),
                                        label='Metadata',
                                    ),
                                    dbc.Tab(
                                        create_alignment_tab(),
                                        label='Alignment',
                                    ),
                                ]
                            )
                        ],
                        id='main-panel',
                        className='pp-main',
                        style={
                            'flex': '1',
                            'padding': '20px',
                            'overflow': 'auto',
                        },
                    ),
                ],
                style={
                    'display': 'flex',
                    'height': '90vh',
                    'overflow': 'hidden',
                },
            ),
            # Transient confirmations. Errors stay inline in their panel.
            dbc.Toast(
                id='app-toast',
                header='',
                is_open=False,
                dismissable=True,
                duration=TOAST_DURATION_MS,
                className='pp-toast',
            ),
            # Hidden stores for data
            dcc.Store(id='alignment-store'),
            dcc.Store(id='metadata-store'),
            dcc.Store(id='network-store'),
            dcc.Store(id='layout-store'),
            dcc.Store(id='computation-status'),
            dcc.Store(id='geographic-mode', data=False),
            dcc.Store(id='manual-edit-flag', data=False),
            # Manually dragged node positions, kept apart from
            # layout-store so a drag does not rebuild the whole graph.
            dcc.Store(id='node-positions-store'),
            dcc.Store(id='sidebar-collapsed', data=False),
            # Uncommitted metadata edits: {'rows': [...], 'edited': [[id, col]]}
            dcc.Store(id='metadata-draft-store'),
            # Bumped once metadata-store is written, which is what actually
            # triggers the network computation. See callbacks/metadata.py.
            dcc.Store(id='metadata-commit-token'),
            # Store to trigger window resize handling
            dcc.Store(id='window-size-store'),
        ]
    )


def create_upload_card() -> dbc.Card:
    """
    Create file upload card.

    Returns
    -------
    dbc.Card
        The file upload card.
    """
    return dbc.Card(
        [
            dbc.CardHeader(
                html.H5('1. Upload Data', className='mb-0'),
            ),
            dbc.CardBody(
                [
                    dbc.Label('Sequence File', className='fw-bold'),
                    html.Small(
                        'Upload aligned sequences in FASTA, NEXUS, or PHYLIP format',
                        className='text-muted d-block mb-2',
                    ),
                    dcc.Upload(
                        id='upload-data',
                        children=html.Div(
                            [
                                dbc.Button(
                                    'Select Sequence File',
                                    color='primary',
                                    className='w-100',
                                ),
                                html.Small(
                                    'or drop a file here',
                                    className='pp-dropzone-hint',
                                ),
                            ]
                        ),
                        multiple=False,
                        className='pp-dropzone',
                        className_active='pp-dropzone--active',
                    ),
                    html.Div(id='upload-status', className='mt-2'),
                    html.Hr(),
                    dbc.Label('Metadata File (Optional)', className='fw-bold'),
                    html.Small(
                        'CSV file with population, location, or trait data',
                        className='text-muted d-block mb-2',
                    ),
                    dcc.Upload(
                        id='upload-metadata',
                        children=html.Div(
                            [
                                dbc.Button(
                                    'Select Metadata File',
                                    color='secondary',
                                    outline=True,
                                    className='w-100',
                                ),
                                html.Small(
                                    'or drop a CSV here',
                                    className='pp-dropzone-hint',
                                ),
                            ]
                        ),
                        multiple=False,
                        className='pp-dropzone',
                        className_active='pp-dropzone--active',
                    ),
                    html.Div(id='metadata-status', className='mt-2'),
                    html.Div(
                        id='metadata-template-section',
                        children=[
                            html.Hr(),
                            dbc.Button(
                                'Download Metadata Template',
                                id='download-template-button',
                                color='info',
                                outline=True,
                                size='sm',
                                className='w-100',
                                disabled=True,
                            ),
                            dcc.Download(id='download-template'),
                            html.Small(
                                'Get a CSV template pre-filled with your sequence IDs',
                                className='text-muted d-block mt-1',
                            ),
                        ],
                    ),
                ]
            ),
        ]
    )


def create_algorithm_card() -> dbc.Card:
    """
    Create algorithm selection and parameter card.

    Returns
    -------
    dbc.Card
        The algorithm selection and parameter card.
    """
    return dbc.Card(
        [
            dbc.CardHeader(
                html.H5('2. Configure Algorithm', className='mb-0'),
            ),
            dbc.CardBody(
                [
                    dbc.Label('Network Algorithm', className='fw-bold'),
                    html.Small(
                        'Choose the method for constructing the haplotype network',
                        className='text-muted d-block mb-2',
                    ),
                    dcc.Dropdown(
                        id='algorithm-select',
                        options=[
                            {
                                'label': 'MST - Minimum Spanning Tree',
                                'value': 'mst',
                            },
                            {
                                'label': 'MSN - Minimum Spanning Network',
                                'value': 'msn',
                            },
                            {
                                'label': 'TCS - Statistical Parsimony',
                                'value': 'tcs',
                            },
                            {
                                'label': 'MJN - Median-Joining Network',
                                'value': 'mjn',
                            },
                            {
                                'label': 'PN - Parsimony Network',
                                'value': 'pn',
                            },
                            {
                                'label': 'TSW - Tight Span Walker (Parsimony Network)',
                                'value': 'tsw',
                            },
                        ],
                        value='mst',
                        style={'whiteSpace': 'nowrap'},
                    ),
                    html.Br(),
                    html.Div(id='algorithm-parameters'),
                    html.Br(),
                    html.Div(id='metadata-draft-badge-sidebar', className='mb-2'),
                    dbc.Button(
                        'Compute Network',
                        id='compute-button',
                        color='success',
                        className='w-100',
                        disabled=True,
                    ),
                    html.Div(id='computation-feedback', className='mt-2'),
                ]
            ),
        ]
    )


def create_layout_card() -> dbc.Card:
    """
    Create layout configuration card.

    Returns
    -------
    dbc.Card
        The layout configuration card.
    """
    return dbc.Card(
        [
            dbc.CardHeader(
                html.H5('3. Layout Options', className='mb-0'),
            ),
            dbc.CardBody(
                [
                    dbc.Label('Layout Algorithm', className='fw-bold'),
                    html.Small(
                        'Choose how to position nodes in the visualization',
                        className='text-muted d-block mb-2',
                    ),
                    dcc.Dropdown(
                        id='layout-select',
                        options=[
                            {
                                'label': 'Hierarchical (Fast)',
                                'value': 'hierarchical',
                            },
                            {
                                'label': 'Spring (Force-directed)',
                                'value': 'spring',
                            },
                            {
                                'label': 'Spring - Proportional Edge Length',
                                'value': 'spring_proportional',
                            },
                            {
                                'label': 'Spectral (Fast, Large networks)',
                                'value': 'spectral',
                            },
                            {'label': 'Circular', 'value': 'circular'},
                            {'label': 'Radial', 'value': 'radial'},
                            {
                                'label': 'Kamada-Kawai (High quality, slow)',
                                'value': 'kamada_kawai',
                            },
                            {
                                'label': 'Kamada-Kawai - Proportional Edge Length',
                                'value': 'kamada_kawai_proportional',
                            },
                            #                                {
                            #'label': 'Geographic (requires coordinates)',
                            #'value': 'geographic',
                            #                                },
                        ],
                        value='spring',
                        style={'whiteSpace': 'nowrap'},
                    ),
                    html.Br(),
                    html.Div(
                        id='geographic-options',
                        children=[
                            dbc.Label('Map Projection'),
                            dcc.Dropdown(
                                id='map-projection',
                                options=[
                                    {'label': 'Mercator', 'value': 'mercator'},
                                    {
                                        'label': 'PlateCarree',
                                        'value': 'platecarree',
                                    },
                                    {
                                        'label': 'Orthographic',
                                        'value': 'orthographic',
                                    },
                                ],
                                value='mercator',
                            ),
                            html.Br(),
                            dbc.Label('Zoom Level'),
                            dcc.Slider(
                                id='map-zoom',
                                min=1,
                                max=10,
                                step=1,
                                value=2,
                                marks={i: str(i) for i in range(1, 11)},
                            ),
                        ],
                        style={'display': 'none'},
                    ),
                    html.Br(),
                    dbc.Label('Node Spacing', className='fw-bold'),
                    html.Small(
                        'Adjust the spacing between nodes',
                        className='text-muted d-block mb-2',
                    ),
                    dcc.Slider(
                        id='spacing-slider',
                        min=0.5,
                        max=3.0,
                        step=0.1,
                        value=2.0,
                        marks={0.5: '0.5x', 1.0: '1.0x', 2.0: '2.0x', 3.0: '3.0x'},
                        tooltip={'placement': 'bottom', 'always_visible': False},
                    ),
                    html.Div(
                        id='repel-options',
                        children=[
                            html.Br(),
                            dbc.Label('Repel Force', className='fw-bold'),
                            html.Small(
                                'How strongly nodes push apart in the spring layouts',
                                className='text-muted d-block mb-2',
                            ),
                            dcc.Slider(
                                id='repel-slider',
                                min=0.5,
                                max=4.0,
                                step=0.1,
                                value=1.0,
                                marks={0.5: '0.5x', 1.0: '1x', 2.0: '2x', 4.0: '4x'},
                                tooltip={
                                    'placement': 'bottom',
                                    'always_visible': False,
                                },
                            ),
                        ],
                        style={'display': 'none'},
                    ),
                    html.Br(),
                    dbc.Label('Node Size', className='fw-bold'),
                    html.Small(
                        'Adjust the size of nodes',
                        className='text-muted d-block mb-2',
                    ),
                    dcc.Slider(
                        id='node-size-slider',
                        min=10,
                        max=100,
                        step=5,
                        value=40,
                        marks={10: '10', 40: '40', 70: '70', 100: '100'},
                        tooltip={'placement': 'bottom', 'always_visible': False},
                    ),
                    html.Br(),
                    dbc.Label('Edge Width', className='fw-bold'),
                    html.Small(
                        'Adjust the width of edges',
                        className='text-muted d-block mb-2',
                    ),
                    dcc.Slider(
                        id='edge-width-slider',
                        min=1,
                        max=10,
                        step=0.5,
                        value=3,
                        marks={1: '1', 3: '3', 6: '6', 10: '10'},
                        tooltip={'placement': 'bottom', 'always_visible': False},
                    ),
                    html.Br(),
                    dbc.Label('Mutation Marks', className='fw-bold'),
                    html.Small(
                        'Draw one dash across an edge per mutation',
                        className='text-muted d-block mb-2',
                    ),
                    dbc.Switch(
                        id='edge-tick-toggle',
                        label='Show tick marks',
                        value=True,
                        className='mb-2',
                    ),
                    dbc.Label('Use numerals above', className='fw-bold'),
                    html.Small(
                        'Edges with more mutations than this show a number',
                        className='text-muted d-block mb-2',
                    ),
                    dcc.Slider(
                        id='edge-tick-threshold',
                        min=1,
                        max=MAX_TICK_MARKS,
                        step=1,
                        value=DEFAULT_TICK_THRESHOLD,
                        marks={1: '1', 10: '10', 20: '20', MAX_TICK_MARKS: '30'},
                        tooltip={'placement': 'bottom', 'always_visible': False},
                    ),
                    html.Br(),
                    dbc.Label('Selection', className='fw-bold'),
                    html.Small(
                        'Click a node to select it, Cmd/Ctrl+click to add more, '
                        'click the background to clear. Dragging any selected '
                        'node moves the whole selection. Shift+drag draws a '
                        'selection box; switch on box select to draw one with '
                        'a plain drag instead of panning.',
                        className='text-muted d-block mb-2',
                    ),
                    dbc.Switch(
                        id='box-select-toggle',
                        label='Box select (drag to select, not pan)',
                        value=False,
                        className='mb-2',
                    ),
                    html.Br(),
                    dbc.Label('Snap to Grid', className='fw-bold'),
                    html.Small(
                        'Align nodes to a grid when dragged or laid out',
                        className='text-muted d-block mb-2',
                    ),
                    dbc.Switch(
                        id='snap-to-grid-toggle',
                        label='Snap to grid',
                        value=False,
                        className='mb-2',
                    ),
                    dcc.Slider(
                        id='grid-size',
                        min=25,
                        max=200,
                        step=25,
                        value=50,
                        marks={25: '25', 50: '50', 100: '100', 200: '200'},
                        tooltip={'placement': 'bottom', 'always_visible': False},
                    ),
                    html.Br(),
                    dbc.Button(
                        'Apply Layout',
                        id='apply-layout-button',
                        color='info',
                        className='w-100',
                        disabled=True,
                    ),
                ]
            ),
        ]
    )


def create_export_card() -> dbc.Card:
    """
    Create export options card.

    Returns
    -------
    dbc.Card
        The export options card.
    """
    return dbc.Card(
        [
            dbc.CardHeader(
                html.H5('4. Export', className='mb-0'),
            ),
            dbc.CardBody(
                [
                    dbc.Label('Export Format', className='fw-bold'),
                    html.Small(
                        'Save your network for further analysis or publication',
                        className='text-muted d-block mb-2',
                    ),
                    dcc.Dropdown(
                        id='export-format',
                        options=[
                            {
                                'label': 'GraphML (Cytoscape/Gephi)',
                                'value': 'graphml',
                            },
                            {'label': 'GML (Graph Format)', 'value': 'gml'},
                            {'label': 'JSON', 'value': 'json'},
                            {'label': 'PNG Image (current view)', 'value': 'png'},
                            {'label': 'SVG Figure (publication)', 'value': 'svg'},
                        ],
                        value='png',
                        style={'whiteSpace': 'nowrap'},
                    ),
                    html.Br(),
                    dbc.Checklist(
                        id='export-legend',
                        options=[
                            {
                                'label': 'Include population legend',
                                'value': 'legend',
                            }
                        ],
                        value=['legend'],
                        switch=True,
                        className='mb-2',
                    ),
                    html.Small(
                        'Drawn into SVG figures only; the PNG captures the '
                        'on-screen legend already.',
                        className='text-muted d-block mb-2',
                    ),
                    dbc.Button(
                        'Download',
                        id='export-button',
                        color='secondary',
                        className='w-100',
                        disabled=True,
                    ),
                    dcc.Download(id='download-data'),
                ]
            ),
        ]
    )


def create_network_tab() -> html.Div:
    """
    Create network visualization tab.

    Returns
    -------
    html.Div
        The network visualization tab.
    """
    return html.Div(
        [
            # Search bar
            html.Div(
                [
                    html.Label(
                        'Search Haplotype:',
                        style={'marginRight': '10px', 'fontWeight': 'bold'},
                    ),
                    dcc.Dropdown(
                        id='haplotype-search',
                        placeholder='Select H number(s) to highlight...',
                        style={'width': '300px', 'display': 'inline-block'},
                        clearable=True,
                        multi=True,
                    ),
                ],
                className='pp-overlay pp-search',
            ),
            # Legend display
            html.Div(
                id='network-legend',
                className='pp-overlay pp-legend',
            ),
            # Tooltip display on hover. Presentation lives in theme.css so
            # the clientside renderer only has to toggle 'display'.
            html.Div(id='node-tooltip'),
            dcc.Loading(
                id='loading-network',
                type='default',
                children=[
                    cyto.Cytoscape(
                        id='network-graph',
                        layout={'name': 'preset'},
                        style={'width': '100%', 'height': '85vh'},
                        elements=[],
                        stylesheet=[],
                        minZoom=0.1,
                        maxZoom=5,
                        wheelSensitivity=0.2,
                        zoom=1,
                        autoungrabify=False,
                        # Shift+drag draws a selection box; the box-select
                        # switch turns panning off so a plain drag does too.
                        boxSelectionEnabled=True,
                        userPanningEnabled=True,
                    )
                ],
            ),
        ],
        style={'position': 'relative'},
    )


def create_statistics_tab() -> html.Div:
    """
    Create statistics display tab.

    Returns
    -------
    html.Div
        The statistics display tab.
    """
    return html.Div(
        [
            dcc.Loading(
                html.Div(
                    id='statistics-display',
                    style={'padding': '20px', 'height': '85vh', 'overflow-y': 'auto'},
                )
            )
        ]
    )


def create_alignment_tab() -> html.Div:
    """
    Create alignment viewer tab.

    Returns
    -------
    html.Div
        The alignment viewer tab.
    """
    return html.Div(
        [
            html.Div(
                id='alignment-display',
                className='pp-alignment',
                style={
                    'padding': '20px',
                    'height': '85vh',
                    'overflow': 'auto',
                    'whiteSpace': 'pre',
                },
            )
        ]
    )


def create_haplotype_summary_tab() -> html.Div:
    """
    Create haplotype summary tab showing H number to sequence name mapping.

    Returns
    -------
    html.Div
        The haplotype summary tab showing H number to sequence name mapping.
    """
    return html.Div(
        [
            html.Div(
                [
                    dbc.Button(
                        'Download Summary CSV',
                        id='download-haplotype-csv-button',
                        color='primary',
                        size='sm',
                    ),
                    dcc.Download(id='download-haplotype-csv'),
                    dbc.Button(
                        'Download Label Template',
                        id='download-h-number-template-button',
                        color='info',
                        outline=True,
                        size='sm',
                    ),
                    dcc.Download(id='download-h-number-template'),
                    dcc.Upload(
                        id='upload-h-number-mapping',
                        children=dbc.Button(
                            'Upload Label Mapping',
                            color='warning',
                            outline=True,
                            size='sm',
                        ),
                    ),
                ],
                # Flex row rather than inline-block plus margins, so the
                # buttons wrap cleanly when the panel is narrow.
                className='d-flex flex-wrap align-items-center gap-2',
                style={'padding': '20px 20px 10px 20px'},
            ),
            html.Div(
                id='h-number-feedback',
                style={'padding': '0 20px 10px 20px'},
            ),
            html.Div(
                id='haplotype-summary-display',
                style={
                    'padding': '0 20px 20px 20px',
                    'height': '75vh',
                    'overflow-y': 'auto',
                },
            ),
            dcc.Store(id='h-number-mapping-store'),
        ]
    )


def create_metadata_tab() -> html.Div:
    """
    Create metadata tab showing imported metadata and alignment IDs.

    Returns
    -------
    html.Div
        The metadata tab showing imported metadata and alignment IDs.
    """
    return html.Div(
        [
            html.Div(
                [
                    html.Small(
                        'Population, colour and coordinate cells are editable. '
                        'Edits stay as a draft until you click Compute Network.',
                        className='text-muted d-block mb-2',
                    ),
                    html.Div(id='metadata-draft-badge'),
                    html.Div(id='metadata-edit-feedback'),
                ],
                style={'padding': '20px 20px 0 20px'},
            ),
            html.Div(
                id='metadata-warnings',
                style={'padding': '10px 20px'},
            ),
            html.Div(
                id='population-colors',
                style={'padding': '0 20px 10px 20px'},
            ),
            html.Div(
                id='metadata-display',
                style={
                    'padding': '0 20px 20px 20px',
                    'height': '75vh',
                    'overflow-y': 'auto',
                },
            ),
        ]
    )
