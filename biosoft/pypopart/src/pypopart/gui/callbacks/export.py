"""Network and CSV export downloads."""

import io
from pathlib import Path
import tempfile
import traceback
from typing import Dict, Optional, Tuple

import dash
from dash import Input, Output, State
from dash.exceptions import PreventUpdate

from pypopart.core.graph import HaplotypeNetwork
from pypopart.gui.callbacks.feedback import task_running, toast
from pypopart.gui.serialization import merge_node_positions
from pypopart.io.network_export import GMLExporter, GraphMLExporter, JSONExporter
from pypopart.visualization.cytoscape_plot import (
    DEFAULT_TICK_THRESHOLD,
    resolve_population_counts,
)
from pypopart.visualization.style import POP_PAPER

#: Formats written server-side, as ``{value: (exporter, suffix, mimetype)}``.
TEXT_EXPORTERS = {
    'graphml': (GraphMLExporter, 'graphml', 'text/xml'),
    'gml': (GMLExporter, 'gml', 'text/plain'),
    'json': (JSONExporter, 'json', 'application/json'),
}

#: Rendered in the browser by dash-cytoscape, straight off the canvas.
CANVAS_FORMATS = ('png',)

#: Rendered server-side with matplotlib. SVG is here rather than on the
#: canvas because Cytoscape.js only emits vector output through the
#: cytoscape-svg extension, which dash-cytoscape does not bundle -- asking
#: the canvas for an SVG silently produced no file at all.
FIGURE_FORMATS = ('svg',)

#: Every image format the callback can serve.
IMAGE_FORMATS = CANVAS_FORMATS + FIGURE_FORMATS


def _default_h_numbers(network: HaplotypeNetwork) -> Dict[str, str]:
    """
    Label haplotypes the way the interactive view numbers them.

    Mirrors the fallback the display callbacks use when no custom label
    mapping has been uploaded: H1, H2, ... over the sorted node IDs.

    Parameters
    ----------
    network : HaplotypeNetwork
        Network being labelled.

    Returns
    -------
    Dict[str, str]
        Node ID to H number.
    """
    return {
        node: f'H{i}' for i, node in enumerate(sorted(network.graph.nodes()), start=1)
    }


def resolve_node_colors(
    network: HaplotypeNetwork,
    population_colors: Optional[Dict],
    population_mapping: Optional[Dict],
) -> Dict[str, str]:
    """
    Colour each node the way the interactive view colours it.

    The static plotter resolves colours from the haplotype's own
    population counts, but a network rebuilt from the store often has
    none -- the metadata is uploaded separately and lives in its own
    store. Resolving here from the same sample-to-population mapping the
    canvas uses keeps an exported figure in step with the screen.

    Mixed-population nodes are drawn as a pie on the canvas; a figure
    gets the dominant population's colour instead.

    Parameters
    ----------
    network : HaplotypeNetwork
        Network being drawn.
    population_colors : Dict, optional
        Population to hex colour mapping.
    population_mapping : Dict, optional
        Sample ID to population mapping.

    Returns
    -------
    Dict[str, str]
        Node ID to hex colour, for nodes with a resolvable population.
    """
    if not population_colors:
        return {}

    colors = {}
    for node in network.graph.nodes():
        if network.is_median_vector(node):
            continue
        hap = network.get_haplotype(node)
        if hap is None:
            continue
        counts = resolve_population_counts(hap, population_mapping)
        counts.pop('Unassigned', None)
        if not counts:
            continue
        dominant = max(counts.items(), key=lambda item: item[1])[0]
        if dominant in population_colors:
            colors[node] = population_colors[dominant]
    return colors


def _export_figure(
    network: HaplotypeNetwork,
    positions: Dict,
    population_colors: Optional[Dict],
    export_format: str,
    node_labels: Optional[Dict] = None,
    population_mapping: Optional[Dict] = None,
    show_edge_ticks: bool = True,
    edge_tick_threshold: int = DEFAULT_TICK_THRESHOLD,
    include_legend: bool = False,
) -> Dict:
    """
    Render the network to a vector figure with matplotlib.

    Parameters
    ----------
    network : HaplotypeNetwork
        Network to draw.
    positions : Dict
        Node positions, so the figure matches what is on screen.
    population_colors : Dict, optional
        Population to hex colour mapping.
    export_format : str
        A matplotlib-supported format, currently only ``'svg'``.
    node_labels : Dict, optional
        H number labels, so the figure is labelled like the screen.
    population_mapping : Dict, optional
        Sample ID to population mapping, used to colour nodes.
    show_edge_ticks : bool, default=True
        Draw mutation counts as tick marks rather than numerals.
    edge_tick_threshold : int, default=DEFAULT_TICK_THRESHOLD
        Above this many mutations a numeral is drawn instead.
    include_legend : bool, default=False
        Draw a key mapping colours to population names.

    Returns
    -------
    Dict
        The ``dcc.Download`` payload.
    """
    import matplotlib

    # There is no display attached to the server process.
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    from pypopart.visualization.static_plot import StaticNetworkPlotter

    plotter = StaticNetworkPlotter(network)
    figure, _ = plotter.plot(
        layout=positions or None,
        population_colors=population_colors or None,
        node_color_map=resolve_node_colors(
            network, population_colors, population_mapping
        ),
        node_labels=node_labels or None,
        show_edge_ticks=show_edge_ticks,
        edge_tick_threshold=edge_tick_threshold,
        show_title=False,
    )
    if include_legend and population_colors:
        plotter.add_legend(
            population_colors=population_colors,
            show_size_scale=False,
        )
    try:
        buffer = io.StringIO()
        # Exports get a true white ground even though the app sits on
        # parchment, since figures end up on a printed page.
        figure.savefig(
            buffer,
            format=export_format,
            bbox_inches='tight',
            facecolor=POP_PAPER,
        )
        content = buffer.getvalue()
    finally:
        plt.close(figure)

    return {
        'content': content,
        'filename': f'network.{export_format}',
        'type': 'image/svg+xml',
    }


def _export_text(network: HaplotypeNetwork, export_format: str) -> Dict:
    """
    Render a network to text and wrap it for ``dcc.Download``.

    Parameters
    ----------
    network : HaplotypeNetwork
        Network to export.
    export_format : str
        A key of :data:`TEXT_EXPORTERS`.

    Returns
    -------
    Dict
        The ``dcc.Download`` payload: content, filename and MIME type.
    """
    exporter_cls, suffix, mimetype = TEXT_EXPORTERS[export_format]
    filename = f'network.{suffix}'

    # The exporters write to a path, so stage in a temporary directory.
    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / filename
        exporter_cls(path).export(network)
        content = path.read_text()

    return {'content': content, 'filename': filename, 'type': mimetype}


def register(app, logger) -> None:
    """
    Register export callbacks on the Dash app.

    Parameters
    ----------
    app : dash.Dash
        The Dash application.
    logger : logging.Logger
        Application logger.
    """

    @app.callback(
        [
            Output('download-data', 'data'),
            Output('network-graph', 'generateImage', allow_duplicate=True),
            Output('app-toast', 'children', allow_duplicate=True),
            Output('app-toast', 'header', allow_duplicate=True),
            Output('app-toast', 'is_open', allow_duplicate=True),
        ],
        Input('export-button', 'n_clicks'),
        [
            State('network-store', 'data'),
            State('export-format', 'value'),
            State('layout-store', 'data'),
            State('node-positions-store', 'data'),
            State('metadata-store', 'data'),
            State('h-number-mapping-store', 'data'),
            State('edge-tick-toggle', 'value'),
            State('edge-tick-threshold', 'value'),
            State('export-legend', 'value'),
        ],
        prevent_initial_call=True,
        running=task_running('Exporting'),
    )
    def export_network(
        n_clicks: int,
        network_data: Dict,
        export_format: str,
        layout_data: Optional[Dict] = None,
        dragged_positions: Optional[Dict] = None,
        metadata_data: Optional[Dict] = None,
        h_number_mapping: Optional[Dict] = None,
        show_edge_ticks: Optional[bool] = True,
        edge_tick_threshold: Optional[int] = None,
        include_legend: Optional[bool] = False,
    ) -> Tuple:
        """
        Export the network in the selected format.

        Parameters
        ----------
        n_clicks : int
            Button click count from Dash.
        network_data : Dict
            Serialized network from the network store.
        export_format : str
            Selected export format.
        layout_data : Dict, optional
            Computed node positions.
        dragged_positions : Dict, optional
            Manually dragged node positions, which win over the layout.
        metadata_data : Dict, optional
            Serialized metadata, used for population colours.
        h_number_mapping : Dict, optional
            Custom haplotype labels, if uploaded.
        show_edge_ticks : bool, optional
            Whether mutation counts render as tick marks.
        edge_tick_threshold : int, optional
            Mutation count above which an edge shows a numeral.
        include_legend : bool, optional
            Whether to draw a population key on the figure.

        Returns
        -------
        tuple
            Download payload for the text formats, the Cytoscape image
            request for PNG/SVG, and a toast reporting the outcome.
        """
        if not network_data:
            raise PreventUpdate

        try:
            network = HaplotypeNetwork.from_serialized(network_data)

            if export_format in TEXT_EXPORTERS:
                payload = _export_text(network, export_format)
                return (
                    payload,
                    dash.no_update,
                    *toast(f'Saved {payload["filename"]}.', header='Network exported'),
                )

            if export_format in FIGURE_FORMATS:
                payload = _export_figure(
                    network,
                    merge_node_positions(layout_data, dragged_positions),
                    (metadata_data or {}).get('population_colors'),
                    export_format,
                    node_labels=h_number_mapping or _default_h_numbers(network),
                    population_mapping=(metadata_data or {}).get('populations'),
                    show_edge_ticks=bool(show_edge_ticks),
                    edge_tick_threshold=edge_tick_threshold or DEFAULT_TICK_THRESHOLD,
                    include_legend=bool(include_legend),
                )
                return (
                    payload,
                    dash.no_update,
                    *toast(f'Saved {payload["filename"]}.', header='Figure exported'),
                )

            if export_format in CANVAS_FORMATS:
                # Rendered by dash-cytoscape in the browser from what is
                # currently on screen, so it never reaches Python.
                image_config = {
                    'type': export_format,
                    'action': 'download',
                    # No extension: dash-cytoscape appends the type itself,
                    # so 'network.png' arrived as 'network.png.png'.
                    'filename': 'network',
                    'options': {
                        'output': 'base64uri',
                        'bg': 'white',
                        'full': True,
                    },
                }
                return (
                    dash.no_update,
                    image_config,
                    *toast(f'Saved network.{export_format}.', header='Image exported'),
                )

            # Without this the function fell off the end returning None,
            # and Dash errored trying to unpack it into the outputs.
            raise ValueError(f'Unsupported export format: {export_format!r}')

        except Exception as e:
            logger.error(f'Error exporting as {export_format!r}: {e}')
            logger.error(traceback.format_exc())
            # Surfaced rather than swallowed: a bare PreventUpdate here
            # made a failed export look like a dead button.
            return (
                dash.no_update,
                dash.no_update,
                *toast(str(e), header='Export failed'),
            )

    # New callbacks for enhanced features

    @app.callback(
        Output('download-haplotype-csv', 'data'),
        Input('download-haplotype-csv-button', 'n_clicks'),
        [
            State('network-store', 'data'),
            State('alignment-store', 'data'),
            State('metadata-store', 'data'),
            State('h-number-mapping-store', 'data'),
        ],
        prevent_initial_call=True,
    )
    def download_haplotype_csv(
        n_clicks: Optional[int],
        network_data: Optional[Dict],
        alignment_data: Optional[Dict],
        metadata_data: Optional[Dict],
        h_number_mapping: Optional[Dict],
    ) -> Optional[Dict]:
        """
        Download haplotype summary as CSV.

        Parameters
        ----------
        n_clicks : int, optional
            Button click count from Dash.
        network_data : Dict, optional
            Serialized network from the network store.
        alignment_data : Dict, optional
            Serialized alignment from the alignment store.
        metadata_data : Dict, optional
            Serialized metadata from the metadata store.
        h_number_mapping : Dict, optional
            Custom haplotype label mapping, if uploaded.

        Returns
        -------
        Optional[Dict]
            Download payload for the haplotype summary.
        """
        if not network_data or not alignment_data:
            raise PreventUpdate

        try:
            import csv
            import io

            # Reconstruct network
            network = HaplotypeNetwork.from_serialized(network_data)

            # Build CSV content
            output = io.StringIO()
            writer = csv.writer(output)

            # Write header
            headers = ['H_Number', 'Type', 'Frequency', 'Sample_IDs']
            if metadata_data and metadata_data.get('populations'):
                headers.append('Populations')
            writer.writerow(headers)

            # Write data rows
            for i, node_id in enumerate(sorted(network.graph.nodes()), start=1):
                node_data = network.graph.nodes[node_id]
                sample_ids = node_data.get('sample_ids', [])
                is_median = node_data.get('median_vector', False)
                frequency = node_data.get('frequency', len(sample_ids))

                # Use custom label if available, otherwise default H number
                if h_number_mapping and node_id in h_number_mapping:
                    h_label = h_number_mapping[node_id]
                else:
                    h_label = f'H{i}'

                if is_median or len(sample_ids) == 0:
                    haplotype_type = 'Inferred'
                    sample_display = 'None'
                else:
                    haplotype_type = 'Observed'
                    sample_display = '; '.join(sample_ids) if sample_ids else 'Unknown'

                row = [h_label, haplotype_type, frequency, sample_display]

                # Add populations if metadata available
                if metadata_data and metadata_data.get('populations'):
                    populations = set()
                    for sid in sample_ids:
                        if sid in metadata_data['populations']:
                            populations.add(metadata_data['populations'][sid])
                    pop_display = '; '.join(sorted(populations)) if populations else ''
                    row.append(pop_display)

                writer.writerow(row)

            return {
                'content': output.getvalue(),
                'filename': 'haplotype_summary.csv',
            }

        except Exception as e:
            logger.error(f'Error generating CSV: {e}')
            raise PreventUpdate from e
