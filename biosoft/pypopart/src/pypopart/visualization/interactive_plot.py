"""
Interactive network visualization using Plotly for PyPopART.

Provides Plotly-based interactive plotting for haplotype networks
with hover information, zoom/pan, and clickable elements.
"""

from typing import Any, Dict, List, Optional, Tuple

import networkx as nx
import numpy as np
import plotly.graph_objects as go
from plotly.graph_objects import Figure

from ..core.graph import HaplotypeNetwork


class InteractiveNetworkPlotter:
    """
    Interactive network plotter using Plotly.

    Creates interactive visualizations of haplotype networks with
    hover information, zoom/pan controls, and clickable nodes.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object to visualize.
    """

    def __init__(self, network: HaplotypeNetwork):
        """
        Initialize interactive plotter with a haplotype network.

        Parameters
        ----------
        network : HaplotypeNetwork
            HaplotypeNetwork object to visualize.
        """
        self.network = network
        self.figure = None

    def plot(
        self,
        layout: Optional[Dict[str, Tuple[float, float]]] = None,
        layout_algorithm: str = 'spring',
        node_size_scale: float = 20.0,
        node_color_map: Optional[Dict[str, str]] = None,
        population_colors: Optional[Dict[str, str]] = None,
        edge_width_scale: float = 2.0,
        show_labels: bool = True,
        show_edge_labels: bool = True,
        median_vector_color: str = 'lightgray',
        title: Optional[str] = None,
        width: int = 1000,
        height: int = 800,
        **kwargs,
    ) -> Figure:
        """
        Create an interactive network plot.

        Parameters
        ----------
        layout : Dict[str, Tuple[float, float]], optional
            Pre-computed node positions {node_id: (x, y)}.
        layout_algorithm : str, default='spring'
            NetworkX layout algorithm ('spring', 'circular', 'kamada_kawai').
        node_size_scale : float, default=20.0
            Scaling factor for node sizes.
        node_color_map : Dict[str, str], optional
            Custom color mapping {node_id: color}.
        population_colors : Dict[str, str], optional
            Color mapping for populations {pop_name: color}.
        edge_width_scale : float, default=2.0
            Scaling factor for edge widths.
        show_labels : bool, default=True
            Whether to show node labels.
        show_edge_labels : bool, default=True
            Whether to show edge labels with mutation counts.
        median_vector_color : str, default='lightgray'
            Color for median vector nodes.
        title : str, optional
            Plot title.
        width : int, default=1000
            Figure width in pixels.
        height : int, default=800
            Figure height in pixels.
        **kwargs : dict
            Additional layout arguments.

        Returns
        -------
        Figure
            Plotly Figure object.
        """
        # Get graph and compute layout if not provided
        graph = self.network._graph
        if layout is None:
            layout = self._compute_layout(graph, layout_algorithm)

        # Create figure
        self.figure = go.Figure()

        # Add edges first (so they appear below nodes)
        self._add_edges(graph, layout, edge_width_scale, show_edge_labels)

        # Add nodes
        self._add_nodes(
            graph,
            layout,
            node_size_scale,
            node_color_map,
            population_colors,
            median_vector_color,
            show_labels,
        )

        # Update layout
        plot_title = title if title else self.network.name
        self.figure.update_layout(
            title={
                'text': plot_title,
                'x': 0.5,
                'xanchor': 'center',
                'font': {'size': 20},
            },
            showlegend=True,
            hovermode='closest',
            width=width,
            height=height,
            xaxis={'showgrid': False, 'zeroline': False, 'showticklabels': False},
            yaxis={'showgrid': False, 'zeroline': False, 'showticklabels': False},
            plot_bgcolor='white',
            **kwargs,
        )

        return self.figure

    def add_population_legend(self, population_colors: Dict[str, str]) -> None:
        """
        Add a legend for population colors.

        Parameters
        ----------
        population_colors : Dict[str, str]
            Color mapping for populations.
        """
        if self.figure is None:
            raise ValueError('No plot exists. Call plot() first.')

        # Add invisible traces for legend
        for pop_name, color in sorted(population_colors.items()):
            self.figure.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode='markers',
                    marker={'size': 10, 'color': color},
                    showlegend=True,
                    name=pop_name,
                    hoverinfo='skip',
                )
            )

    def save_html(self, filename: str, auto_open: bool = False, **kwargs) -> None:
        """
        Save the interactive plot as an HTML file.

        Parameters
        ----------
        filename : str
            Output filename (should end with .html).
        auto_open : bool, default=False
            Whether to automatically open in browser.
        **kwargs : dict
            Additional arguments passed to write_html().
        """
        if self.figure is None:
            raise ValueError('No plot exists. Call plot() first.')

        self.figure.write_html(filename, auto_open=auto_open, **kwargs)

    def show(self) -> None:
        """Display the interactive plot in a browser or notebook."""
        if self.figure is None:
            raise ValueError('No plot exists. Call plot() first.')

        self.figure.show()

    def _compute_layout(
        self, graph: nx.Graph, algorithm: str
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compute node layout via the shared LayoutManager.

        Parameters
        ----------
        graph : nx.Graph
            NetworkX graph (unused; the manager works on self.network).
        algorithm : str
            Layout algorithm name.

        Returns
        -------
        dict
            Dictionary mapping node IDs to (x, y) positions.
        """
        from ..layout.algorithms import LayoutManager

        return LayoutManager(self.network).compute_layout(algorithm)

    def _add_edges(
        self,
        graph: nx.Graph,
        layout: Dict[str, Tuple[float, float]],
        width_scale: float,
        show_edge_labels: bool = True,
    ) -> None:
        """
        Add edges to the plot.

        Parameters
        ----------
        graph : nx.Graph
            NetworkX graph.
        layout : Dict[str, Tuple[float, float]]
            Node positions.
        width_scale : float
            Edge width scaling factor.
        show_edge_labels : bool, default=True
            Whether to show edge labels with mutation counts.
        """
        edge_traces = []
        edge_annotations = []

        for u, v in graph.edges():
            x0, y0 = layout[u]
            x1, y1 = layout[v]

            # Get edge weight
            weight = graph[u][v].get('weight', 1)

            # Calculate width (inverse of distance)
            width = width_scale * max(0.5, 3.0 / max(weight, 1))

            # Create hover text
            hover_text = f'{u} ↔ {v}<br>Distance: {weight} mutation(s)<br>'

            # Create edge trace
            edge_trace = go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                mode='lines',
                line={'width': width, 'color': 'rgba(150, 150, 150, 0.6)'},
                hoverinfo='text',
                hovertext=hover_text,
                showlegend=False,
            )

            edge_traces.append(edge_trace)

            # Add edge label annotation if requested
            if show_edge_labels and weight > 0:
                # Position label at midpoint of edge
                x_mid = (x0 + x1) / 2
                y_mid = (y0 + y1) / 2

                edge_annotations.append(
                    {
                        'x': x_mid,
                        'y': y_mid,
                        'text': str(int(weight)),
                        'showarrow': False,
                        'font': {'size': 10, 'color': 'rgba(100, 100, 100, 0.8)'},
                        'bgcolor': 'rgba(255, 255, 255, 0.7)',
                        'borderpad': 2,
                    }
                )

        # Add all edge traces
        for trace in edge_traces:
            self.figure.add_trace(trace)

        # Add edge label annotations
        if edge_annotations:
            self.figure.update_layout(annotations=edge_annotations)

    def _add_nodes(
        self,
        graph: nx.Graph,
        layout: Dict[str, Tuple[float, float]],
        size_scale: float,
        node_color_map: Optional[Dict[str, str]],
        population_colors: Optional[Dict[str, str]],
        median_vector_color: str,
        show_labels: bool,
    ) -> None:
        """
        Add nodes to the plot.

        Parameters
        ----------
        graph : nx.Graph
            NetworkX graph.
        layout : Dict[str, Tuple[float, float]]
            Node positions.
        size_scale : float
            Node size scaling factor.
        node_color_map : Dict[str, str], optional
            Custom node color mapping.
        population_colors : Dict[str, str], optional
            Population color mapping.
        median_vector_color : str
            Color for median vectors.
        show_labels : bool
            Whether to show node labels.
        """
        # Separate haplotypes and median vectors
        haplotype_nodes = [
            n for n in graph.nodes() if not self.network.is_median_vector(n)
        ]
        median_nodes = self.network.median_vector_ids

        # Add haplotype nodes
        if haplotype_nodes:
            self._add_node_trace(
                haplotype_nodes,
                layout,
                size_scale,
                node_color_map,
                population_colors,
                'circle',
                show_labels,
                is_median=False,
            )

        # Add median vector nodes
        if median_nodes:
            self._add_node_trace(
                median_nodes,
                layout,
                size_scale,
                dict.fromkeys(median_nodes, median_vector_color),
                None,
                'square',
                show_labels,
                is_median=True,
            )

    def _add_node_trace(
        self,
        nodes: List[str],
        layout: Dict[str, Tuple[float, float]],
        size_scale: float,
        node_color_map: Optional[Dict[str, str]],
        population_colors: Optional[Dict[str, str]],
        symbol: str,
        show_labels: bool,
        is_median: bool,
    ) -> None:
        """
        Add a trace for a group of nodes.

        Parameters
        ----------
        nodes : List[str]
            List of node IDs.
        layout : Dict[str, Tuple[float, float]]
            Node positions.
        size_scale : float
            Size scaling factor.
        node_color_map : Dict[str, str], optional
            Custom color mapping.
        population_colors : Dict[str, str], optional
            Population color mapping.
        symbol : str
            Marker symbol ('circle' or 'square').
        show_labels : bool
            Whether to show labels.
        is_median : bool
            Whether these are median vectors.
        """
        x_coords = []
        y_coords = []
        sizes = []
        colors = []
        texts = []
        hover_texts = []

        for node in nodes:
            # Position
            x, y = layout[node]
            x_coords.append(x)
            y_coords.append(y)

            # Get haplotype data
            hap = self.network.get_haplotype(node)

            # Size
            if is_median:
                sizes.append(size_scale * 0.8)
            elif hap:
                sizes.append(size_scale * np.sqrt(hap.frequency))
            else:
                sizes.append(size_scale * 0.5)

            # Color
            color = self._get_node_color(
                node, hap, node_color_map, population_colors, is_median
            )
            colors.append(color)

            # Label
            if show_labels:
                texts.append(node)
            else:
                texts.append('')

            # Hover text
            hover_text = self._create_hover_text(node, hap, is_median)
            hover_texts.append(hover_text)

        # Create node trace
        node_trace = go.Scatter(
            x=x_coords,
            y=y_coords,
            mode='markers+text' if show_labels else 'markers',
            marker={
                'size': sizes,
                'color': colors,
                'symbol': symbol,
                'line': {'width': 2, 'color': 'black'},
            },
            text=texts,
            textposition='top center',
            textfont={'size': 10, 'family': 'Arial Black'},
            hoverinfo='text',
            hovertext=hover_texts,
            showlegend=False,
        )

        self.figure.add_trace(node_trace)

    def _get_node_color(
        self,
        node: str,
        hap: Any,
        node_color_map: Optional[Dict[str, str]],
        population_colors: Optional[Dict[str, str]],
        is_median: bool,
    ) -> str:
        """
        Determine node color based on priority.

        Parameters
        ----------
        node : str
            Node ID.
        hap : Any
            Haplotype object or None.
        node_color_map : Dict[str, str], optional
            Custom color mapping.
        population_colors : Dict[str, str], optional
            Population color mapping.
        is_median : bool
            Whether this is a median vector.

        Returns
        -------
        str
            Color string.
        """
        from .style import node_color

        return node_color(
            node,
            hap,
            is_median,
            node_color_map=node_color_map,
            population_colors=population_colors,
        )

    def _create_hover_text(self, node: str, hap: Any, is_median: bool) -> str:
        """
        Create hover text for a node.

        Parameters
        ----------
        node : str
            Node ID.
        hap : Any
            Haplotype object or None.
        is_median : bool
            Whether this is a median vector.

        Returns
        -------
        str
            Formatted hover text string.
        """
        lines = [f'<b>{node}</b>']

        if is_median:
            lines.append('Type: Median Vector')
        elif hap:
            lines.append(f'Frequency: {hap.frequency}')

            # Add population information
            pop_counts = hap.get_frequency_by_population()
            if pop_counts:
                lines.append('<br><b>Populations:</b>')
                for pop, count in sorted(pop_counts.items()):
                    lines.append(f'  {pop}: {count}')

            # Add sample IDs if not too many
            sample_ids = hap.sample_ids
            if len(sample_ids) <= 10:
                lines.append('<br><b>Samples:</b>')
                lines.append(', '.join(sample_ids))
            else:
                lines.append(f'<br><b>Samples:</b> {len(sample_ids)} total')

        return '<br>'.join(lines)


def plot_interactive_network(network: HaplotypeNetwork, **kwargs) -> Figure:
    """
    Plot an interactive haplotype network.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object to visualize.
    **kwargs : dict
        Arguments passed to InteractiveNetworkPlotter.plot().

    Returns
    -------
    Figure
        Plotly Figure object.

    Examples
    --------
    >>> from pypopart.core.graph import HaplotypeNetwork
    >>> from pypopart.visualization.interactive_plot import plot_interactive_network
    >>> network = HaplotypeNetwork()
    >>> # ... build network ...
    >>> fig = plot_interactive_network(network, layout_algorithm='spring')
    >>> fig.show()
    """
    plotter = InteractiveNetworkPlotter(network)
    return plotter.plot(**kwargs)


def create_interactive_figure(
    network: HaplotypeNetwork,
    population_colors: Optional[Dict[str, str]] = None,
    filename: Optional[str] = None,
    auto_open: bool = False,
    **kwargs,
) -> Figure:
    """
    Create an interactive figure with legend.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object to visualize.
    population_colors : Dict[str, str], optional
        Color mapping for populations.
    filename : str, optional
        Optional filename to save HTML file.
    auto_open : bool, default=False
        Whether to open the file in browser.
    **kwargs : dict
        Additional arguments passed to plot().

    Returns
    -------
    Figure
        Plotly Figure object.

    Examples
    --------
    >>> fig = create_interactive_figure(
    ...     network,
    ...     population_colors={'PopA': 'red', 'PopB': 'blue'},
    ...     filename='network.html'
    ... )
    """
    plotter = InteractiveNetworkPlotter(network)

    # Create plot
    fig = plotter.plot(population_colors=population_colors, **kwargs)

    # Add legend if population colors provided
    if population_colors:
        plotter.add_population_legend(population_colors)

    # Save if filename provided
    if filename:
        plotter.save_html(filename, auto_open=auto_open)

    return fig
