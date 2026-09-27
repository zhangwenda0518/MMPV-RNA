"""
Network export module for PyPopART.

Provides exporters for various network file formats.
"""

import csv
import json
from pathlib import Path
from typing import Any, Dict, Union

import networkx as nx

from pypopart.core.graph import HaplotypeNetwork
from pypopart.core.haplotype import Haplotype


def _as_nx_graph(network: Union[HaplotypeNetwork, nx.Graph]) -> nx.Graph:
    """
    Return the underlying NetworkX graph for a network object.

    A plain ``nx.Graph`` also has a ``.graph`` attribute (its attribute
    dict), so this must type-check rather than duck-type.

    Parameters
    ----------
    network : HaplotypeNetwork or nx.Graph
        HaplotypeNetwork or NetworkX graph.

    Returns
    -------
    nx.Graph
        The NetworkX graph itself, or the network's underlying graph.
    """
    if isinstance(network, nx.Graph):
        return network
    return network.graph


def _sanitize_graph_for_export(
    graph: nx.Graph, format_type: str = 'generic'
) -> nx.Graph:
    """
    Create a copy of the graph with serializable attributes.

    Converts Haplotype objects and other non-serializable types to
    basic Python types for export.

    Parameters
    ----------
    graph : nx.Graph
        NetworkX graph with potentially non-serializable attributes.
    format_type : str, default='generic'
        Export format type. 'graphml' uses stricter serialization (only primitives),
        'generic' allows lists and dicts for formats like JSON and GML.

    Returns
    -------
    nx.Graph
        New graph with sanitized attributes.
    """
    # Create a deep copy to avoid modifying original
    sanitized = nx.Graph()

    # GraphML only supports str, int, float, bool (no lists or dicts)
    strict_mode = format_type == 'graphml'

    # Copy nodes with sanitized attributes
    for node, attrs in graph.nodes(data=True):
        clean_attrs = {}
        for key, value in attrs.items():
            if key == 'haplotype':
                # Skip the Haplotype object - we have the data in other attributes
                continue
            elif isinstance(value, Haplotype):
                # Convert Haplotype to string
                clean_attrs[key] = str(value)
            elif isinstance(value, (list, tuple, set)):
                if strict_mode:
                    # Convert to comma-separated string for GraphML
                    clean_attrs[key] = ','.join(str(item) for item in value)
                else:
                    # Convert to list for formats that support it
                    clean_attrs[key] = [
                        str(item)
                        if not isinstance(item, (str, int, float, bool, type(None)))
                        else item
                        for item in value
                    ]
            elif isinstance(value, dict):
                if strict_mode:
                    # Convert to JSON string for GraphML
                    import json

                    clean_attrs[key] = json.dumps(value)
                else:
                    # Keep as dict for formats that support it
                    clean_attrs[key] = {
                        k: str(v)
                        if not isinstance(v, (str, int, float, bool, type(None)))
                        else v
                        for k, v in value.items()
                    }
            elif isinstance(value, (str, int, float, bool, type(None))):
                # Already serializable
                clean_attrs[key] = value
            else:
                # Convert anything else to string
                clean_attrs[key] = str(value)

        sanitized.add_node(node, **clean_attrs)

    # Copy edges with sanitized attributes
    for source, target, attrs in graph.edges(data=True):
        clean_attrs = {}
        for key, value in attrs.items():
            if isinstance(value, (str, int, float, bool, type(None))):
                clean_attrs[key] = value
            elif isinstance(value, (list, tuple)):
                if strict_mode:
                    clean_attrs[key] = ','.join(str(item) for item in value)
                else:
                    clean_attrs[key] = list(value)
            elif isinstance(value, dict):
                if strict_mode:
                    import json

                    clean_attrs[key] = json.dumps(value)
                else:
                    clean_attrs[key] = dict(value)
            else:
                clean_attrs[key] = str(value)

        sanitized.add_edge(source, target, **clean_attrs)

    # Copy graph-level attributes
    if hasattr(graph, 'graph'):
        sanitized.graph.update(graph.graph)

    return sanitized


class GraphMLExporter:
    """
    Export haplotype networks to GraphML format.

    Parameters
    ----------
    filepath : str or Path
        Output file path.
    """

    def __init__(self, filepath: Union[str, Path]):
        """
        Initialize GraphML exporter.

        Parameters
        ----------
        filepath : str or Path
            Output file path.
        """
        self.filepath = Path(filepath)

    def export(self, network: HaplotypeNetwork) -> None:
        """
        Export network to GraphML format.

        Parameters
        ----------
        network : HaplotypeNetwork
            HaplotypeNetwork object.
        """
        # Convert network to NetworkX graph if needed
        graph = _as_nx_graph(network)

        # Sanitize graph data for export (strict mode for GraphML)
        sanitized_graph = _sanitize_graph_for_export(graph, format_type='graphml')

        # Write to GraphML
        nx.write_graphml(sanitized_graph, self.filepath)


class GMLExporter:
    """
    Export haplotype networks to GML format.

    Parameters
    ----------
    filepath : str or Path
        Output file path.
    """

    def __init__(self, filepath: Union[str, Path]):
        """
        Initialize GML exporter.

        Parameters
        ----------
        filepath : str or Path
            Output file path.
        """
        self.filepath = Path(filepath)

    def export(self, network: HaplotypeNetwork) -> None:
        """
        Export network to GML format.

        Parameters
        ----------
        network : HaplotypeNetwork
            HaplotypeNetwork object.
        """
        graph = _as_nx_graph(network)

        # Sanitize graph data for export (GML requires strings like GraphML)
        sanitized_graph = _sanitize_graph_for_export(graph, format_type='graphml')

        # Write to GML
        nx.write_gml(sanitized_graph, self.filepath)


class CytoscapeExporter:
    """
    Export haplotype networks to Cytoscape JSON format.

    Parameters
    ----------
    filepath : str or Path
        Output file path.
    """

    def __init__(self, filepath: Union[str, Path]):
        """
        Initialize Cytoscape exporter.

        Parameters
        ----------
        filepath : str or Path
            Output file path.
        """
        self.filepath = Path(filepath)

    def export(self, network: HaplotypeNetwork) -> None:
        """
        Export network to Cytoscape JSON format.

        Parameters
        ----------
        network : HaplotypeNetwork
            HaplotypeNetwork object.
        """
        graph = _as_nx_graph(network)

        # Sanitize graph data for export
        sanitized_graph = _sanitize_graph_for_export(graph)

        # Convert to Cytoscape JSON format
        cytoscape_data = nx.cytoscape_data(sanitized_graph)

        with open(self.filepath, 'w') as f:
            json.dump(cytoscape_data, f, indent=2)


class JSONExporter:
    """
    Export haplotype networks to JSON format.

    Parameters
    ----------
    filepath : str or Path
        Output file path.
    """

    def __init__(self, filepath: Union[str, Path]):
        """
        Initialize JSON exporter.

        Parameters
        ----------
        filepath : str or Path
            Output file path.
        """
        self.filepath = Path(filepath)

    def export(self, network: HaplotypeNetwork, include_layout: bool = True) -> None:
        """
        Export network to JSON format.

        Parameters
        ----------
        network : HaplotypeNetwork or nx.Graph
            Network to export.
        include_layout : bool
            Whether to include node layout positions.
        """
        graph = _as_nx_graph(network)

        # Sanitize graph data for export
        sanitized_graph = _sanitize_graph_for_export(graph)

        # Build JSON structure
        data = {'nodes': [], 'edges': [], 'metadata': {}}

        # Add nodes
        for node, attrs in sanitized_graph.nodes(data=True):
            node_data = {'id': str(node), 'attributes': dict(attrs)}
            data['nodes'].append(node_data)

        # Add edges
        for source, target, attrs in sanitized_graph.edges(data=True):
            edge_data = {
                'source': str(source),
                'target': str(target),
                'attributes': dict(attrs),
            }
            data['edges'].append(edge_data)

        # Add graph metadata
        if hasattr(sanitized_graph, 'graph'):
            data['metadata'] = dict(sanitized_graph.graph)

        with open(self.filepath, 'w') as f:
            json.dump(data, f, indent=2)


class CSVExporter:
    """
    Export haplotype network statistics to CSV format.

    Parameters
    ----------
    filepath : str or Path
        Output file path.
    """

    def __init__(self, filepath: Union[str, Path]):
        """
        Initialize CSV exporter.

        Parameters
        ----------
        filepath : str or Path
            Output file path.
        """
        self.filepath = Path(filepath)

    def export_nodes(self, network: HaplotypeNetwork) -> None:
        """
        Export node attributes to CSV.

        Parameters
        ----------
        network : HaplotypeNetwork or nx.Graph
            Network to export.
        """
        graph = _as_nx_graph(network)

        # Collect all attribute keys
        all_keys = set()
        for _node, attrs in graph.nodes(data=True):
            all_keys.update(attrs.keys())

        fieldnames = ['node_id'] + sorted(all_keys)

        with open(self.filepath, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()

            for node, attrs in graph.nodes(data=True):
                row = {'node_id': str(node)}
                row.update(attrs)
                writer.writerow(row)

    def export_edges(self, network: HaplotypeNetwork) -> None:
        """
        Export edge attributes to CSV.

        Parameters
        ----------
        network : HaplotypeNetwork or nx.Graph
            Network to export.
        """
        graph = _as_nx_graph(network)

        # Collect all attribute keys
        all_keys = set()
        for _source, _target, attrs in graph.edges(data=True):
            all_keys.update(attrs.keys())

        fieldnames = ['source', 'target'] + sorted(all_keys)

        with open(self.filepath, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()

            for source, target, attrs in graph.edges(data=True):
                row = {'source': str(source), 'target': str(target)}
                row.update(attrs)
                writer.writerow(row)

    def export_statistics(
        self, network: HaplotypeNetwork, statistics: Dict[str, Any]
    ) -> None:
        """
        Export network statistics to CSV.

        Parameters
        ----------
        network : HaplotypeNetwork or nx.Graph
            Network the statistics belong to.
        statistics : dict
            Dictionary of statistics to export.
        """
        with open(self.filepath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Statistic', 'Value'])

            for key, value in sorted(statistics.items()):
                writer.writerow([key, value])
