"""
Network serialization for dcc.Store round trips in the GUI.

The GUI keeps the computed network in a dcc.Store as plain JSON; these
helpers are the single place that defines that wire format, shared by
the compute and rendering callbacks and by HaplotypeNetwork's
from_serialized reader.
"""

from typing import Dict, Optional, Tuple

import networkx as nx

from pypopart.core.graph import HaplotypeNetwork


def network_to_store(network: HaplotypeNetwork) -> Dict:
    """
    Serialise a network into the dcc.Store JSON format.

    Parameters
    ----------
    network : HaplotypeNetwork
        Computed network.

    Returns
    -------
    dict
        Store payload with 'nodes' and 'edges' lists, compatible with
        HaplotypeNetwork.from_serialized.
    """
    graph = network.graph
    return {
        'nodes': [
            {
                'id': node,
                'sequence': graph.nodes[node].get('sequence', ''),
                'frequency': graph.nodes[node].get('frequency', 1),
                'is_median': graph.nodes[node].get('median_vector', False),
                'sample_ids': graph.nodes[node].get('sample_ids', []),
            }
            for node in graph.nodes()
        ],
        'edges': [
            {
                'source': u,
                'target': v,
                'distance': graph[u][v].get('distance', 0),
                'weight': graph[u][v].get('weight', 1.0),
            }
            for u, v in graph.edges()
        ],
    }


def store_to_networkx(network_data: Dict) -> nx.Graph:
    """
    Rebuild a NetworkX graph from the dcc.Store payload.

    Parameters
    ----------
    network_data : dict
        Store payload from network_to_store.

    Returns
    -------
    nx.Graph
        Graph with node attributes (sequence, frequency, is_median,
        sample_ids) and edge attributes (distance, weight).
    """
    graph = nx.Graph()
    for node in network_data.get('nodes', []):
        graph.add_node(
            node['id'],
            sequence=node.get('sequence', ''),
            frequency=node.get('frequency', 1),
            is_median=node.get('is_median', False),
            median_vector=node.get('is_median', False),
            sample_ids=node.get('sample_ids', []),
        )
    for edge in network_data.get('edges', []):
        graph.add_edge(
            edge['source'],
            edge['target'],
            distance=edge.get('distance', 0),
            weight=edge.get('weight', 1.0),
        )
    return graph


def merge_node_positions(
    layout_data: Optional[Dict], dragged: Optional[Dict]
) -> Dict[str, Tuple[float, float]]:
    """
    Overlay manually dragged node positions onto a computed layout.

    Dragged positions live in their own store rather than in
    ``layout-store``. Writing them back into the layout store would make
    every drag re-trigger the graph rebuild that reads it, which
    regenerates all elements, re-fits the view and pushes the position
    through a lossy float round trip.

    Parameters
    ----------
    layout_data : dict, optional
        Positions from ``layout-store``, as ``{node: [x, y]}``.
    dragged : dict, optional
        Positions from ``node-positions-store``, same shape. These win.

    Returns
    -------
    dict
        Node positions as ``{node: (x, y)}``, ready for the plotters.
    """
    positions = {
        node: (float(pos[0]), float(pos[1]))
        for node, pos in (layout_data or {}).items()
    }
    for node, pos in (dragged or {}).items():
        # A drag for a node the current layout no longer has is stale.
        if node in positions:
            positions[node] = (float(pos[0]), float(pos[1]))
    return positions
