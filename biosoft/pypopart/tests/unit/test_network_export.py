"""Tests for network export/import round trips (io.network_export)."""

import json

from pypopart.io import load_network, save_network
from pypopart.io.network_export import (
    CSVExporter,
    CytoscapeExporter,
    GMLExporter,
    GraphMLExporter,
    JSONExporter,
)


class TestExporters:
    """Each exporter writes a readable file for a HaplotypeNetwork."""

    def test_graphml_round_trip(self, tiny_network, tmp_path):
        """GraphML export loads back with nodes, edges, and medians."""
        path = tmp_path / 'net.graphml'
        GraphMLExporter(path).export(tiny_network)

        loaded = load_network(str(path))
        assert loaded.num_nodes == 3
        assert loaded.num_edges == 2
        assert loaded.is_median_vector('Median_0')
        assert loaded.get_edge_distance('H1', 'Median_0') == 1

    def test_gml_export(self, tiny_network, tmp_path):
        """GML export produces a loadable file."""
        path = tmp_path / 'net.gml'
        GMLExporter(path).export(tiny_network)
        loaded = load_network(str(path))
        assert loaded.num_nodes == 3

    def test_json_round_trip(self, tiny_network, tmp_path):
        """JSON export uses the nodes/edges structure and loads back."""
        path = tmp_path / 'net.json'
        JSONExporter(path).export(tiny_network)

        data = json.loads(path.read_text())
        assert {'nodes', 'edges', 'metadata'} <= set(data)

        loaded = load_network(str(path))
        assert loaded.num_nodes == 3
        assert loaded.num_edges == 2

    def test_cytoscape_export(self, tiny_network, tmp_path):
        """Cytoscape JSON export is valid JSON with elements."""
        path = tmp_path / 'net.cyjs'
        CytoscapeExporter(path).export(tiny_network)
        data = json.loads(path.read_text())
        assert 'elements' in data

    def test_csv_exports(self, tiny_network, tmp_path):
        """Node and edge CSV exports carry one row per element."""
        node_path = tmp_path / 'nodes.csv'
        CSVExporter(node_path).export_nodes(tiny_network)
        assert len(node_path.read_text().strip().splitlines()) == 4  # header + 3

        edge_path = tmp_path / 'edges.csv'
        CSVExporter(edge_path).export_edges(tiny_network)
        assert len(edge_path.read_text().strip().splitlines()) == 3  # header + 2

    def test_save_network_accepts_raw_graph(self, tiny_network, tmp_path):
        """save_network works for both HaplotypeNetwork and nx.Graph."""
        path_a = tmp_path / 'a.graphml'
        path_b = tmp_path / 'b.graphml'
        save_network(tiny_network, path_a)
        save_network(tiny_network.graph, path_b)
        assert load_network(str(path_a)).num_nodes == 3
        assert load_network(str(path_b)).num_nodes == 3
