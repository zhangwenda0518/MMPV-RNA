"""Golden-network validation against hand-traced PopART expectations."""

from .conftest import canonical_edges, run_golden


def test_golden_network(golden):
    """The constructed network matches the golden node and edge sets."""
    network = run_golden(golden)

    expected_nodes = set(golden['expected']['nodes'])
    actual_nodes = {str(n) for n in network.graph.nodes()}
    assert actual_nodes == expected_nodes

    expected_edges = {
        (min(u, v), max(u, v), float(d)) for u, v, d in golden['expected']['edges']
    }
    assert canonical_edges(network) == expected_edges
