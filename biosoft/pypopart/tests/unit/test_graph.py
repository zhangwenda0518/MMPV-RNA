"""Unit tests for HaplotypeNetwork class."""

import networkx as nx
import pytest

from pypopart.core.graph import HaplotypeNetwork
from pypopart.core.haplotype import Haplotype
from pypopart.core.sequence import Sequence


class TestHaplotypeNetwork:
    """Test cases for HaplotypeNetwork class."""

    def test_empty_network_creation(self):
        """Test creating an empty network."""
        network = HaplotypeNetwork()
        assert network.num_nodes == 0
        assert network.num_edges == 0
        assert network.name == 'HaplotypeNetwork'
        assert len(network) == 0

    def test_network_with_name(self):
        """Test creating a network with a custom name."""
        network = HaplotypeNetwork(name='MyNetwork')
        assert network.name == 'MyNetwork'

    def test_add_haplotype(self):
        """Test adding a haplotype to the network."""
        network = HaplotypeNetwork()
        seq = Sequence('H1', 'ATCG')
        hap = Haplotype(seq, sample_ids=['s1', 's2'])

        network.add_haplotype(hap)

        assert network.num_nodes == 1
        assert network.has_node('H1')
        assert network.get_haplotype('H1') == hap

    def test_add_haplotype_duplicate_id(self):
        """Test adding haplotype with duplicate ID fails."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        hap1 = Haplotype(seq1, sample_ids=['s1'])
        seq2 = Sequence('H1', 'GTCA')
        hap2 = Haplotype(seq2, sample_ids=['s2'])

        network.add_haplotype(hap1)

        with pytest.raises(ValueError, match='already exists'):
            network.add_haplotype(hap2)

    def test_add_median_vector(self):
        """Test adding a median vector node."""
        network = HaplotypeNetwork()
        seq = Sequence('MV1', 'ATCG')
        hap = Haplotype(seq)

        network.add_haplotype(hap, median_vector=True)

        assert network.has_node('MV1')
        assert network.is_median_vector('MV1')
        assert 'MV1' in network.median_vector_ids

    def test_remove_haplotype(self):
        """Test removing a haplotype from network."""
        network = HaplotypeNetwork()
        seq = Sequence('H1', 'ATCG')
        hap = Haplotype(seq, sample_ids=['s1'])

        network.add_haplotype(hap)
        assert network.has_node('H1')

        network.remove_haplotype('H1')
        assert not network.has_node('H1')
        assert network.num_nodes == 0

    def test_remove_nonexistent_haplotype(self):
        """Test removing non-existent haplotype fails."""
        network = HaplotypeNetwork()

        with pytest.raises(KeyError, match='not found'):
            network.remove_haplotype('H1')

    def test_add_edge(self):
        """Test adding an edge between haplotypes."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_edge('H1', 'H2', distance=1.0)

        assert network.num_edges == 1
        assert network.has_edge('H1', 'H2')
        assert network.get_edge_distance('H1', 'H2') == 1.0

    def test_add_edge_with_attributes(self):
        """Test adding edge with custom attributes."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_edge('H1', 'H2', distance=2.0, mutation='G>C')

        assert network.has_edge('H1', 'H2')
        assert network.get_edge_distance('H1', 'H2') == 2.0

    def test_add_edge_nonexistent_source(self):
        """Test adding edge with non-existent source fails."""
        network = HaplotypeNetwork()
        seq = Sequence('H2', 'ATCC')
        hap = Haplotype(seq)
        network.add_haplotype(hap)

        with pytest.raises(KeyError, match='Source node'):
            network.add_edge('H1', 'H2', distance=1.0)

    def test_add_edge_nonexistent_target(self):
        """Test adding edge with non-existent target fails."""
        network = HaplotypeNetwork()
        seq = Sequence('H1', 'ATCG')
        hap = Haplotype(seq)
        network.add_haplotype(hap)

        with pytest.raises(KeyError, match='Target node'):
            network.add_edge('H1', 'H2', distance=1.0)

    def test_remove_edge(self):
        """Test removing an edge."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_edge('H1', 'H2', distance=1.0)

        assert network.has_edge('H1', 'H2')
        network.remove_edge('H1', 'H2')
        assert not network.has_edge('H1', 'H2')

    def test_remove_nonexistent_edge(self):
        """Test removing non-existent edge fails."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        with pytest.raises(KeyError, match='not found'):
            network.remove_edge('H1', 'H2')

    def test_get_haplotype(self):
        """Test retrieving haplotype by ID."""
        network = HaplotypeNetwork()
        seq = Sequence('H1', 'ATCG')
        hap = Haplotype(seq, sample_ids=['s1'])

        network.add_haplotype(hap)
        retrieved = network.get_haplotype('H1')

        assert retrieved == hap
        assert retrieved.id == 'H1'

    def test_get_nonexistent_haplotype(self):
        """Test retrieving non-existent haplotype fails."""
        network = HaplotypeNetwork()

        with pytest.raises(KeyError, match='not found'):
            network.get_haplotype('H1')

    def test_get_neighbors(self):
        """Test getting neighboring nodes."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        seq3 = Sequence('H3', 'GTCG')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)
        hap3 = Haplotype(seq3)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_haplotype(hap3)
        network.add_edge('H1', 'H2', distance=1.0)
        network.add_edge('H1', 'H3', distance=1.0)

        neighbors = network.get_neighbors('H1')
        assert set(neighbors) == {'H2', 'H3'}

    def test_get_degree(self):
        """Test getting node degree."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        seq3 = Sequence('H3', 'GTCG')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)
        hap3 = Haplotype(seq3)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_haplotype(hap3)
        network.add_edge('H1', 'H2', distance=1.0)
        network.add_edge('H1', 'H3', distance=1.0)

        assert network.get_degree('H1') == 2
        assert network.get_degree('H2') == 1
        assert network.get_degree('H3') == 1

    def test_nodes_property(self):
        """Test getting list of nodes."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        nodes = network.nodes
        assert set(nodes) == {'H1', 'H2'}

    def test_edges_property(self):
        """Test getting list of edges."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_edge('H1', 'H2', distance=1.0)

        edges = network.edges
        assert len(edges) == 1
        assert ('H1', 'H2') in edges or ('H2', 'H1') in edges

    def test_haplotypes_property(self):
        """Test getting list of haplotypes excluding median vectors."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('MV1', 'ATCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)

        network.add_haplotype(hap1, median_vector=False)
        network.add_haplotype(hap2, median_vector=True)

        haplotypes = network.haplotypes
        assert len(haplotypes) == 1
        assert haplotypes[0].id == 'H1'

    def test_is_connected_true(self):
        """Test connected network detection."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_edge('H1', 'H2', distance=1.0)

        assert network.is_connected()

    def test_is_connected_false(self):
        """Test disconnected network detection."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        # No edge - disconnected

        assert not network.is_connected()

    def test_get_connected_components(self):
        """Test getting connected components."""
        network = HaplotypeNetwork()
        # Create two disconnected components
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        seq3 = Sequence('H3', 'GTCG')
        seq4 = Sequence('H4', 'GTCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)
        hap3 = Haplotype(seq3)
        hap4 = Haplotype(seq4)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_haplotype(hap3)
        network.add_haplotype(hap4)
        network.add_edge('H1', 'H2', distance=1.0)
        network.add_edge('H3', 'H4', distance=1.0)

        components = network.get_connected_components()
        assert len(components) == 2
        assert {'H1', 'H2'} in components
        assert {'H3', 'H4'} in components

    def test_calculate_diameter_connected(self):
        """Test diameter calculation for connected network."""
        network = HaplotypeNetwork()
        # Create linear network: H1 -- H2 -- H3
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        seq3 = Sequence('H3', 'GTCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)
        hap3 = Haplotype(seq3)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_haplotype(hap3)
        network.add_edge('H1', 'H2', distance=1.0)
        network.add_edge('H2', 'H3', distance=1.0)

        assert network.calculate_diameter() == 2

    def test_calculate_diameter_disconnected(self):
        """Test diameter calculation for disconnected network."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        assert network.calculate_diameter() == -1

    def test_get_shortest_path(self):
        """Test finding shortest path."""
        network = HaplotypeNetwork()
        # Create: H1 -- H2 -- H3
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        seq3 = Sequence('H3', 'GTCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)
        hap3 = Haplotype(seq3)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_haplotype(hap3)
        network.add_edge('H1', 'H2', distance=1.0)
        network.add_edge('H2', 'H3', distance=1.0)

        path = network.get_shortest_path('H1', 'H3')
        assert path == ['H1', 'H2', 'H3']

    def test_get_shortest_path_length(self):
        """Test getting shortest path length."""
        network = HaplotypeNetwork()
        # Create: H1 -- H2 -- H3
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        seq3 = Sequence('H3', 'GTCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)
        hap3 = Haplotype(seq3)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_haplotype(hap3)
        network.add_edge('H1', 'H2', distance=1.0)
        network.add_edge('H2', 'H3', distance=1.0)

        length = network.get_shortest_path_length('H1', 'H3')
        assert length == 2

    def test_calculate_centrality(self):
        """Test centrality calculation."""
        network = HaplotypeNetwork()
        # Create: H1 -- H2 -- H3, where H2 is central
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        seq3 = Sequence('H3', 'GTCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)
        hap3 = Haplotype(seq3)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_haplotype(hap3)
        network.add_edge('H1', 'H2', distance=1.0)
        network.add_edge('H2', 'H3', distance=1.0)

        centrality = network.calculate_centrality()
        assert centrality['H2'] > centrality['H1']
        assert centrality['H2'] > centrality['H3']

    def test_get_total_samples(self):
        """Test calculating total samples."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        hap1 = Haplotype(seq1, sample_ids=['s1', 's2', 's3'])
        hap2 = Haplotype(seq2, sample_ids=['s4', 's5'])

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        assert network.get_total_samples() == 5

    def test_calculate_stats(self):
        """Test comprehensive statistics calculation."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        seq3 = Sequence('MV1', 'ATCT')
        hap1 = Haplotype(seq1, sample_ids=['s1', 's2'])
        hap2 = Haplotype(seq2, sample_ids=['s3'])
        hap3 = Haplotype(seq3)

        network.add_haplotype(hap1, median_vector=False)
        network.add_haplotype(hap2, median_vector=False)
        network.add_haplotype(hap3, median_vector=True)
        network.add_edge('H1', 'MV1', distance=1.0)
        network.add_edge('MV1', 'H2', distance=1.0)

        stats = network.calculate_stats()

        assert stats.num_nodes == 3
        assert stats.num_edges == 2
        assert stats.num_haplotypes == 2
        assert stats.num_median_vectors == 1
        assert stats.total_samples == 3
        assert stats.diameter == 2
        assert stats.avg_degree == pytest.approx(4 / 3)
        assert stats.num_components == 1

    def test_validate_valid_network(self):
        """Test validation of valid network."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_edge('H1', 'H2', distance=1.0)

        network.validate()  # Should not raise

    def test_validate_negative_distance(self):
        """Test validation catches negative distances."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_edge('H1', 'H2', distance=-1.0)

        with pytest.raises(ValueError, match='negative distance'):
            network.validate()

    def test_to_networkx(self):
        """Test converting to NetworkX graph."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_edge('H1', 'H2', distance=1.0)

        nx_graph = network.to_networkx()

        assert isinstance(nx_graph, nx.Graph)
        assert nx_graph.number_of_nodes() == 2
        assert nx_graph.number_of_edges() == 1

    def test_to_dict(self):
        """Test converting network to dictionary."""
        network = HaplotypeNetwork(name='TestNet')
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        hap1 = Haplotype(seq1, sample_ids=['s1'])
        hap2 = Haplotype(seq2, sample_ids=['s2'])

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_edge('H1', 'H2', distance=1.0)

        net_dict = network.to_dict()

        assert net_dict['name'] == 'TestNet'
        assert len(net_dict['nodes']) == 2
        assert len(net_dict['edges']) == 1
        assert net_dict['edges'][0]['distance'] == 1.0

    def test_string_representation(self):
        """Test string representation."""
        network = HaplotypeNetwork(name='MyNetwork')
        seq1 = Sequence('H1', 'ATCG')
        hap1 = Haplotype(seq1)

        network.add_haplotype(hap1)

        str_repr = str(network)
        assert 'MyNetwork' in str_repr
        assert '1 haplotype' in str_repr

    def test_repr(self):
        """Test detailed representation."""
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCC')
        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)

        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        repr_str = repr(network)
        assert 'HaplotypeNetwork' in repr_str
        assert 'nodes=2' in repr_str

    def test_from_serialized(self):
        """Test reconstructing network from serialized data."""
        # Create original network
        network = HaplotypeNetwork()
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'GTCA')
        hap1 = Haplotype(seq1, sample_ids=['s1', 's2'])
        hap2 = Haplotype(seq2, sample_ids=['s3'])

        network.add_haplotype(hap1)
        network.add_haplotype(hap2, median_vector=False)
        network.add_edge('H1', 'H2', distance=2)

        # Serialize to dict (mimicking GUI storage format)
        network_data = {
            'nodes': [
                {
                    'id': 'H1',
                    'sequence': 'ATCG',
                    'frequency': 2,
                    'is_median': False,
                    'sample_ids': ['s1', 's2'],
                },
                {
                    'id': 'H2',
                    'sequence': 'GTCA',
                    'frequency': 1,
                    'is_median': False,
                    'sample_ids': ['s3'],
                },
            ],
            'edges': [{'source': 'H1', 'target': 'H2', 'weight': 2}],
        }

        # Reconstruct network
        reconstructed = HaplotypeNetwork.from_serialized(network_data)

        # Verify structure
        assert reconstructed.num_nodes == 2
        assert reconstructed.num_edges == 1
        assert reconstructed.has_node('H1')
        assert reconstructed.has_node('H2')

        # Verify haplotypes are accessible
        hap1_reconstructed = reconstructed.get_haplotype('H1')
        hap2_reconstructed = reconstructed.get_haplotype('H2')
        assert hap1_reconstructed.id == 'H1'
        assert hap1_reconstructed.data == 'ATCG'
        assert hap2_reconstructed.id == 'H2'
        assert hap2_reconstructed.data == 'GTCA'

        # Verify edge
        assert reconstructed.has_edge('H1', 'H2')
        assert reconstructed.get_edge_distance('H1', 'H2') == 2

    def test_from_serialized_with_samples_key(self):
        """Test reconstructing network using legacy 'samples' key."""
        # Some older code might use 'samples' instead of 'sample_ids'
        network_data = {
            'nodes': [
                {
                    'id': 'H1',
                    'sequence': 'ATCG',
                    'frequency': 2,
                    'is_median': False,
                    'samples': ['s1', 's2'],
                }
            ],
            'edges': [],
        }

        # Should still work with 'samples' key
        reconstructed = HaplotypeNetwork.from_serialized(network_data)
        assert reconstructed.num_nodes == 1
        hap = reconstructed.get_haplotype('H1')
        assert hap.frequency == 2
