"""Unit tests for topology analysis module."""

from pypopart.core.graph import HaplotypeNetwork
from pypopart.core.haplotype import Haplotype
from pypopart.core.sequence import Sequence
from pypopart.stats.topology import (
    calculate_node_centrality,
    calculate_topology_summary,
    detect_bridges,
    detect_network_partitions,
    find_central_hub_nodes,
    identify_ancestral_nodes,
    identify_bottleneck_nodes,
    identify_star_patterns,
)


class TestIdentifyStarPatterns:
    """Test star pattern identification."""

    def test_empty_network(self):
        """Test with empty network."""
        network = HaplotypeNetwork()
        patterns = identify_star_patterns(network)

        assert patterns == []

    def test_perfect_star(self):
        """Test perfect star pattern."""
        center = Haplotype(Sequence('center', 'ATCG'), sample_ids=['s0'])
        leaves = [
            Haplotype(Sequence(f'leaf{i}', data), sample_ids=[f's{i}'])
            for i, data in enumerate(['ATCC', 'GTCG', 'GTCC', 'TTAA'], start=1)
        ]

        network = HaplotypeNetwork()
        network.add_haplotype(center)
        for leaf in leaves:
            network.add_haplotype(leaf)
            network.add_edge(center.id, leaf.id, distance=1)

        patterns = identify_star_patterns(network, min_leaves=3)

        assert len(patterns) == 1
        assert patterns[0].center == center.id
        assert len(patterns[0].leaves) == 4
        assert patterns[0].is_perfect_star is True

    def test_partial_star(self):
        """Test partial star pattern (center has non-leaf connections)."""
        center = Haplotype(Sequence('center', 'ATCG'), sample_ids=['s0'])
        leaf1 = Haplotype(Sequence('leaf1', 'ATCC'), sample_ids=['s1'])
        leaf2 = Haplotype(Sequence('leaf2', 'GTCG'), sample_ids=['s2'])
        leaf3 = Haplotype(Sequence('leaf3', 'GTCC'), sample_ids=['s3'])
        hub = Haplotype(Sequence('hub', 'TTAA'), sample_ids=['s4'])

        network = HaplotypeNetwork()
        for hap in [center, leaf1, leaf2, leaf3, hub]:
            network.add_haplotype(hap)

        # Center connects to leaves and hub
        network.add_edge(center.id, leaf1.id, distance=1)
        network.add_edge(center.id, leaf2.id, distance=1)
        network.add_edge(center.id, leaf3.id, distance=1)
        network.add_edge(center.id, hub.id, distance=1)

        # Hub connects to another node
        extra = Haplotype(Sequence('extra', 'GGGG'), sample_ids=['s5'])
        network.add_haplotype(extra)
        network.add_edge(hub.id, extra.id, distance=1)

        patterns = identify_star_patterns(network, min_leaves=2)

        # Should find center with 3 leaves
        assert len(patterns) >= 1
        center_pattern = [p for p in patterns if p.center == center.id][0]
        assert len(center_pattern.leaves) == 3
        assert center_pattern.is_perfect_star is False


class TestDetectNetworkPartitions:
    """Test network partition detection."""

    def test_connected_network(self):
        """Test fully connected network."""
        haps = [
            Haplotype(Sequence(f'hap{i}', data), sample_ids=[f's{i}'])
            for i, data in enumerate(['ATCG', 'ATCC', 'GTCG'])
        ]

        network = HaplotypeNetwork()
        for hap in haps:
            network.add_haplotype(hap)

        network.add_edge(haps[0].id, haps[1].id, distance=1)
        network.add_edge(haps[1].id, haps[2].id, distance=1)

        partitions = detect_network_partitions(network)

        assert len(partitions) == 1
        assert partitions[0].num_nodes == 3
        assert partitions[0].is_connected is True

    def test_disconnected_network(self):
        """Test disconnected network."""
        # Component 1
        hap1 = Haplotype(Sequence('hap1', 'ATCG'), sample_ids=['s1'])
        hap2 = Haplotype(Sequence('hap2', 'ATCC'), sample_ids=['s2'])

        # Component 2
        hap3 = Haplotype(Sequence('hap3', 'GTCG'), sample_ids=['s3'])
        hap4 = Haplotype(Sequence('hap4', 'GTCC'), sample_ids=['s4'])

        network = HaplotypeNetwork()
        for hap in [hap1, hap2, hap3, hap4]:
            network.add_haplotype(hap)

        network.add_edge(hap1.id, hap2.id, distance=1)
        network.add_edge(hap3.id, hap4.id, distance=1)

        partitions = detect_network_partitions(network)

        assert len(partitions) == 2
        assert all(p.num_nodes == 2 for p in partitions)

    def test_singleton_partition(self):
        """Test network with isolated node."""
        hap1 = Haplotype(Sequence('hap1', 'ATCG'), sample_ids=['s1'])
        hap2 = Haplotype(Sequence('hap2', 'ATCC'), sample_ids=['s2'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        # No edges - two singletons

        partitions = detect_network_partitions(network)

        assert len(partitions) == 2
        assert all(p.num_nodes == 1 for p in partitions)
        assert all(p.diameter == 0 for p in partitions)


class TestCalculateNodeCentrality:
    """Test node centrality calculations."""

    def test_empty_network(self):
        """Test with empty network."""
        network = HaplotypeNetwork()
        centrality = calculate_node_centrality(network)

        assert centrality == {}

    def test_degree_centrality(self):
        """Test degree centrality calculation."""
        center = Haplotype(Sequence('center', 'ATCG'), sample_ids=['s0'])
        leaves = [
            Haplotype(Sequence(f'leaf{i}', data), sample_ids=[f's{i}'])
            for i, data in enumerate(['ATCC', 'GTCG'], start=1)
        ]

        network = HaplotypeNetwork()
        network.add_haplotype(center)
        for leaf in leaves:
            network.add_haplotype(leaf)
            network.add_edge(center.id, leaf.id, distance=1)

        centrality = calculate_node_centrality(network, methods=['degree'])

        assert 'degree' in centrality[center.id]
        # Center should have highest degree
        assert centrality[center.id]['degree'] > centrality[leaves[0].id]['degree']

    def test_multiple_methods(self):
        """Test calculating multiple centrality methods."""
        haps = [
            Haplotype(Sequence(f'hap{i}', data), sample_ids=[f's{i}'])
            for i, data in enumerate(['ATCG', 'ATCC', 'GTCG', 'GTCC', 'TTAA'])
        ]

        network = HaplotypeNetwork()
        for hap in haps:
            network.add_haplotype(hap)

        # Create linear chain
        for i in range(len(haps) - 1):
            network.add_edge(haps[i].id, haps[i + 1].id, distance=1)

        centrality = calculate_node_centrality(network)

        assert 'degree' in centrality[haps[0].id]
        assert 'betweenness' in centrality[haps[0].id]
        assert 'closeness' in centrality[haps[0].id]


class TestIdentifyAncestralNodes:
    """Test ancestral node identification."""

    def test_empty_network(self):
        """Test with empty network."""
        network = HaplotypeNetwork()
        ancestral = identify_ancestral_nodes(network)

        assert ancestral == []

    def test_high_frequency_center(self):
        """Test that high frequency central node is identified."""
        center = Haplotype(
            Sequence('center', 'ATCG'), sample_ids=[f's{i}' for i in range(10)]
        )
        leaves = [
            Haplotype(Sequence(f'leaf{i}', data), sample_ids=[f's{i + 10}'])
            for i, data in enumerate(['ATCC', 'GTCG', 'GTCC'])
        ]

        network = HaplotypeNetwork()
        network.add_haplotype(center)
        for leaf in leaves:
            network.add_haplotype(leaf)
            network.add_edge(center.id, leaf.id, distance=1)

        ancestral = identify_ancestral_nodes(network, top_n=3)

        assert len(ancestral) >= 1
        # Center should have highest score
        assert ancestral[0].node_id == center.id
        assert ancestral[0].frequency == 10

    def test_score_components(self):
        """Test that score includes frequency, degree, and centrality."""
        hap1 = Haplotype(
            Sequence('hap1', 'ATCG'), sample_ids=[f's{i}' for i in range(5)]
        )
        hap2 = Haplotype(Sequence('hap2', 'ATCC'), sample_ids=['s5'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_edge(hap1.id, hap2.id, distance=1)

        ancestral = identify_ancestral_nodes(network)

        assert ancestral[0].degree > 0
        assert ancestral[0].frequency > 0
        assert ancestral[0].centrality >= 0
        assert ancestral[0].score > 0


class TestCalculateTopologySummary:
    """Test topology summary calculation."""

    def test_empty_network(self):
        """Test with empty network."""
        network = HaplotypeNetwork()
        summary = calculate_topology_summary(network)

        assert summary['num_nodes'] == 0
        assert summary['num_edges'] == 0
        assert summary['star_patterns'] == []
        assert summary['partitions'] == []

    def test_simple_network(self):
        """Test with simple network."""
        haps = [
            Haplotype(Sequence(f'hap{i}', data), sample_ids=[f's{i}'])
            for i, data in enumerate(['ATCG', 'ATCC', 'GTCG'])
        ]

        network = HaplotypeNetwork()
        for hap in haps:
            network.add_haplotype(hap)

        network.add_edge(haps[0].id, haps[1].id, distance=1)
        network.add_edge(haps[1].id, haps[2].id, distance=1)

        summary = calculate_topology_summary(network)

        assert summary['num_nodes'] == 3
        assert summary['num_edges'] == 2
        assert summary['is_connected'] is True
        assert summary['num_components'] == 1
        assert 'degree_distribution' in summary
        assert 'ancestral_candidates' in summary


class TestFindCentralHubNodes:
    """Test hub node identification."""

    def test_empty_network(self):
        """Test with empty network."""
        network = HaplotypeNetwork()
        hubs = find_central_hub_nodes(network)

        assert hubs == []

    def test_star_network(self):
        """Test identifying hub in star network."""
        center = Haplotype(Sequence('center', 'ATCG'), sample_ids=['s0'])
        leaves = [
            Haplotype(Sequence(f'leaf{i}', data), sample_ids=[f's{i}'])
            for i, data in enumerate(['ATCC', 'GTCG', 'GTCC'], start=1)
        ]

        network = HaplotypeNetwork()
        network.add_haplotype(center)
        for leaf in leaves:
            network.add_haplotype(leaf)
            network.add_edge(center.id, leaf.id, distance=1)

        hubs = find_central_hub_nodes(network)

        # Center should be identified as hub
        assert len(hubs) >= 1
        assert hubs[0][0] == center.id
        assert hubs[0][1] == 3  # Degree of 3

    def test_custom_threshold(self):
        """Test with custom degree threshold."""
        haps = [
            Haplotype(Sequence(f'hap{i}', data), sample_ids=[f's{i}'])
            for i, data in enumerate(['ATCG', 'ATCC', 'GTCG'])
        ]

        network = HaplotypeNetwork()
        for hap in haps:
            network.add_haplotype(hap)

        network.add_edge(haps[0].id, haps[1].id, distance=1)
        network.add_edge(haps[1].id, haps[2].id, distance=1)

        # With threshold=2, only hap1 (degree 2) should be a hub
        hubs = find_central_hub_nodes(network, degree_threshold=2)

        assert len(hubs) == 1
        assert hubs[0][0] == haps[1].id


class TestDetectBridges:
    """Test bridge edge detection."""

    def test_no_bridges(self):
        """Test network with no bridges (cycle)."""
        haps = [
            Haplotype(Sequence(f'hap{i}', data), sample_ids=[f's{i}'])
            for i, data in enumerate(['ATCG', 'ATCC', 'GTCG'])
        ]

        network = HaplotypeNetwork()
        for hap in haps:
            network.add_haplotype(hap)

        # Create triangle
        network.add_edge(haps[0].id, haps[1].id, distance=1)
        network.add_edge(haps[1].id, haps[2].id, distance=1)
        network.add_edge(haps[2].id, haps[0].id, distance=1)

        bridges = detect_bridges(network)

        assert bridges == []

    def test_with_bridges(self):
        """Test network with bridges."""
        haps = [
            Haplotype(Sequence(f'hap{i}', data), sample_ids=[f's{i}'])
            for i, data in enumerate(['ATCG', 'ATCC', 'GTCG'])
        ]

        network = HaplotypeNetwork()
        for hap in haps:
            network.add_haplotype(hap)

        # Create linear chain (all edges are bridges)
        network.add_edge(haps[0].id, haps[1].id, distance=1)
        network.add_edge(haps[1].id, haps[2].id, distance=1)

        bridges = detect_bridges(network)

        assert len(bridges) == 2


class TestIdentifyBottleneckNodes:
    """Test bottleneck node identification."""

    def test_no_bottlenecks(self):
        """Test network with no bottlenecks (triangle)."""
        haps = [
            Haplotype(Sequence(f'hap{i}', data), sample_ids=[f's{i}'])
            for i, data in enumerate(['ATCG', 'ATCC', 'GTCG'])
        ]

        network = HaplotypeNetwork()
        for hap in haps:
            network.add_haplotype(hap)

        # Create triangle (no articulation points)
        network.add_edge(haps[0].id, haps[1].id, distance=1)
        network.add_edge(haps[1].id, haps[2].id, distance=1)
        network.add_edge(haps[2].id, haps[0].id, distance=1)

        bottlenecks = identify_bottleneck_nodes(network)

        assert bottlenecks == []

    def test_with_bottleneck(self):
        """Test network with bottleneck node."""
        # Create dumbbell: hap1 - hap2 - hap3
        # where hap2 is the bottleneck
        haps = [
            Haplotype(Sequence(f'hap{i}', data), sample_ids=[f's{i}'])
            for i, data in enumerate(['ATCG', 'ATCC', 'GTCG'])
        ]

        network = HaplotypeNetwork()
        for hap in haps:
            network.add_haplotype(hap)

        network.add_edge(haps[0].id, haps[1].id, distance=1)
        network.add_edge(haps[1].id, haps[2].id, distance=1)

        bottlenecks = identify_bottleneck_nodes(network)

        # hap1 (middle node) is articulation point
        assert len(bottlenecks) == 1
        assert bottlenecks[0][0] == haps[1].id
