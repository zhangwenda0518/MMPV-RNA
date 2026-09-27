"""Unit tests for Minimum Spanning Tree (MST) algorithm."""

import numpy as np
import pytest

from pypopart.algorithms.mst import MinimumSpanningTree
from pypopart.core.alignment import Alignment
from pypopart.core.distance import DistanceMatrix
from pypopart.core.sequence import Sequence


class TestMinimumSpanningTree:
    """Test cases for MST algorithm."""

    def test_mst_initialization(self):
        """Test MST algorithm initialization."""
        mst = MinimumSpanningTree()
        assert mst.distance_method == 'hamming'
        assert mst.algorithm == 'prim'

        mst_kruskal = MinimumSpanningTree(algorithm='kruskal')
        assert mst_kruskal.algorithm == 'kruskal'

    def test_mst_invalid_algorithm(self):
        """Test MST with invalid algorithm name."""
        with pytest.raises(ValueError, match='Unknown MST algorithm'):
            MinimumSpanningTree(algorithm='invalid')

    def test_mst_empty_alignment(self):
        """Test MST with empty alignment."""
        mst = MinimumSpanningTree()
        alignment = Alignment()
        network = mst.construct_network(alignment)

        assert len(network) == 0

    def test_mst_single_sequence(self):
        """Test MST with single sequence."""
        mst = MinimumSpanningTree()
        alignment = Alignment([Sequence('seq1', 'ATCG')])
        network = mst.construct_network(alignment)

        assert len(network) == 1
        assert len(network.edges) == 0

    def test_mst_two_sequences(self):
        """Test MST with two sequences."""
        mst = MinimumSpanningTree()
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCC'),  # 1 difference
            ]
        )
        network = mst.construct_network(alignment)

        assert len(network) == 2
        assert len(network.edges) == 1
        # Check edge exists between the two haplotypes
        hap_ids = [h.id for h in network.haplotypes]
        assert network.has_edge(hap_ids[0], hap_ids[1])

    def test_mst_three_sequences_prim(self):
        """Test MST with three sequences using Prim's algorithm."""
        mst = MinimumSpanningTree(algorithm='prim')
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCC'),  # 1 diff from seq1
                Sequence('seq3', 'GTCG'),  # 1 diff from seq1, 2 diff from seq2
            ]
        )
        network = mst.construct_network(alignment)

        assert len(network) == 3
        assert len(network.edges) == 2
        assert network.is_connected()

    def test_mst_three_sequences_kruskal(self):
        """Test MST with three sequences using Kruskal's algorithm."""
        mst = MinimumSpanningTree(algorithm='kruskal')
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCC'),  # 1 diff from seq1
                Sequence('seq3', 'GTCG'),  # 1 diff from seq1, 2 diff from seq2
            ]
        )
        network = mst.construct_network(alignment)

        assert len(network) == 3
        assert len(network.edges) == 2
        assert network.is_connected()

    def test_mst_prim_kruskal_same_result(self):
        """Test that Prim and Kruskal produce same total weight."""
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCC'),
                Sequence('seq3', 'GTCG'),
                Sequence('seq4', 'GTCC'),
            ]
        )

        mst_prim = MinimumSpanningTree(algorithm='prim')
        network_prim = mst_prim.construct_network(alignment)

        mst_kruskal = MinimumSpanningTree(algorithm='kruskal')
        network_kruskal = mst_kruskal.construct_network(alignment)

        # Both should produce connected networks
        assert network_prim.is_connected()
        assert network_kruskal.is_connected()

        # Both should have n-1 edges
        assert len(network_prim.edges) == 3
        assert len(network_kruskal.edges) == 3

        # Calculate total weight (should be same for both)
        weight_prim = sum(
            network_prim.get_edge_distance(u, v) or 0 for u, v in network_prim.edges
        )
        weight_kruskal = sum(
            network_kruskal.get_edge_distance(u, v) or 0
            for u, v in network_kruskal.edges
        )
        assert weight_prim == weight_kruskal

    def test_mst_with_identical_sequences(self):
        """Test MST with identical sequences (should be single haplotype)."""
        mst = MinimumSpanningTree()
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCG'),
                Sequence('seq3', 'ATCG'),
            ]
        )
        network = mst.construct_network(alignment)

        # Should collapse to single haplotype with frequency 3
        assert len(network) == 1
        assert len(network.edges) == 0
        hap = network.haplotypes[0]
        assert hap.frequency == 3

    def test_mst_with_provided_distance_matrix(self):
        """Test MST with pre-computed distance matrix."""
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG'),
                Sequence('seq2', 'ATCC'),
                Sequence('seq3', 'GTCG'),
            ]
        )

        # Create custom distance matrix
        labels = ['seq1', 'seq2', 'seq3']
        matrix = np.array([[0, 1, 1], [1, 0, 2], [1, 2, 0]])
        dist_matrix = DistanceMatrix(labels, matrix)

        mst = MinimumSpanningTree()
        network = mst.construct_network(alignment, distance_matrix=dist_matrix)

        assert len(network) == 3
        assert len(network.edges) == 2

    def test_mst_star_topology(self):
        """Test MST creates star topology for equidistant sequences."""
        mst = MinimumSpanningTree()
        alignment = Alignment(
            [
                Sequence('seq1', 'AAAA'),  # Central
                Sequence('seq2', 'AAAT'),  # 1 diff
                Sequence('seq3', 'AAAC'),  # 1 diff
                Sequence('seq4', 'AAAG'),  # 1 diff
            ]
        )
        network = mst.construct_network(alignment)

        assert len(network) == 4
        assert len(network.edges) == 3

        # Check if one haplotype has degree 3 (central node)
        degrees = [network.get_degree(hap.id) for hap in network.haplotypes]
        assert max(degrees) == 3

    def test_mst_parameters(self):
        """Test getting MST parameters."""
        mst = MinimumSpanningTree(distance_method='k2p', algorithm='kruskal')
        params = mst.get_parameters()

        assert params['distance_method'] == 'k2p'
        assert params['algorithm'] == 'kruskal'

    def test_mst_string_representation(self):
        """Test string representation of MST."""
        mst = MinimumSpanningTree(distance_method='hamming')
        assert 'MinimumSpanningTree' in str(mst)
        assert 'hamming' in str(mst)
