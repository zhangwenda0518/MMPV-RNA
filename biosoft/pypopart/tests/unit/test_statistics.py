"""Unit tests for statistics module."""

import numpy as np
import pytest

from pypopart.core.alignment import Alignment
from pypopart.core.graph import HaplotypeNetwork
from pypopart.core.haplotype import Haplotype
from pypopart.core.sequence import Sequence
from pypopart.stats.statistics import (
    calculate_diversity_metrics,
    calculate_haplotype_frequencies,
    calculate_network_metrics,
    calculate_reticulation_index,
    calculate_summary_statistics,
    get_frequency_distribution,
    identify_central_haplotypes,
)


class TestCalculateHaplotypeFrequencies:
    """Test haplotype frequency calculations."""

    def test_empty_network(self):
        """Test with empty network."""
        network = HaplotypeNetwork()
        freq_info = calculate_haplotype_frequencies(network)

        assert freq_info['overall'] == {}
        assert freq_info['by_population'] == {}
        assert freq_info['total_samples'] == 0

    def test_single_haplotype(self):
        """Test with single haplotype."""
        seq = Sequence('seq1', 'ATCG')
        hap = Haplotype(seq, sample_ids=['s1', 's2', 's3'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap)

        freq_info = calculate_haplotype_frequencies(network)

        assert freq_info['overall'][hap.id] == 3
        assert freq_info['total_samples'] == 3

    def test_multiple_haplotypes(self):
        """Test with multiple haplotypes."""
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATCC')
        hap1 = Haplotype(seq1, sample_ids=['s1', 's2'])
        hap2 = Haplotype(seq2, sample_ids=['s3'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        freq_info = calculate_haplotype_frequencies(network)

        assert freq_info['overall'][hap1.id] == 2
        assert freq_info['overall'][hap2.id] == 1
        assert freq_info['total_samples'] == 3

    def test_with_populations(self):
        """Test frequency calculation by population."""
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATCC')

        hap1 = Haplotype(
            seq1, sample_ids=['s1', 's2'], populations={'s1': 'PopA', 's2': 'PopA'}
        )
        hap2 = Haplotype(
            seq2, sample_ids=['s3', 's4'], populations={'s3': 'PopB', 's4': 'PopB'}
        )

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        freq_info = calculate_haplotype_frequencies(network)

        assert freq_info['by_population']['PopA'][hap1.id] == 2
        assert freq_info['by_population']['PopB'][hap2.id] == 2

    def test_normalized_frequencies(self):
        """Test normalized frequency calculation."""
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATCC')
        hap1 = Haplotype(seq1, sample_ids=['s1', 's2', 's3'])
        hap2 = Haplotype(seq2, sample_ids=['s4'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        freq_info = calculate_haplotype_frequencies(network, normalize=True)

        assert freq_info['overall'][hap1.id] == pytest.approx(0.75)
        assert freq_info['overall'][hap2.id] == pytest.approx(0.25)


class TestCalculateDiversityMetrics:
    """Test diversity metrics calculation."""

    def test_empty_network(self):
        """Test with empty network."""
        network = HaplotypeNetwork()
        metrics = calculate_diversity_metrics(network)

        assert metrics.num_haplotypes == 0
        assert metrics.num_samples == 0
        assert metrics.haplotype_diversity == 0.0
        assert metrics.shannon_index == 0.0

    def test_single_haplotype(self):
        """Test with single haplotype."""
        seq = Sequence('seq1', 'ATCG')
        hap = Haplotype(seq, sample_ids=['s1', 's2'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap)

        metrics = calculate_diversity_metrics(network)

        assert metrics.num_haplotypes == 1
        assert metrics.num_samples == 2
        assert metrics.haplotype_diversity == 0.0  # All same haplotype
        assert metrics.shannon_index == 0.0

    def test_two_equal_haplotypes(self):
        """Test with two haplotypes of equal frequency."""
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATCC')
        hap1 = Haplotype(seq1, sample_ids=['s1'])
        hap2 = Haplotype(seq2, sample_ids=['s2'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        metrics = calculate_diversity_metrics(network)

        assert metrics.num_haplotypes == 2
        assert metrics.num_samples == 2
        # H = (n/(n-1)) * (1 - sum(p^2)) = (2/1) * (1 - 0.5) = 1.0
        assert metrics.haplotype_diversity == pytest.approx(1.0)
        assert metrics.shannon_index == pytest.approx(np.log(2))

    def test_with_alignment(self):
        """Test with alignment for nucleotide diversity."""
        sequences = [Sequence('seq1', 'ATCG'), Sequence('seq2', 'ATCC')]
        alignment = Alignment(sequences)

        hap1 = Haplotype(sequences[0], sample_ids=['s1'])
        hap2 = Haplotype(sequences[1], sample_ids=['s2'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        metrics = calculate_diversity_metrics(network, alignment)

        # π = 2 * p1 * p2 * d / L = 2 * 0.5 * 0.5 * 1 / 4 = 0.125
        assert metrics.nucleotide_diversity == pytest.approx(0.125)

    def test_all_unique_haplotypes(self):
        """Test maximum diversity (all unique)."""
        seqs = [
            Sequence(f'seq{i}', data)
            for i, data in enumerate(['ATCG', 'ATCC', 'GTCG', 'GTCC', 'TTAA'])
        ]
        haps = [Haplotype(seq, sample_ids=[f's{i}']) for i, seq in enumerate(seqs)]

        network = HaplotypeNetwork()
        for hap in haps:
            network.add_haplotype(hap)

        metrics = calculate_diversity_metrics(network)

        assert metrics.num_haplotypes == 5
        assert metrics.num_samples == 5
        # With all unique: H = (n/(n-1)) * (1 - 1/n) = (n/(n-1)) * ((n-1)/n) = 1
        assert metrics.haplotype_diversity == pytest.approx(1.0)


class TestCalculateNetworkMetrics:
    """Test network topology metrics."""

    def test_empty_network(self):
        """Test with empty network."""
        network = HaplotypeNetwork()
        metrics = calculate_network_metrics(network)

        assert metrics.diameter == 0
        assert metrics.num_components == 0
        assert metrics.central_haplotypes == []

    def test_single_node(self):
        """Test with single node."""
        seq = Sequence('seq1', 'ATCG')
        hap = Haplotype(seq, sample_ids=['s1'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap)

        metrics = calculate_network_metrics(network)

        assert metrics.diameter == 0
        assert metrics.num_components == 1
        assert metrics.central_haplotypes == [hap.id]

    def test_linear_network(self):
        """Test linear chain of nodes."""
        seqs = [
            Sequence(f'seq{i}', data)
            for i, data in enumerate(['ATCG', 'ATCC', 'GTCG', 'GTCC'])
        ]
        haps = [Haplotype(seq, sample_ids=[f's{i}']) for i, seq in enumerate(seqs)]

        network = HaplotypeNetwork()
        for hap in haps:
            network.add_haplotype(hap)

        # Create linear chain: hap0 - hap1 - hap2 - hap3
        for i in range(len(haps) - 1):
            network.add_edge(haps[i].id, haps[i + 1].id, distance=1)

        metrics = calculate_network_metrics(network)

        assert metrics.diameter == 3  # Longest path
        assert metrics.num_components == 1
        assert metrics.reticulation_index == 0.0  # It's a tree

    def test_star_network(self):
        """Test star topology."""
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

        metrics = calculate_network_metrics(network)

        assert metrics.diameter == 2  # Leaf to leaf through center
        assert metrics.num_components == 1
        assert (
            center.id in metrics.central_haplotypes[:1]
        )  # Center should be most central

    def test_network_with_reticulations(self):
        """Test network with alternative connections."""
        seqs = [
            Sequence(f'seq{i}', data) for i, data in enumerate(['ATCG', 'ATCC', 'GTCG'])
        ]
        haps = [Haplotype(seq, sample_ids=[f's{i}']) for i, seq in enumerate(seqs)]

        network = HaplotypeNetwork()
        for hap in haps:
            network.add_haplotype(hap)

        # Create triangle (tree would have 2 edges, this has 3)
        network.add_edge(haps[0].id, haps[1].id, distance=1)
        network.add_edge(haps[1].id, haps[2].id, distance=1)
        network.add_edge(haps[2].id, haps[0].id, distance=1)

        metrics = calculate_network_metrics(network)

        # Reticulation index = (actual_edges - tree_edges) / tree_edges
        # = (3 - 2) / 2 = 0.5
        assert metrics.reticulation_index == pytest.approx(0.5)

    def test_disconnected_network(self):
        """Test network with multiple components."""
        # Component 1
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATCC')
        hap1 = Haplotype(seq1, sample_ids=['s1'])
        hap2 = Haplotype(seq2, sample_ids=['s2'])

        # Component 2
        seq3 = Sequence('seq3', 'GTCG')
        hap3 = Haplotype(seq3, sample_ids=['s3'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_haplotype(hap3)
        network.add_edge(hap1.id, hap2.id, distance=1)

        metrics = calculate_network_metrics(network)

        assert metrics.num_components == 2


class TestIdentifyCentralHaplotypes:
    """Test central haplotype identification."""

    def test_degree_centrality(self):
        """Test identification by degree centrality."""
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

        central = identify_central_haplotypes(network, method='degree')

        assert central[0][0] == center.id  # Center should be most central
        assert central[0][1] > central[1][1]  # Center should have higher score

    def test_betweenness_centrality(self):
        """Test identification by betweenness centrality."""
        # Create linear chain where middle node has highest betweenness
        seqs = [
            Sequence(f'seq{i}', data)
            for i, data in enumerate(['ATCG', 'ATCC', 'GTCG', 'GTCC', 'TTAA'])
        ]
        haps = [Haplotype(seq, sample_ids=[f's{i}']) for i, seq in enumerate(seqs)]

        network = HaplotypeNetwork()
        for hap in haps:
            network.add_haplotype(hap)

        for i in range(len(haps) - 1):
            network.add_edge(haps[i].id, haps[i + 1].id, distance=1)

        central = identify_central_haplotypes(network, method='betweenness')

        # Middle node (index 2) should have highest betweenness
        assert central[0][0] == haps[2].id

    def test_empty_network(self):
        """Test with empty network."""
        network = HaplotypeNetwork()
        central = identify_central_haplotypes(network)

        assert central == []

    def test_invalid_method(self):
        """Test with invalid centrality method."""
        network = HaplotypeNetwork()
        seq = Sequence('seq1', 'ATCG')
        hap = Haplotype(seq, sample_ids=['s1'])
        network.add_haplotype(hap)

        with pytest.raises(ValueError, match='Unknown centrality method'):
            identify_central_haplotypes(network, method='invalid')


class TestCalculateReticulationIndex:
    """Test reticulation index calculation."""

    def test_empty_network(self):
        """Test with empty network."""
        network = HaplotypeNetwork()
        ri = calculate_reticulation_index(network)

        assert ri == 0.0

    def test_tree_network(self):
        """Test with tree (no reticulations)."""
        seqs = [
            Sequence(f'seq{i}', data) for i, data in enumerate(['ATCG', 'ATCC', 'GTCG'])
        ]
        haps = [Haplotype(seq, sample_ids=[f's{i}']) for i, seq in enumerate(seqs)]

        network = HaplotypeNetwork()
        for hap in haps:
            network.add_haplotype(hap)

        # Tree with 3 nodes should have 2 edges
        network.add_edge(haps[0].id, haps[1].id, distance=1)
        network.add_edge(haps[1].id, haps[2].id, distance=1)

        ri = calculate_reticulation_index(network)

        assert ri == 0.0

    def test_network_with_one_reticulation(self):
        """Test network with one extra edge."""
        seqs = [
            Sequence(f'seq{i}', data) for i, data in enumerate(['ATCG', 'ATCC', 'GTCG'])
        ]
        haps = [Haplotype(seq, sample_ids=[f's{i}']) for i, seq in enumerate(seqs)]

        network = HaplotypeNetwork()
        for hap in haps:
            network.add_haplotype(hap)

        # Create triangle
        network.add_edge(haps[0].id, haps[1].id, distance=1)
        network.add_edge(haps[1].id, haps[2].id, distance=1)
        network.add_edge(haps[2].id, haps[0].id, distance=1)

        ri = calculate_reticulation_index(network)

        # (3 edges - 2 tree edges) / 2 = 0.5
        assert ri == pytest.approx(0.5)


class TestGetFrequencyDistribution:
    """Test frequency distribution calculation."""

    def test_empty_network(self):
        """Test with empty network."""
        network = HaplotypeNetwork()
        dist = get_frequency_distribution(network)

        assert dist == {}

    def test_uniform_frequencies(self):
        """Test with all same frequency."""
        seqs = [
            Sequence(f'seq{i}', data) for i, data in enumerate(['ATCG', 'ATCC', 'GTCG'])
        ]
        haps = [
            Haplotype(seq, sample_ids=[f's{i}_{j}' for j in range(2)])
            for i, seq in enumerate(seqs)
        ]

        network = HaplotypeNetwork()
        for hap in haps:
            network.add_haplotype(hap)

        dist = get_frequency_distribution(network)

        assert dist[2] == 3  # Three haplotypes with frequency 2

    def test_varied_frequencies(self):
        """Test with varied frequencies."""
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATCC')
        seq3 = Sequence('seq3', 'GTCG')

        hap1 = Haplotype(seq1, sample_ids=['s1'])
        hap2 = Haplotype(seq2, sample_ids=['s2', 's3'])
        hap3 = Haplotype(seq3, sample_ids=['s4', 's5', 's6'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_haplotype(hap3)

        dist = get_frequency_distribution(network)

        assert dist[1] == 1  # One haplotype with frequency 1
        assert dist[2] == 1  # One haplotype with frequency 2
        assert dist[3] == 1  # One haplotype with frequency 3


class TestCalculateSummaryStatistics:
    """Test comprehensive summary statistics."""

    def test_empty_network(self):
        """Test with empty network."""
        network = HaplotypeNetwork()
        summary = calculate_summary_statistics(network)

        assert 'diversity' in summary
        assert 'network' in summary
        assert 'frequencies' in summary
        assert 'frequency_distribution' in summary

        assert summary['diversity']['num_haplotypes'] == 0
        assert summary['network']['num_components'] == 0

    def test_simple_network(self):
        """Test with simple network."""
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATCC')
        hap1 = Haplotype(seq1, sample_ids=['s1', 's2'])
        hap2 = Haplotype(seq2, sample_ids=['s3'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_edge(hap1.id, hap2.id, distance=1)

        summary = calculate_summary_statistics(network)

        assert summary['diversity']['num_haplotypes'] == 2
        assert summary['diversity']['num_samples'] == 3
        assert summary['network']['num_components'] == 1
        assert summary['frequencies']['total_samples'] == 3

    def test_with_alignment(self):
        """Test with alignment for additional metrics."""
        sequences = [Sequence('seq1', 'ATCG'), Sequence('seq2', 'ATCC')]
        alignment = Alignment(sequences)

        hap1 = Haplotype(sequences[0], sample_ids=['s1'])
        hap2 = Haplotype(sequences[1], sample_ids=['s2'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        summary = calculate_summary_statistics(network, alignment)

        assert summary['diversity']['nucleotide_diversity'] > 0
