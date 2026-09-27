"""Unit tests for population genetics module."""

from pypopart.core.alignment import Alignment
from pypopart.core.graph import HaplotypeNetwork
from pypopart.core.haplotype import Haplotype
from pypopart.core.sequence import Sequence
from pypopart.stats.popgen import (
    calculate_amova,
    calculate_fst_matrix,
    calculate_fu_fs,
    calculate_mismatch_distribution,
    calculate_pairwise_fst,
    calculate_tajimas_d,
)


class TestCalculateTajimasD:
    """Test Tajima's D calculation."""

    def test_empty_alignment(self):
        """Test with empty alignment."""
        alignment = Alignment()
        result = calculate_tajimas_d(alignment)

        assert result.D == 0.0
        assert result.n_samples == 0

    def test_single_sequence(self):
        """Test with single sequence."""
        seq = Sequence('seq1', 'ATCG')
        alignment = Alignment([seq])
        result = calculate_tajimas_d(alignment)

        assert result.D == 0.0
        assert result.n_samples == 1

    def test_no_variation(self):
        """Test with no segregating sites."""
        sequences = [
            Sequence('seq1', 'ATCG'),
            Sequence('seq2', 'ATCG'),
            Sequence('seq3', 'ATCG'),
        ]
        alignment = Alignment(sequences)
        result = calculate_tajimas_d(alignment)

        assert result.n_segregating_sites == 0
        assert result.D == 0.0

    def test_with_variation(self):
        """Test with some variation."""
        sequences = [
            Sequence('seq1', 'ATCG'),
            Sequence('seq2', 'ATCC'),
            Sequence('seq3', 'ATCG'),
        ]
        alignment = Alignment(sequences)
        result = calculate_tajimas_d(alignment)

        assert result.n_segregating_sites == 1
        assert result.n_samples == 3
        assert result.pi > 0
        assert result.theta_w > 0

    def test_all_different(self):
        """Test with all sequences different."""
        sequences = [
            Sequence('seq1', 'ATCG'),
            Sequence('seq2', 'GGGG'),
            Sequence('seq3', 'TTTT'),
        ]
        alignment = Alignment(sequences)
        result = calculate_tajimas_d(alignment)

        assert result.n_segregating_sites == 4
        assert result.pi > 0


class TestCalculateFuFs:
    """Test Fu's Fs calculation."""

    def test_empty_network(self):
        """Test with empty network."""
        network = HaplotypeNetwork()
        alignment = Alignment()
        result = calculate_fu_fs(network, alignment)

        assert result.Fs == 0.0
        assert result.n_haplotypes == 0
        assert result.n_samples == 0

    def test_single_haplotype(self):
        """Test with single haplotype."""
        seq = Sequence('seq1', 'ATCG')
        hap = Haplotype(seq, sample_ids=['s1', 's2'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap)

        sequences = [Sequence(f's{i}', 'ATCG') for i in range(1, 3)]
        alignment = Alignment(sequences)
        result = calculate_fu_fs(network, alignment)

        assert result.n_haplotypes == 1
        assert result.n_samples == 2

    def test_multiple_haplotypes(self):
        """Test with multiple haplotypes."""
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATCC')
        hap1 = Haplotype(seq1, sample_ids=['s1', 's2'])
        hap2 = Haplotype(seq2, sample_ids=['s3'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        sequences = [
            Sequence('s1', 'ATCG'),
            Sequence('s2', 'ATCG'),
            Sequence('s3', 'ATCC'),
        ]
        alignment = Alignment(sequences)
        result = calculate_fu_fs(network, alignment)

        assert result.n_haplotypes == 2
        assert result.n_samples == 3
        assert result.theta_pi >= 0


class TestCalculatePairwiseFst:
    """Test pairwise FST calculation."""

    def test_no_differentiation(self):
        """Test with no population differentiation."""
        seq1 = Sequence('seq1', 'ATCG')
        hap1 = Haplotype(
            seq1,
            sample_ids=['s1', 's2', 's3', 's4'],
            populations={'s1': 'PopA', 's2': 'PopA', 's3': 'PopB', 's4': 'PopB'},
        )

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)

        result = calculate_pairwise_fst(network, 'PopA', 'PopB')

        assert result.fst == 0.0
        assert result.pop1 == 'PopA'
        assert result.pop2 == 'PopB'

    def test_complete_differentiation(self):
        """Test with complete population differentiation."""
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'GGGG')

        hap1 = Haplotype(
            seq1, sample_ids=['s1', 's2'], populations={'s1': 'PopA', 's2': 'PopA'}
        )
        hap2 = Haplotype(
            seq2, sample_ids=['s3', 's4'], populations={'s3': 'PopB', 's4': 'PopB'}
        )

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        result = calculate_pairwise_fst(network, 'PopA', 'PopB')

        # Complete differentiation should give FST close to 1
        assert result.fst > 0.5

    def test_partial_differentiation(self):
        """Test with partial differentiation."""
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATCC')

        hap1 = Haplotype(
            seq1,
            sample_ids=['s1', 's2', 's3'],
            populations={'s1': 'PopA', 's2': 'PopA', 's3': 'PopB'},
        )
        hap2 = Haplotype(seq2, sample_ids=['s4'], populations={'s4': 'PopB'})

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        result = calculate_pairwise_fst(network, 'PopA', 'PopB')

        assert 0.0 < result.fst < 1.0

    def test_missing_population(self):
        """Test with missing population."""
        seq1 = Sequence('seq1', 'ATCG')
        hap1 = Haplotype(seq1, sample_ids=['s1'], populations={'s1': 'PopA'})

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)

        result = calculate_pairwise_fst(network, 'PopA', 'PopB')

        # Should return 0 when one population is missing
        assert result.fst == 0.0


class TestCalculateFstMatrix:
    """Test FST matrix calculation."""

    def test_single_population(self):
        """Test with single population."""
        seq1 = Sequence('seq1', 'ATCG')
        hap1 = Haplotype(
            seq1, sample_ids=['s1', 's2'], populations={'s1': 'PopA', 's2': 'PopA'}
        )

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)

        matrix = calculate_fst_matrix(network)

        assert matrix[('PopA', 'PopA')] == 0.0

    def test_multiple_populations(self):
        """Test with multiple populations."""
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'GGGG')
        seq3 = Sequence('seq3', 'TTTT')

        hap1 = Haplotype(seq1, sample_ids=['s1'], populations={'s1': 'PopA'})
        hap2 = Haplotype(seq2, sample_ids=['s2'], populations={'s2': 'PopB'})
        hap3 = Haplotype(seq3, sample_ids=['s3'], populations={'s3': 'PopC'})

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_haplotype(hap3)

        matrix = calculate_fst_matrix(network)

        # Check diagonal is 0
        assert matrix[('PopA', 'PopA')] == 0.0
        assert matrix[('PopB', 'PopB')] == 0.0
        assert matrix[('PopC', 'PopC')] == 0.0

        # Check symmetry
        assert matrix[('PopA', 'PopB')] == matrix[('PopB', 'PopA')]
        assert matrix[('PopA', 'PopC')] == matrix[('PopC', 'PopA')]
        assert matrix[('PopB', 'PopC')] == matrix[('PopC', 'PopB')]


class TestCalculateAMOVA:
    """Test AMOVA calculation."""

    def test_single_population(self):
        """Test with single population."""
        seq1 = Sequence('seq1', 'ATCG')
        hap1 = Haplotype(
            seq1, sample_ids=['s1', 's2'], populations={'s1': 'PopA', 's2': 'PopA'}
        )

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)

        sequences = [Sequence(f's{i}', 'ATCG') for i in range(1, 3)]
        alignment = Alignment(sequences)
        result = calculate_amova(network, alignment)

        # With single population, phi_st should be 0
        assert result.phi_st == 0.0

    def test_two_populations_no_variation(self):
        """Test with two populations but no variation."""
        seq1 = Sequence('seq1', 'ATCG')
        hap1 = Haplotype(
            seq1,
            sample_ids=['s1', 's2', 's3', 's4'],
            populations={'s1': 'PopA', 's2': 'PopA', 's3': 'PopB', 's4': 'PopB'},
        )

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)

        sequences = [Sequence(f's{i}', 'ATCG') for i in range(1, 5)]
        alignment = Alignment(sequences)
        result = calculate_amova(network, alignment)

        # No variation means phi_st should be 0
        assert result.phi_st == 0.0

    def test_two_populations_with_variation(self):
        """Test with two populations and variation."""
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

        sequences = [
            Sequence('s1', 'ATCG'),
            Sequence('s2', 'ATCG'),
            Sequence('s3', 'ATCC'),
            Sequence('s4', 'ATCC'),
        ]
        alignment = Alignment(sequences)
        result = calculate_amova(network, alignment)

        # With complete differentiation, phi_st should be high
        assert result.phi_st > 0.0
        assert result.variance_among_pops >= 0
        assert result.variance_within_pops >= 0


class TestCalculateMismatchDistribution:
    """Test mismatch distribution calculation."""

    def test_single_haplotype(self):
        """Test with single haplotype."""
        seq1 = Sequence('seq1', 'ATCG')
        hap1 = Haplotype(seq1, sample_ids=['s1', 's2', 's3'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)

        dist = calculate_mismatch_distribution(network)

        # All comparisons have 0 differences
        assert dist[0] == 3  # 3 samples means 3 comparisons with self (0 diff)

    def test_two_haplotypes(self):
        """Test with two haplotypes."""
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATCC')
        hap1 = Haplotype(seq1, sample_ids=['s1'])
        hap2 = Haplotype(seq2, sample_ids=['s2'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        dist = calculate_mismatch_distribution(network)

        # Two haplotypes with 1 difference
        assert 1 in dist
        assert dist[1] == 1  # One comparison with 1 difference

    def test_multiple_haplotypes(self):
        """Test with multiple haplotypes."""
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATCC')
        seq3 = Sequence('seq3', 'GGGG')

        hap1 = Haplotype(seq1, sample_ids=['s1'])
        hap2 = Haplotype(seq2, sample_ids=['s2'])
        hap3 = Haplotype(seq3, sample_ids=['s3'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)
        network.add_haplotype(hap3)

        dist = calculate_mismatch_distribution(network)

        # Should have various numbers of differences
        assert len(dist) > 0
        assert all(k >= 0 for k in dist.keys())
        assert all(v >= 0 for v in dist.values())

    def test_with_frequency_weighting(self):
        """Test that frequencies are properly weighted."""
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATCC')

        hap1 = Haplotype(seq1, sample_ids=['s1', 's2'])
        hap2 = Haplotype(seq2, sample_ids=['s3'])

        network = HaplotypeNetwork()
        network.add_haplotype(hap1)
        network.add_haplotype(hap2)

        dist = calculate_mismatch_distribution(network)

        # hap1 vs hap2: freq 2 * freq 1 = 2 comparisons with 1 difference
        assert dist[1] == 2
