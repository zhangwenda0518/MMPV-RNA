"""Unit tests for Haplotype class and related functions."""

import pytest

from pypopart.core.alignment import Alignment
from pypopart.core.haplotype import (
    Haplotype,
    HaplotypeFrequency,
    calculate_haplotype_diversity,
    identify_haplotypes_from_alignment,
)
from pypopart.core.sequence import Sequence


class TestHaplotype:
    """Test cases for Haplotype class."""

    def test_basic_haplotype_creation(self):
        """Test basic haplotype creation."""
        seq = Sequence('H1', 'ATCG')
        haplotype = Haplotype(seq)

        assert haplotype.id == 'H1'
        assert haplotype.data == 'ATCG'
        assert haplotype.frequency == 0
        assert haplotype.sample_ids == []

    def test_haplotype_with_samples(self):
        """Test haplotype with sample IDs."""
        seq = Sequence('H1', 'ATCG')
        samples = ['sample1', 'sample2', 'sample3']
        haplotype = Haplotype(seq, sample_ids=samples)

        assert haplotype.frequency == 3
        assert set(haplotype.sample_ids) == set(samples)

    def test_haplotype_with_populations(self):
        """Test haplotype with population assignments."""
        seq = Sequence('H1', 'ATCG')
        samples = ['s1', 's2', 's3']
        populations = {'s1': 'PopA', 's2': 'PopA', 's3': 'PopB'}

        haplotype = Haplotype(seq, sample_ids=samples, populations=populations)

        assert haplotype.get_population('s1') == 'PopA'
        assert haplotype.get_population('s3') == 'PopB'
        assert haplotype.get_populations() == {'PopA', 'PopB'}

    def test_add_sample(self):
        """Test adding samples to haplotype."""
        seq = Sequence('H1', 'ATCG')
        haplotype = Haplotype(seq)

        haplotype.add_sample('s1', 'PopA')
        assert haplotype.frequency == 1
        assert haplotype.get_population('s1') == 'PopA'

        haplotype.add_sample('s2', 'PopB')
        assert haplotype.frequency == 2
        assert haplotype.get_populations() == {'PopA', 'PopB'}

    def test_add_duplicate_sample(self):
        """Test adding duplicate sample (should not increase frequency)."""
        seq = Sequence('H1', 'ATCG')
        haplotype = Haplotype(seq)

        haplotype.add_sample('s1')
        haplotype.add_sample('s1')  # Duplicate

        assert haplotype.frequency == 1  # Set removes duplicates

    def test_remove_sample(self):
        """Test removing samples from haplotype."""
        seq = Sequence('H1', 'ATCG')
        haplotype = Haplotype(seq, sample_ids=['s1', 's2'], populations={'s1': 'PopA'})

        haplotype.remove_sample('s1')
        assert haplotype.frequency == 1
        assert 's1' not in haplotype.sample_ids
        assert haplotype.get_population('s1') is None

    def test_remove_nonexistent_sample(self):
        """Test removing non-existent sample fails."""
        seq = Sequence('H1', 'ATCG')
        haplotype = Haplotype(seq)

        with pytest.raises(KeyError, match='not in haplotype'):
            haplotype.remove_sample('nonexistent')

    def test_get_frequency_by_population(self):
        """Test per-population frequency calculation."""
        seq = Sequence('H1', 'ATCG')
        samples = ['s1', 's2', 's3', 's4']
        populations = {'s1': 'PopA', 's2': 'PopA', 's3': 'PopB', 's4': 'PopB'}

        haplotype = Haplotype(seq, sample_ids=samples, populations=populations)
        freq_by_pop = haplotype.get_frequency_by_population()

        assert freq_by_pop['PopA'] == 2
        assert freq_by_pop['PopB'] == 2

    def test_get_frequency_by_population_unassigned(self):
        """Test frequency with unassigned samples."""
        seq = Sequence('H1', 'ATCG')
        samples = ['s1', 's2', 's3']
        populations = {'s1': 'PopA'}  # Only s1 assigned

        haplotype = Haplotype(seq, sample_ids=samples, populations=populations)
        freq_by_pop = haplotype.get_frequency_by_population()

        assert freq_by_pop['PopA'] == 1
        assert freq_by_pop['Unassigned'] == 2

    def test_get_frequency_info(self):
        """Test comprehensive frequency information."""
        seq = Sequence('H1', 'ATCG')
        samples = ['s1', 's2', 's3']
        populations = {'s1': 'PopA', 's2': 'PopA', 's3': 'PopB'}

        haplotype = Haplotype(seq, sample_ids=samples, populations=populations)
        freq_info = haplotype.get_frequency_info()

        assert isinstance(freq_info, HaplotypeFrequency)
        assert freq_info.total == 3
        assert freq_info.by_population['PopA'] == 2
        assert freq_info.by_population['PopB'] == 1

    def test_haplotype_equality(self):
        """Test haplotype equality based on sequence."""
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCG')  # Same data, different ID
        seq3 = Sequence('H3', 'GCTA')  # Different data

        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)
        hap3 = Haplotype(seq3)

        assert hap1 == hap2  # Same sequence data
        assert hap1 != hap3

    def test_haplotype_hash(self):
        """Test haplotype can be used in sets/dicts."""
        seq1 = Sequence('H1', 'ATCG')
        seq2 = Sequence('H2', 'ATCG')  # Same data
        seq3 = Sequence('H3', 'GCTA')  # Different data

        hap1 = Haplotype(seq1)
        hap2 = Haplotype(seq2)
        hap3 = Haplotype(seq3)

        hap_set = {hap1, hap2, hap3}
        assert len(hap_set) == 2  # hap1 and hap2 are equal

    def test_to_dict(self):
        """Test conversion to dictionary."""
        seq = Sequence('H1', 'ATCG')
        samples = ['s1', 's2']
        populations = {'s1': 'PopA', 's2': 'PopB'}

        haplotype = Haplotype(seq, sample_ids=samples, populations=populations)
        hap_dict = haplotype.to_dict()

        assert hap_dict['id'] == 'H1'
        assert hap_dict['frequency'] == 2
        assert set(hap_dict['sample_ids']) == set(samples)
        assert hap_dict['populations'] == populations
        assert 'frequency_by_population' in hap_dict

    def test_from_dict(self):
        """Test creation from dictionary."""
        hap_dict = {
            'sequence': {
                'id': 'H1',
                'data': 'ATCG',
                'metadata': {},
                'description': None,
            },
            'sample_ids': ['s1', 's2'],
            'populations': {'s1': 'PopA', 's2': 'PopB'},
        }

        haplotype = Haplotype.from_dict(hap_dict)

        assert haplotype.id == 'H1'
        assert haplotype.data == 'ATCG'
        assert haplotype.frequency == 2

    def test_string_representation(self):
        """Test string representation."""
        seq = Sequence('H1', 'ATCG')
        samples = ['s1', 's2', 's3']
        populations = {'s1': 'PopA', 's2': 'PopA', 's3': 'PopB'}

        haplotype = Haplotype(seq, sample_ids=samples, populations=populations)
        str_repr = str(haplotype)

        assert 'H1' in str_repr
        assert '3 samples' in str_repr
        assert '2 populations' in str_repr


class TestIdentifyHaplotypes:
    """Test haplotype identification from alignment."""

    def test_identify_haplotypes_basic(self):
        """Test basic haplotype identification."""
        sequences = [
            Sequence('seq1', 'ATCG'),
            Sequence('seq2', 'ATCG'),  # Same as seq1
            Sequence('seq3', 'GCTA'),
            Sequence('seq4', 'ATCG'),  # Same as seq1
        ]
        alignment = Alignment(sequences)

        haplotypes = identify_haplotypes_from_alignment(alignment)

        assert len(haplotypes) == 2
        # Find haplotype with 3 samples
        hap_atcg = next(h for h in haplotypes if h.frequency == 3)
        assert hap_atcg.data == 'ATCG'
        assert set(hap_atcg.sample_ids) == {'seq1', 'seq2', 'seq4'}

    def test_identify_haplotypes_with_gaps(self):
        """Gapped sequences dedup by exact aligned string (PopART semantics)."""
        sequences = [
            Sequence('seq1', 'AT-CG'),
            Sequence('seq2', 'ATC-G'),  # Same ungapped, different aligned string
            Sequence('seq3', 'AT-CG'),  # Identical aligned string to seq1
        ]
        alignment = Alignment(sequences)

        haplotypes = identify_haplotypes_from_alignment(alignment)

        # seq1/seq3 group; seq2 stays distinct despite equal ungapped data
        assert len(haplotypes) == 2
        hap = next(h for h in haplotypes if 'seq1' in h.sample_ids)
        assert hap.frequency == 2
        assert set(hap.sample_ids) == {'seq1', 'seq3'}
        # Haplotype sequences keep the alignment length (gaps included)
        assert all(len(h.data) == 5 for h in haplotypes)

    def test_identify_haplotypes_with_populations(self):
        """Test haplotype identification with population map."""
        sequences = [
            Sequence('seq1', 'ATCG'),
            Sequence('seq2', 'ATCG'),
            Sequence('seq3', 'GCTA'),
        ]
        alignment = Alignment(sequences)

        pop_map = {'seq1': 'PopA', 'seq2': 'PopA', 'seq3': 'PopB'}

        haplotypes = identify_haplotypes_from_alignment(alignment, pop_map)

        hap_atcg = next(h for h in haplotypes if h.frequency == 2)
        assert hap_atcg.get_populations() == {'PopA'}

        hap_gcta = next(h for h in haplotypes if h.frequency == 1)
        assert hap_gcta.get_populations() == {'PopB'}

    def test_identify_haplotypes_all_unique(self):
        """Test with all unique sequences."""
        sequences = [
            Sequence('seq1', 'ATCG'),
            Sequence('seq2', 'GCTA'),
            Sequence('seq3', 'CGAT'),
        ]
        alignment = Alignment(sequences)

        haplotypes = identify_haplotypes_from_alignment(alignment)

        assert len(haplotypes) == 3
        assert all(h.frequency == 1 for h in haplotypes)

    def test_identify_haplotypes_all_same(self):
        """Test with all identical sequences."""
        sequences = [
            Sequence('seq1', 'ATCG'),
            Sequence('seq2', 'ATCG'),
            Sequence('seq3', 'ATCG'),
        ]
        alignment = Alignment(sequences)

        haplotypes = identify_haplotypes_from_alignment(alignment)

        assert len(haplotypes) == 1
        assert haplotypes[0].frequency == 3

    def test_identify_haplotypes_numbering(self):
        """Test haplotypes are numbered correctly."""
        sequences = [
            Sequence('seq1', 'GCTA'),  # Will be H1 (alphabetically first)
            Sequence('seq2', 'ATCG'),  # Will be H2
        ]
        alignment = Alignment(sequences)

        haplotypes = identify_haplotypes_from_alignment(alignment)

        # Check IDs are H1, H2
        hap_ids = [h.id for h in haplotypes]
        assert 'H1' in hap_ids
        assert 'H2' in hap_ids


class TestHaplotypeDiversity:
    """Test haplotype diversity calculations."""

    def test_diversity_empty(self):
        """Test diversity with no haplotypes."""
        diversity = calculate_haplotype_diversity([])

        assert diversity['num_haplotypes'] == 0
        assert diversity['total_samples'] == 0
        assert diversity['haplotype_diversity'] == 0.0
        assert diversity['singleton_count'] == 0

    def test_diversity_single_haplotype(self):
        """Test diversity with single haplotype."""
        seq = Sequence('H1', 'ATCG')
        haplotype = Haplotype(seq, sample_ids=['s1', 's2', 's3'])

        diversity = calculate_haplotype_diversity([haplotype])

        assert diversity['num_haplotypes'] == 1
        assert diversity['total_samples'] == 3
        assert diversity['haplotype_diversity'] == 0.0  # No diversity
        assert diversity['singleton_count'] == 0

    def test_diversity_all_unique(self):
        """Test diversity with all unique haplotypes (maximum diversity)."""
        haplotypes = [
            Haplotype(Sequence('H1', 'ATCG'), sample_ids=['s1']),
            Haplotype(Sequence('H2', 'GCTA'), sample_ids=['s2']),
            Haplotype(Sequence('H3', 'CGAT'), sample_ids=['s3']),
        ]

        diversity = calculate_haplotype_diversity(haplotypes)

        assert diversity['num_haplotypes'] == 3
        assert diversity['total_samples'] == 3
        assert diversity['haplotype_diversity'] == 1.0  # Maximum diversity
        assert diversity['singleton_count'] == 3  # All singletons

    def test_diversity_mixed_frequencies(self):
        """Test diversity with mixed frequencies."""
        haplotypes = [
            Haplotype(Sequence('H1', 'ATCG'), sample_ids=['s1', 's2', 's3']),
            Haplotype(Sequence('H2', 'GCTA'), sample_ids=['s4']),
            Haplotype(Sequence('H3', 'CGAT'), sample_ids=['s5']),
        ]

        diversity = calculate_haplotype_diversity(haplotypes)

        assert diversity['num_haplotypes'] == 3
        assert diversity['total_samples'] == 5
        assert 0 < diversity['haplotype_diversity'] < 1.0
        assert diversity['singleton_count'] == 2  # H2 and H3

    def test_diversity_calculation_formula(self):
        """Test diversity formula with known values."""
        # Two haplotypes: 3 and 2 samples (total 5)
        # H = (5/4) * (1 - ((3/5)^2 + (2/5)^2))
        # H = 1.25 * (1 - (0.36 + 0.16))
        # H = 1.25 * 0.48 = 0.6
        haplotypes = [
            Haplotype(Sequence('H1', 'ATCG'), sample_ids=['s1', 's2', 's3']),
            Haplotype(Sequence('H2', 'GCTA'), sample_ids=['s4', 's5']),
        ]

        diversity = calculate_haplotype_diversity(haplotypes)

        assert diversity['num_haplotypes'] == 2
        assert diversity['total_samples'] == 5
        assert abs(diversity['haplotype_diversity'] - 0.6) < 0.001
        assert diversity['singleton_count'] == 0
