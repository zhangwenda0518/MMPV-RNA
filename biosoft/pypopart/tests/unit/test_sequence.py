"""Unit tests for Sequence class."""

import pytest

from pypopart.core.sequence import Sequence


class TestSequence:
    """Test cases for Sequence class."""

    def test_basic_sequence_creation(self):
        """Test basic sequence creation."""
        seq = Sequence(id='test1', data='ATCG')
        assert seq.id == 'test1'
        assert seq.data == 'ATCG'
        assert len(seq) == 4
        assert seq.metadata == {}
        assert seq.description is None

    def test_sequence_with_metadata(self):
        """Test sequence creation with metadata."""
        metadata = {'population': 'A', 'location': 'site1'}
        seq = Sequence(
            id='test2', data='ATCG', metadata=metadata, description='Test sequence'
        )
        assert seq.metadata == metadata
        assert seq.description == 'Test sequence'

    def test_sequence_normalization(self):
        """Test sequence data is normalized (uppercase, stripped)."""
        seq = Sequence(id='test3', data='  atcg  ')
        assert seq.data == 'ATCG'

    def test_sequence_validation_valid(self):
        """Test validation passes for valid sequences."""
        # All IUPAC nucleotide codes should be valid
        valid_seq = 'ACGTRYSWKMBDHVN-?'
        seq = Sequence(id='valid', data=valid_seq)
        assert seq.data == valid_seq

    def test_sequence_validation_invalid(self):
        """Test validation fails for invalid characters."""
        with pytest.raises(ValueError, match='Invalid characters'):
            Sequence(id='invalid', data='ATCGX')

    def test_sequence_equality(self):
        """Test sequence equality comparison."""
        seq1 = Sequence(id='seq1', data='ATCG')
        seq2 = Sequence(id='seq2', data='ATCG')  # different ID, same data
        seq3 = Sequence(id='seq1', data='GCTA')  # same ID, different data

        assert seq1 == seq2  # equality based on data
        assert seq1 != seq3
        assert seq2 != seq3

    def test_sequence_hash(self):
        """Test sequence can be used in sets/dicts."""
        seq1 = Sequence(id='seq1', data='ATCG')
        seq2 = Sequence(id='seq2', data='ATCG')  # same data
        seq3 = Sequence(id='seq3', data='GCTA')  # different data

        seq_set = {seq1, seq2, seq3}
        assert len(seq_set) == 2  # seq1 and seq2 should be considered equal

    def test_reverse_complement(self):
        """Test reverse complement calculation."""
        seq = Sequence(id='test', data='ATCG')
        rev_comp = seq.reverse_complement()

        assert rev_comp.data == 'CGAT'
        assert rev_comp.id == 'test_rev_comp'
        assert 'Reverse complement' in rev_comp.description

    def test_reverse_complement_ambiguous(self):
        """Test reverse complement with ambiguous characters."""
        seq = Sequence(id='test', data='ATCGRYSW')
        rev_comp = seq.reverse_complement()

        assert rev_comp.data == 'WSRYCGAT'

    def test_gc_content(self):
        """Test GC content calculation."""
        # 50% GC content
        seq1 = Sequence(id='test1', data='ATCG')
        assert seq1.gc_content() == 0.5

        # 100% GC content
        seq2 = Sequence(id='test2', data='GCGC')
        assert seq2.gc_content() == 1.0

        # 0% GC content
        seq3 = Sequence(id='test3', data='ATAT')
        assert seq3.gc_content() == 0.0

        # With gaps (should ignore gaps)
        seq4 = Sequence(id='test4', data='AT-CG')
        assert seq4.gc_content() == 0.5

    def test_gap_counting(self):
        """Test gap counting."""
        seq = Sequence(id='test', data='AT-C-G')
        assert seq.count_gaps() == 2

    def test_ambiguous_counting(self):
        """Test ambiguous character counting."""
        seq = Sequence(id='test', data='ATCGN?')
        assert seq.count_ambiguous() == 2

    def test_remove_gaps(self):
        """Test gap removal."""
        seq = Sequence(id='test', data='AT-C-G')
        ungapped = seq.remove_gaps()

        assert ungapped.data == 'ATCG'
        assert ungapped.id == 'test'  # ID should be preserved
        assert ungapped.metadata == seq.metadata

    def test_sequence_slice(self):
        """Test sequence slicing."""
        seq = Sequence(id='test', data='ATCGTA')

        # Slice with start and end
        slice1 = seq.slice(1, 4)
        assert slice1.data == 'TCG'
        assert 'slice' in slice1.id.lower()

        # Slice with start only
        slice2 = seq.slice(2)
        assert slice2.data == 'CGTA'

    def test_to_dict(self):
        """Test conversion to dictionary."""
        seq = Sequence(
            id='test', data='ATCG', metadata={'pop': 'A'}, description='Test sequence'
        )

        seq_dict = seq.to_dict()

        assert seq_dict['id'] == 'test'
        assert seq_dict['data'] == 'ATCG'
        assert seq_dict['metadata'] == {'pop': 'A'}
        assert seq_dict['description'] == 'Test sequence'
        assert seq_dict['length'] == 4
        assert seq_dict['gc_content'] == 0.5
        assert seq_dict['gaps'] == 0
        assert seq_dict['ambiguous'] == 0

    def test_from_dict(self):
        """Test creation from dictionary."""
        seq_dict = {
            'id': 'test',
            'data': 'ATCG',
            'metadata': {'pop': 'A'},
            'description': 'Test sequence',
        }

        seq = Sequence.from_dict(seq_dict)

        assert seq.id == 'test'
        assert seq.data == 'ATCG'
        assert seq.metadata == {'pop': 'A'}
        assert seq.description == 'Test sequence'

    def test_string_representation(self):
        """Test string representation (FASTA format)."""
        seq = Sequence(id='test', data='ATCG')
        str_repr = str(seq)

        assert str_repr.startswith('>test')
        assert 'ATCG' in str_repr
