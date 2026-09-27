"""Unit tests for Alignment class."""

import numpy as np
import pytest

from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence


class TestAlignment:
    """Test cases for Alignment class."""

    def test_empty_alignment(self):
        """Test empty alignment creation."""
        alignment = Alignment()
        assert len(alignment) == 0
        assert alignment.length == 0
        assert alignment.sequence_ids == []
        assert alignment.is_valid()

    def test_alignment_with_sequences(self):
        """Test alignment creation with sequences."""
        sequences = [
            Sequence('seq1', 'ATCG'),
            Sequence('seq2', 'ATCC'),
            Sequence('seq3', 'GTCG'),
        ]
        alignment = Alignment(sequences)

        assert len(alignment) == 3
        assert alignment.length == 4
        assert alignment.sequence_ids == ['seq1', 'seq2', 'seq3']

    def test_add_sequence(self):
        """Test adding sequences to alignment."""
        alignment = Alignment()
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'GTCC')

        alignment.add_sequence(seq1)
        assert len(alignment) == 1

        alignment.add_sequence(seq2)
        assert len(alignment) == 2

    def test_add_sequence_length_mismatch(self):
        """Test adding sequence with wrong length fails."""
        alignment = Alignment()
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATCGAA')  # Different length

        alignment.add_sequence(seq1)

        with pytest.raises(ValueError, match="doesn't match alignment length"):
            alignment.add_sequence(seq2)

    def test_add_duplicate_id(self):
        """Test adding sequence with duplicate ID fails."""
        alignment = Alignment()
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq1', 'GTCC')  # Same ID

        alignment.add_sequence(seq1)

        with pytest.raises(ValueError, match='already exists'):
            alignment.add_sequence(seq2)

    def test_get_sequence(self):
        """Test retrieving sequences by ID."""
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'GTCC')
        alignment = Alignment([seq1, seq2])

        retrieved = alignment.get_sequence('seq1')
        assert retrieved == seq1
        assert retrieved.data == 'ATCG'

    def test_get_sequence_not_found(self):
        """Test retrieving non-existent sequence fails."""
        alignment = Alignment([Sequence('seq1', 'ATCG')])

        with pytest.raises(KeyError, match='not found'):
            alignment.get_sequence('nonexistent')

    def test_remove_sequence(self):
        """Test removing sequences."""
        sequences = [
            Sequence('seq1', 'ATCG'),
            Sequence('seq2', 'GTCC'),
            Sequence('seq3', 'TTCG'),
        ]
        alignment = Alignment(sequences)

        alignment.remove_sequence('seq2')
        assert len(alignment) == 2
        assert 'seq2' not in alignment.sequence_ids
        assert alignment.sequence_ids == ['seq1', 'seq3']

    def test_remove_sequence_not_found(self):
        """Test removing non-existent sequence fails."""
        alignment = Alignment([Sequence('seq1', 'ATCG')])

        with pytest.raises(KeyError, match='not found'):
            alignment.remove_sequence('nonexistent')

    def test_indexing_by_int(self):
        """Test indexing alignment by integer."""
        sequences = [Sequence('seq1', 'ATCG'), Sequence('seq2', 'GTCC')]
        alignment = Alignment(sequences)

        assert alignment[0] == sequences[0]
        assert alignment[1] == sequences[1]
        assert alignment[-1] == sequences[1]

    def test_indexing_by_id(self):
        """Test indexing alignment by sequence ID."""
        seq1 = Sequence('seq1', 'ATCG')
        alignment = Alignment([seq1])

        assert alignment['seq1'] == seq1

    def test_slicing(self):
        """Test slicing alignment."""
        sequences = [
            Sequence('seq1', 'ATCG'),
            Sequence('seq2', 'GTCC'),
            Sequence('seq3', 'TTCG'),
        ]
        alignment = Alignment(sequences)

        subset = alignment[1:3]
        assert len(subset) == 2
        assert subset.sequence_ids == ['seq2', 'seq3']

    def test_iteration(self):
        """Test iterating over alignment."""
        sequences = [Sequence('seq1', 'ATCG'), Sequence('seq2', 'GTCC')]
        alignment = Alignment(sequences)

        iterated_sequences = list(alignment)
        assert iterated_sequences == sequences

    def test_get_column(self):
        """Test getting alignment columns."""
        sequences = [
            Sequence('seq1', 'ATCG'),
            Sequence('seq2', 'GTCC'),
            Sequence('seq3', 'ATCG'),
        ]
        alignment = Alignment(sequences)

        col0 = alignment.get_column(0)
        assert col0 == ['A', 'G', 'A']

        col3 = alignment.get_column(3)
        assert col3 == ['G', 'C', 'G']

    def test_get_column_out_of_range(self):
        """Test getting column with invalid index."""
        alignment = Alignment([Sequence('seq1', 'ATCG')])

        with pytest.raises(IndexError, match='out of range'):
            alignment.get_column(10)

    def test_slice_alignment(self):
        """Test slicing alignment by positions."""
        sequences = [Sequence('seq1', 'ATCGAA'), Sequence('seq2', 'GTCCGG')]
        alignment = Alignment(sequences)

        # Slice positions 1-4
        sliced = alignment.slice_alignment(1, 4)
        assert sliced.length == 3
        assert sliced[0].data == 'TCG'
        assert sliced[1].data == 'TCC'

    def test_validation(self):
        """Test alignment validation."""
        # Valid alignment
        sequences = [Sequence('seq1', 'ATCG'), Sequence('seq2', 'GTCC')]
        alignment = Alignment(sequences)
        assert alignment.is_valid()
        alignment.validate()  # Should not raise

        # Invalid alignment (manually create)
        alignment._sequences = [
            Sequence('seq1', 'ATCG'),
            Sequence('seq2', 'GTCCAA'),  # Different length
        ]
        assert not alignment.is_valid()

        with pytest.raises(ValueError, match='different lengths'):
            alignment.validate()

    def test_remove_gap_columns(self):
        """Test removing columns with gaps."""
        sequences = [
            Sequence('seq1', 'AT-G'),
            Sequence('seq2', 'GT-C'),
            Sequence('seq3', 'AT-G'),
        ]
        alignment = Alignment(sequences)

        # Remove columns with >=50% gaps (100% gaps in this case for column 2)
        no_gaps = alignment.remove_gaps_columns(0.5)
        assert no_gaps.length == 3  # Column 2 (100% gaps) should be removed
        assert no_gaps[0].data == 'ATG'
        assert no_gaps[1].data == 'GTC'

    def test_calculate_stats(self):
        """Test alignment statistics calculation."""
        sequences = [
            Sequence('seq1', 'ATCG'),  # All different at position 0
            Sequence('seq2', 'GTCG'),  # All same at position 3
            Sequence('seq3', 'CTCG'),
            Sequence('seq4', 'TTCG'),
        ]
        alignment = Alignment(sequences)

        stats = alignment.calculate_stats()

        assert stats.length == 4
        assert stats.num_sequences == 4
        assert stats.conserved_sites == 3  # Positions 1,2,3
        assert stats.variable_sites == 1  # Position 0
        assert (
            stats.parsimony_informative_sites == 0
        )  # No position has ≥2 chars appearing ≥2 times

    def test_calculate_stats_empty(self):
        """Test statistics for empty alignment."""
        alignment = Alignment()
        stats = alignment.calculate_stats()

        assert stats.length == 0
        assert stats.num_sequences == 0
        assert stats.total_gaps == 0

    def test_hamming_distance(self):
        """Test Hamming distance calculation."""
        alignment = Alignment()
        seq1 = Sequence('seq1', 'ATCG')
        seq2 = Sequence('seq2', 'ATCC')  # 1 difference

        distance = alignment._hamming_distance(seq1, seq2)
        assert distance == 1

    def test_distance_matrix(self):
        """Test distance matrix calculation."""
        sequences = [
            Sequence('seq1', 'ATCG'),
            Sequence('seq2', 'ATCC'),  # 1 diff from seq1
            Sequence('seq3', 'GTCG'),  # 1 diff from seq1, 2 diff from seq2
        ]
        alignment = Alignment(sequences)

        matrix = alignment.get_distance_matrix()

        assert matrix.shape == (3, 3)
        assert matrix[0, 1] == matrix[1, 0] == 1  # seq1 vs seq2
        assert matrix[0, 2] == matrix[2, 0] == 1  # seq1 vs seq3
        assert matrix[1, 2] == matrix[2, 1] == 2  # seq2 vs seq3
        assert np.all(np.diag(matrix) == 0)

    def test_identify_haplotypes(self):
        """Test haplotype identification."""
        sequences = [
            Sequence('seq1', 'ATCG'),
            Sequence('seq2', 'ATCG'),  # Same as seq1
            Sequence('seq3', 'GTCC'),
            Sequence('seq4', 'ATCG'),  # Same as seq1 and seq2
        ]
        alignment = Alignment(sequences)

        haplotypes = alignment.identify_haplotypes()

        assert len(haplotypes) == 2
        assert set(haplotypes['ATCG']) == {'seq1', 'seq2', 'seq4'}
        assert haplotypes['GTCC'] == ['seq3']

    def test_to_fasta(self):
        """Test FASTA format conversion."""
        sequences = [Sequence('seq1', 'ATCG'), Sequence('seq2', 'GTCC')]
        alignment = Alignment(sequences)

        fasta = alignment.to_fasta()

        assert '>seq1' in fasta
        assert '>seq2' in fasta
        assert 'ATCG' in fasta
        assert 'GTCC' in fasta

    def test_string_representation(self):
        """Test string representation."""
        sequences = [Sequence('seq1', 'ATCG'), Sequence('seq2', 'GTCC')]
        alignment = Alignment(sequences)

        str_repr = str(alignment)
        assert '2 sequences' in str_repr
        assert '4 positions' in str_repr
