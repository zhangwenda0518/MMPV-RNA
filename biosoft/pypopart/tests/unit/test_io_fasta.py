"""Unit tests for FASTA file I/O."""

import gzip

import pytest

from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence
from pypopart.io.fasta import FastaReader, FastaWriter


class TestFastaReader:
    """Test cases for FastaReader."""

    def test_read_simple_fasta(self, tmp_path):
        """Test reading a simple FASTA file."""
        fasta_file = tmp_path / 'test.fasta'
        fasta_file.write_text('>seq1\nATCGATCG\n>seq2\nGCTAGCTA\n')

        reader = FastaReader(fasta_file)
        sequences = list(reader.read_sequences())

        assert len(sequences) == 2
        assert sequences[0].id == 'seq1'
        assert sequences[0].data == 'ATCGATCG'
        assert sequences[1].id == 'seq2'
        assert sequences[1].data == 'GCTAGCTA'

    def test_read_fasta_with_metadata(self, tmp_path):
        """Test reading FASTA with metadata in header."""
        fasta_file = tmp_path / 'test.fasta'
        fasta_file.write_text('>seq1|pop=A|location=NY\nATCG\n')

        reader = FastaReader(fasta_file)
        sequences = list(reader.read_sequences())

        assert len(sequences) == 1
        assert sequences[0].metadata['pop'] == 'A'
        assert sequences[0].metadata['location'] == 'NY'

    def test_read_alignment(self, tmp_path):
        """Test reading alignment from FASTA."""
        fasta_file = tmp_path / 'test.fasta'
        fasta_file.write_text('>seq1\nATCG\n>seq2\nATCC\n')

        reader = FastaReader(fasta_file)
        alignment = reader.read_alignment()

        assert len(alignment) == 2
        assert alignment.length == 4

    def test_read_gzip_compressed(self, tmp_path):
        """Test reading gzip compressed FASTA."""
        fasta_file = tmp_path / 'test.fasta.gz'

        with gzip.open(fasta_file, 'wt') as f:
            f.write('>seq1\nATCG\n')

        reader = FastaReader(fasta_file)
        sequences = list(reader.read_sequences())

        assert len(sequences) == 1
        assert sequences[0].data == 'ATCG'

    def test_file_not_found(self):
        """Test error handling for missing file."""
        with pytest.raises(FileNotFoundError):
            FastaReader('nonexistent.fasta')

    def test_validation_error(self, tmp_path):
        """Test validation of invalid sequences."""
        fasta_file = tmp_path / 'test.fasta'
        fasta_file.write_text('>seq1\nXXXX\n')  # Invalid characters

        reader = FastaReader(fasta_file, validate=True)
        with pytest.raises(ValueError):
            list(reader.read_sequences())


class TestFastaWriter:
    """Test cases for FastaWriter."""

    def test_write_sequences(self, tmp_path):
        """Test writing sequences to FASTA."""
        sequences = [Sequence('seq1', 'ATCG'), Sequence('seq2', 'GCTA')]

        output_file = tmp_path / 'output.fasta'
        writer = FastaWriter(output_file)
        count = writer.write_sequences(iter(sequences))

        assert count == 2
        assert output_file.exists()

        # Read back and verify
        reader = FastaReader(output_file)
        read_seqs = list(reader.read_sequences())
        assert len(read_seqs) == 2
        assert read_seqs[0].data == 'ATCG'

    def test_write_with_metadata(self, tmp_path):
        """Test writing sequences with metadata."""
        seq = Sequence('seq1', 'ATCG', metadata={'pop': 'A'})

        output_file = tmp_path / 'output.fasta'
        writer = FastaWriter(output_file)
        writer.write_sequences([seq])

        # Read back and verify metadata
        reader = FastaReader(output_file)
        read_seqs = list(reader.read_sequences())
        assert read_seqs[0].metadata['pop'] == 'A'

    def test_write_alignment(self, tmp_path):
        """Test writing alignment to FASTA."""
        alignment = Alignment([Sequence('seq1', 'ATCG'), Sequence('seq2', 'ATCC')])

        output_file = tmp_path / 'output.fasta'
        writer = FastaWriter(output_file)
        count = writer.write_alignment(alignment)

        assert count == 2
        assert output_file.exists()

    def test_write_with_line_wrapping(self, tmp_path):
        """Test writing with line length limit."""
        seq = Sequence('seq1', 'A' * 100)

        output_file = tmp_path / 'output.fasta'
        writer = FastaWriter(output_file, line_length=50)
        writer.write_sequences([seq])

        # Verify wrapping
        lines = output_file.read_text().strip().split('\n')
        assert len(lines) == 3  # Header + 2 wrapped lines
        assert len(lines[1]) == 50

    def test_write_gzip_compressed(self, tmp_path):
        """Test writing gzip compressed FASTA."""
        seq = Sequence('seq1', 'ATCG')

        output_file = tmp_path / 'output.fasta'
        writer = FastaWriter(output_file, compress='gzip')
        writer.write_sequences([seq])

        # Verify compression
        assert output_file.with_suffix('.fasta.gz').exists()

        # Read back
        reader = FastaReader(output_file.with_suffix('.fasta.gz'))
        read_seqs = list(reader.read_sequences())
        assert len(read_seqs) == 1


class TestFastaFromString:
    """Reading FASTA content from in-memory text."""

    def test_from_string_roundtrip(self):
        """from_string parses the same content a file would."""
        from pypopart.io.fasta import FastaReader

        alignment = FastaReader.from_string('>s1\nACGT\n>s2\nACGA\n').read_alignment()
        assert len(alignment) == 2
        assert alignment[0].id == 's1'
        assert alignment[0].data == 'ACGT'
