"""Unit tests for PHYLIP file I/O."""

import pytest

from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence
from pypopart.io.phylip import PhylipReader, PhylipWriter


class TestPhylipReader:
    """Test cases for PhylipReader."""

    def test_read_sequential_format(self, tmp_path):
        """Test reading sequential PHYLIP format."""
        phylip_file = tmp_path / 'test.phy'
        phylip_file.write_text(' 2 4\nseq1       ATCG\nseq2       ATCC\n')

        reader = PhylipReader(phylip_file)
        alignment = reader.read_alignment()

        assert len(alignment) == 2
        assert alignment.length == 4
        assert alignment[0].id == 'seq1'
        assert alignment[0].data == 'ATCG'

    def test_read_interleaved_format(self, tmp_path):
        """Test reading interleaved PHYLIP format."""
        phylip_file = tmp_path / 'test.phy'
        phylip_file.write_text(' 2 8\nseq1       ATCG\nseq2       GTCC\n\nATAT\nGCGC\n')

        reader = PhylipReader(phylip_file)
        alignment = reader.read_alignment()

        assert alignment.length == 8
        assert alignment[0].data == 'ATCGATAT'
        assert alignment[1].data == 'GTCCGCGC'

    def test_read_strict_format(self, tmp_path):
        """Test reading strict PHYLIP format (10-char IDs)."""
        phylip_file = tmp_path / 'test.phy'
        phylip_file.write_text(' 2 4\nseq1      ATCG\nseq2      ATCC\n')

        reader = PhylipReader(phylip_file, strict=True)
        alignment = reader.read_alignment()

        assert len(alignment) == 2
        assert alignment[0].id == 'seq1'

    def test_file_not_found(self):
        """Test error handling for missing file."""
        with pytest.raises(FileNotFoundError):
            PhylipReader('nonexistent.phy')


class TestPhylipWriter:
    """Test cases for PhylipWriter."""

    def test_write_sequential_format(self, tmp_path):
        """Test writing sequential PHYLIP format."""
        alignment = Alignment([Sequence('seq1', 'ATCG'), Sequence('seq2', 'ATCC')])

        output_file = tmp_path / 'output.phy'
        writer = PhylipWriter(output_file)
        writer.write_alignment(alignment)

        assert output_file.exists()

        # Read back and verify
        reader = PhylipReader(output_file)
        read_alignment = reader.read_alignment()
        assert len(read_alignment) == 2
        assert read_alignment.length == 4

    def test_write_interleaved_format(self, tmp_path):
        """Test writing interleaved PHYLIP format."""
        alignment = Alignment(
            [Sequence('seq1', 'A' * 100), Sequence('seq2', 'T' * 100)]
        )

        output_file = tmp_path / 'output.phy'
        writer = PhylipWriter(output_file, interleaved=True, line_length=50)
        writer.write_alignment(alignment)

        # Read back and verify
        reader = PhylipReader(output_file)
        read_alignment = reader.read_alignment()
        assert read_alignment.length == 100

    def test_write_strict_format(self, tmp_path):
        """Test writing strict PHYLIP format."""
        alignment = Alignment(
            [Sequence('verylongseqid', 'ATCG'), Sequence('seq2', 'ATCC')]
        )

        output_file = tmp_path / 'output.phy'
        writer = PhylipWriter(output_file, strict=True)
        writer.write_alignment(alignment)

        # Verify format
        content = output_file.read_text()
        lines = content.strip().split('\n')
        # In strict format, IDs should be truncated to 10 chars
        assert len(lines[1].split()[0]) <= 10


class TestPhylipFromString:
    """Reading PHYLIP content from in-memory text."""

    def test_from_string_roundtrip(self):
        """from_string parses the same content a file would."""
        from pypopart.io.phylip import PhylipReader

        alignment = PhylipReader.from_string(
            ' 2 4\ns1  ACGT\ns2  ACGA\n'
        ).read_alignment()
        assert len(alignment) == 2
        assert alignment[1].data == 'ACGA'
