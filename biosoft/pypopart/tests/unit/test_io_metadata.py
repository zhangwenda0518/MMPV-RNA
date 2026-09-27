"""Unit tests for metadata file I/O."""

import pytest

from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence
from pypopart.io.metadata import MetadataReader, MetadataWriter


class TestMetadataReader:
    """Test cases for MetadataReader."""

    def test_read_simple_metadata(self, tmp_path):
        """Test reading simple metadata CSV."""
        csv_file = tmp_path / 'metadata.csv'
        csv_file.write_text('id,population,location\nseq1,A,NY\nseq2,B,CA\n')

        reader = MetadataReader(csv_file)
        metadata = reader.read_metadata()

        assert len(metadata) == 2
        assert metadata['seq1']['population'] == 'A'
        assert metadata['seq1']['location'] == 'NY'
        assert metadata['seq2']['population'] == 'B'

    def test_apply_to_alignment(self, tmp_path):
        """Test applying metadata to alignment."""
        csv_file = tmp_path / 'metadata.csv'
        csv_file.write_text('id,population\nseq1,A\nseq2,B\n')

        alignment = Alignment([Sequence('seq1', 'ATCG'), Sequence('seq2', 'GCTA')])

        reader = MetadataReader(csv_file)
        reader.apply_to_alignment(alignment)

        assert alignment[0].metadata['population'] == 'A'
        assert alignment[1].metadata['population'] == 'B'

    def test_custom_delimiter(self, tmp_path):
        """Test reading metadata with custom delimiter."""
        tsv_file = tmp_path / 'metadata.tsv'
        tsv_file.write_text('id\tpopulation\nseq1\tA\n')

        reader = MetadataReader(tsv_file, delimiter='\t')
        metadata = reader.read_metadata()

        assert metadata['seq1']['population'] == 'A'

    def test_missing_id_column(self, tmp_path):
        """Test error handling for missing ID column."""
        csv_file = tmp_path / 'metadata.csv'
        csv_file.write_text('name,population\nseq1,A\n')

        reader = MetadataReader(csv_file, id_column='id')
        with pytest.raises(ValueError, match='ID column'):
            reader.read_metadata()


class TestMetadataWriter:
    """Test cases for MetadataWriter."""

    def test_write_metadata(self, tmp_path):
        """Test writing metadata to CSV."""
        metadata = {
            'seq1': {'population': 'A', 'location': 'NY'},
            'seq2': {'population': 'B', 'location': 'CA'},
        }

        output_file = tmp_path / 'output.csv'
        writer = MetadataWriter(output_file)
        writer.write_metadata(metadata)

        assert output_file.exists()

        # Read back and verify
        reader = MetadataReader(output_file)
        read_metadata = reader.read_metadata()
        assert read_metadata['seq1']['population'] == 'A'

    def test_write_from_alignment(self, tmp_path):
        """Test extracting and writing metadata from alignment."""
        alignment = Alignment(
            [
                Sequence('seq1', 'ATCG', metadata={'pop': 'A'}),
                Sequence('seq2', 'GCTA', metadata={'pop': 'B'}),
            ]
        )

        output_file = tmp_path / 'output.csv'
        writer = MetadataWriter(output_file)
        writer.write_from_alignment(alignment)

        # Read back and verify
        reader = MetadataReader(output_file)
        metadata = reader.read_metadata()
        assert metadata['seq1']['pop'] == 'A'

    def test_trait_ordering(self, tmp_path):
        """Test writing metadata with specific trait order."""
        metadata = {'seq1': {'location': 'NY', 'population': 'A', 'year': '2020'}}

        output_file = tmp_path / 'output.csv'
        writer = MetadataWriter(output_file)
        writer.write_metadata(metadata, trait_order=['population', 'location'])

        # Verify column order
        content = output_file.read_text()
        header = content.split('\n')[0]
        assert header.index('population') < header.index('location')
