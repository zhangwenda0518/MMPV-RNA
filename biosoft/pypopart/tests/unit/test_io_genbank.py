"""Tests for the GenBank reader."""

import pytest

from pypopart.io.genbank import GenBankReader

GENBANK = """LOCUS       TEST1                     8 bp    DNA     linear   UNA 01-JAN-2020
DEFINITION  Test sequence one.
ACCESSION   TEST1
ORIGIN
        1 atcgatcg
//
LOCUS       TEST2                     8 bp    DNA     linear   UNA 01-JAN-2020
DEFINITION  Test sequence two.
ACCESSION   TEST2
ORIGIN
        1 atcgatca
//
"""


@pytest.fixture
def genbank_file(tmp_path):
    """Write a two-record GenBank file."""
    path = tmp_path / 'seqs.gb'
    path.write_text(GENBANK)
    return str(path)


class TestGenBankReader:
    """GenBank parsing basics."""

    def test_read_sequences(self, genbank_file):
        """Both records parse with upper-case sequence data."""
        sequences = list(GenBankReader(genbank_file).read_sequences())
        assert len(sequences) == 2
        assert sequences[0].id == 'TEST1'
        assert sequences[0].data == 'ATCGATCG'
        assert sequences[1].data == 'ATCGATCA'

    def test_missing_file_raises(self):
        """A missing path fails at construction time."""
        with pytest.raises(FileNotFoundError):
            GenBankReader('no_such_file.gb')

    def test_load_alignment_via_dispatcher(self, genbank_file):
        """load_alignment auto-detects .gb files."""
        from pypopart.io import load_alignment

        alignment = load_alignment(genbank_file)
        assert len(alignment) == 2
