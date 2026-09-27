"""Shared fixtures for the PyPopART test suite."""

import pytest

from pypopart.core.alignment import Alignment
from pypopart.core.graph import HaplotypeNetwork
from pypopart.core.haplotype import Haplotype
from pypopart.core.sequence import Sequence

SMALL_FASTA = """>seq1
ATCGATCGAT
>seq2
ATCGATCGAT
>seq3
ATCGTTCGAT
>seq4
ATCGTTCGTT
"""


@pytest.fixture
def small_alignment() -> Alignment:
    """Four-sequence alignment with three unique haplotypes."""
    return Alignment(
        [
            Sequence('seq1', 'ATCGATCGAT'),
            Sequence('seq2', 'ATCGATCGAT'),
            Sequence('seq3', 'ATCGTTCGAT'),
            Sequence('seq4', 'ATCGTTCGTT'),
        ]
    )


@pytest.fixture
def alignment_with_gaps() -> Alignment:
    """Alignment containing gap and ambiguity characters."""
    return Alignment(
        [
            Sequence('seq1', 'AT-GATNGAT'),
            Sequence('seq2', 'ATCGATCGAT'),
            Sequence('seq3', 'ATCGTTCG-T'),
        ]
    )


@pytest.fixture
def tiny_network() -> HaplotypeNetwork:
    """Three-node network with one median vertex."""
    network = HaplotypeNetwork()
    network.add_haplotype(Haplotype(Sequence('H1', 'AAT'), sample_ids=['s1', 's2']))
    network.add_haplotype(Haplotype(Sequence('H2', 'ATA'), sample_ids=['s3']))
    network.add_haplotype(
        Haplotype(Sequence('Median_0', 'AAA'), sample_ids=[]), median_vector=True
    )
    network.add_edge('H1', 'Median_0', distance=1)
    network.add_edge('H2', 'Median_0', distance=1)
    return network


@pytest.fixture
def fasta_file(tmp_path):
    """Write the small alignment to a temporary FASTA file."""
    path = tmp_path / 'seqs.fasta'
    path.write_text(SMALL_FASTA)
    return str(path)
