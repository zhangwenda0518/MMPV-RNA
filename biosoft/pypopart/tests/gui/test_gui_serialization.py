"""Tests for the GUI's dcc.Store serialization and upload format sniffing."""

from pypopart.core.graph import HaplotypeNetwork
from pypopart.core.haplotype import Haplotype
from pypopart.core.sequence import Sequence
from pypopart.gui.serialization import network_to_store, store_to_networkx


def build_network() -> HaplotypeNetwork:
    """Build a small network with one median vertex."""
    network = HaplotypeNetwork()
    network.add_haplotype(Haplotype(Sequence('H1', 'AAT'), sample_ids=['s1', 's2']))
    network.add_haplotype(Haplotype(Sequence('H2', 'ATA'), sample_ids=['s3']))
    network.add_haplotype(
        Haplotype(Sequence('Median_0', 'AAA'), sample_ids=[]), median_vector=True
    )
    network.add_edge('H1', 'Median_0', distance=1)
    network.add_edge('H2', 'Median_0', distance=1)
    return network


class TestStoreRoundTrip:
    """network_to_store payloads survive the store round trip."""

    def test_round_trip_preserves_medians(self):
        """The is_median flag survives serialize -> from_serialized."""
        store = network_to_store(build_network())
        rebuilt = HaplotypeNetwork.from_serialized(store)

        assert rebuilt.num_nodes == 3
        assert rebuilt.num_edges == 2
        assert rebuilt.is_median_vector('Median_0')
        assert not rebuilt.is_median_vector('H1')

    def test_round_trip_preserves_distances(self):
        """Edge distances survive the round trip."""
        store = network_to_store(build_network())
        rebuilt = HaplotypeNetwork.from_serialized(store)
        assert rebuilt.get_edge_distance('H1', 'Median_0') == 1

    def test_store_to_networkx(self):
        """The raw-graph helper carries node and edge attributes."""
        graph = store_to_networkx(network_to_store(build_network()))
        assert graph.number_of_nodes() == 3
        assert graph.nodes['Median_0']['is_median'] is True
        assert graph['H1']['Median_0']['distance'] == 1


class TestUploadFormatSniffing:
    """The upload handler detects formats by content, not extension."""

    def test_detects_fasta_by_content(self):
        """FASTA content is recognised regardless of file name."""
        from pypopart.gui.callbacks.upload import _detect_reader
        from pypopart.io.fasta import FastaReader

        assert _detect_reader('>s1\nACGT\n', 'data.txt') is FastaReader

    def test_detects_nexus_by_content(self):
        """NEXUS content is recognised by its #NEXUS header."""
        from pypopart.gui.callbacks.upload import _detect_reader
        from pypopart.io.nexus import NexusReader

        assert _detect_reader('#NEXUS\nBEGIN DATA;\n', 'foo.dat') is NexusReader

    def test_detects_phylip_by_content(self):
        """PHYLIP content is recognised by its dimensions header."""
        from pypopart.gui.callbacks.upload import _detect_reader
        from pypopart.io.phylip import PhylipReader

        assert _detect_reader(' 2 4\ns1 ACGT\ns2 ACGA\n', 'x.txt') is PhylipReader

    def test_falls_back_to_extension(self):
        """Ambiguous content falls back to the file extension."""
        from pypopart.gui.callbacks.upload import _detect_reader
        from pypopart.io.fasta import FastaReader

        assert _detect_reader('', 'seqs.fasta') is FastaReader

    def test_unknown_returns_none(self):
        """Undetectable input returns None (handler shows an error)."""
        from pypopart.gui.callbacks.upload import _detect_reader

        assert _detect_reader('garbage', 'file.xyz') is None
