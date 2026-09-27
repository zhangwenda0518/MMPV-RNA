"""End-to-end tests for the pypopart command-line interface."""

import json

from click.testing import CliRunner
import pytest

from pypopart.cli.main import main

FASTA = """>seq1
ATCGATCGAT
>seq2
ATCGATCGAT
>seq3
ATCGTTCGAT
>seq4
ATCGTTCGTT
>seq5
ATGGATCGAT
>seq6
ATGGATCGTT
"""

METADATA_CSV = """id,population
seq1,PopA
seq2,PopA
seq3,PopB
seq4,PopB
seq5,PopC
seq6,PopC
"""


@pytest.fixture
def runner():
    """Create a Click test runner."""
    return CliRunner()


@pytest.fixture
def fasta_file(tmp_path):
    """Write a small alignment to a temporary FASTA file."""
    path = tmp_path / 'seqs.fasta'
    path.write_text(FASTA)
    return str(path)


@pytest.fixture
def metadata_file(tmp_path):
    """Write matching metadata to a temporary CSV file."""
    path = tmp_path / 'meta.csv'
    path.write_text(METADATA_CSV)
    return str(path)


@pytest.fixture
def network_file(runner, fasta_file, tmp_path):
    """Build an MST network file to feed analyze/visualize tests."""
    out = tmp_path / 'net.graphml'
    result = runner.invoke(main, ['network', fasta_file, '-a', 'mst', '-o', str(out)])
    assert result.exit_code == 0, result.output
    return str(out)


class TestLoad:
    """Tests for `pypopart load`."""

    def test_load_basic(self, runner, fasta_file):
        """Load a FASTA file and report statistics."""
        result = runner.invoke(main, ['load', fasta_file])
        assert result.exit_code == 0, result.output
        assert 'Loaded 6 sequences' in result.output
        assert 'Alignment Statistics' in result.output

    def test_load_with_metadata(self, runner, fasta_file, metadata_file):
        """Metadata file is actually read and matched."""
        result = runner.invoke(main, ['load', fasta_file, '-m', metadata_file])
        assert result.exit_code == 0, result.output
        assert '6/6 sequences matched' in result.output

    def test_load_save_output(self, runner, fasta_file, tmp_path):
        """The -o option re-saves the alignment."""
        out = tmp_path / 'resaved.fasta'
        result = runner.invoke(main, ['load', fasta_file, '-o', str(out)])
        assert result.exit_code == 0, result.output
        assert out.exists()
        assert '>seq1' in out.read_text()

    def test_load_missing_file(self, runner):
        """A missing input file fails with a non-zero exit."""
        result = runner.invoke(main, ['load', 'no_such_file.fasta'])
        assert result.exit_code != 0

    def test_load_quiet(self, runner, fasta_file):
        """The -q flag suppresses informational output."""
        result = runner.invoke(main, ['-q', 'load', fasta_file])
        assert result.exit_code == 0, result.output
        assert 'Loaded' not in result.output


class TestNetwork:
    """Tests for `pypopart network`."""

    @pytest.mark.parametrize('algorithm', ['mst', 'msn', 'tcs', 'mjn', 'pn', 'tsw'])
    def test_network_each_algorithm(self, runner, fasta_file, tmp_path, algorithm):
        """Every advertised algorithm builds and saves a network."""
        out = tmp_path / f'{algorithm}.graphml'
        args = ['network', fasta_file, '-a', algorithm, '-o', str(out)]
        if algorithm == 'pn':
            args += ['--seed', '42']
        result = runner.invoke(main, args)
        assert result.exit_code == 0, result.output
        assert out.exists()
        assert 'Network constructed' in result.output

    @pytest.mark.parametrize('fmt', ['graphml', 'gml', 'json'])
    def test_network_output_formats(self, runner, fasta_file, tmp_path, fmt):
        """Each output format writes a file."""
        out = tmp_path / f'net.{fmt}'
        result = runner.invoke(
            main, ['network', fasta_file, '-a', 'mst', '-o', str(out), '--format', fmt]
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    @pytest.mark.parametrize('distance', ['hamming', 'jc', 'k2p', 'tn', 'tamura_nei'])
    def test_network_distance_metrics(self, runner, fasta_file, tmp_path, distance):
        """Every advertised distance metric is accepted."""
        out = tmp_path / 'net.graphml'
        result = runner.invoke(
            main, ['network', fasta_file, '-a', 'mst', '-d', distance, '-o', str(out)]
        )
        assert result.exit_code == 0, result.output

    def test_network_bad_algorithm(self, runner, fasta_file):
        """An unknown algorithm is rejected."""
        result = runner.invoke(main, ['network', fasta_file, '-a', 'bogus'])
        assert result.exit_code != 0

    def test_network_bad_distance(self, runner, fasta_file):
        """An unknown distance metric is rejected."""
        result = runner.invoke(main, ['network', fasta_file, '-d', 'bogus'])
        assert result.exit_code != 0


class TestAnalyze:
    """Tests for `pypopart analyze`."""

    def test_analyze_default_stats(self, runner, network_file):
        """Analyze with no flags prints statistics."""
        result = runner.invoke(main, ['analyze', network_file])
        assert result.exit_code == 0, result.output
        assert 'Network Statistics' in result.output

    def test_analyze_topology(self, runner, network_file):
        """The --topology flag prints topology analysis."""
        result = runner.invoke(main, ['analyze', network_file, '--topology'])
        assert result.exit_code == 0, result.output
        assert 'Topology Analysis' in result.output

    def test_analyze_popgen_with_alignment(self, runner, network_file, fasta_file):
        """The --popgen flag works with an alignment."""
        result = runner.invoke(
            main, ['analyze', network_file, '--popgen', '-a', fasta_file]
        )
        assert result.exit_code == 0, result.output
        assert "Tajima's D" in result.output

    def test_analyze_popgen_without_alignment_fails(self, runner, network_file):
        """The --popgen flag without alignment fails clearly."""
        result = runner.invoke(main, ['analyze', network_file, '--popgen'])
        assert result.exit_code != 0
        assert 'alignment' in result.output.lower()

    def test_analyze_json_output(self, runner, network_file, tmp_path):
        """Results are written as valid JSON."""
        out = tmp_path / 'results.json'
        result = runner.invoke(
            main, ['analyze', network_file, '--stats', '-o', str(out)]
        )
        assert result.exit_code == 0, result.output
        data = json.loads(out.read_text())
        assert 'statistics' in data


class TestVisualize:
    """Tests for `pypopart visualize`."""

    def test_visualize_static_png(self, runner, network_file, tmp_path):
        """Static PNG output is written."""
        pytest.importorskip('matplotlib')
        out = tmp_path / 'net.png'
        result = runner.invoke(main, ['visualize', network_file, '-o', str(out)])
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_visualize_interactive_html(self, runner, network_file, tmp_path):
        """Interactive HTML output is written."""
        pytest.importorskip('plotly')
        out = tmp_path / 'net.html'
        result = runner.invoke(main, ['visualize', network_file, '-o', str(out)])
        assert result.exit_code == 0, result.output
        assert out.exists()
        assert '<html' in out.read_text()[:200].lower()


class TestInfo:
    """Tests for `pypopart info`."""

    def test_list_algorithms_complete(self, runner):
        """All six algorithms are listed."""
        result = runner.invoke(main, ['info', '--list-algorithms'])
        assert result.exit_code == 0
        for name in ('mst', 'msn', 'tcs', 'mjn', 'pn', 'tsw'):
            assert name in result.output

    def test_list_distances(self, runner):
        """Distance metrics are listed."""
        result = runner.invoke(main, ['info', '--list-distances'])
        assert result.exit_code == 0
        for name in ('hamming', 'jc', 'k2p', 'tn'):
            assert name in result.output

    def test_list_formats(self, runner):
        """File formats are listed."""
        result = runner.invoke(main, ['info', '--list-formats'])
        assert result.exit_code == 0
        assert 'fasta' in result.output

    def test_no_geo_visualize_command(self, runner):
        """The removed geo-visualize command is gone."""
        result = runner.invoke(main, ['--help'])
        assert result.exit_code == 0
        assert 'geo-visualize' not in result.output


class TestRoundTrip:
    """Network files written by `network` load back identically."""

    @pytest.mark.parametrize('fmt', ['graphml', 'json'])
    def test_save_load_round_trip(self, runner, fasta_file, tmp_path, fmt):
        """Saved networks load back with usable distances."""
        from pypopart.io import load_network

        out = tmp_path / f'net.{fmt}'
        result = runner.invoke(
            main, ['network', fasta_file, '-a', 'mst', '-o', str(out), '--format', fmt]
        )
        assert result.exit_code == 0, result.output

        network = load_network(str(out))
        assert network.num_nodes > 0
        assert network.num_edges > 0
        # Every edge carries a usable distance
        for source, target in network.graph.edges():
            assert network.get_edge_distance(source, target) >= 0
