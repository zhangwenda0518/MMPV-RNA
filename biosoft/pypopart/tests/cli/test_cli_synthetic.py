"""End-to-end CLI runs on the structured synthetic dataset."""

import csv
import json
import re

from click.testing import CliRunner
import pytest

from pypopart.cli.main import main
from pypopart.io.fasta import FastaWriter
from pypopart.io.nexus import NexusWriter
from pypopart.io.phylip import PhylipWriter
from pypopart.io.synthetic import POPULATIONS, simulate_structured_alignment


@pytest.fixture(scope='module')
def dataset(tmp_path_factory):
    """Write the synthetic alignment in every format, once per module."""
    out = tmp_path_factory.mktemp('synthetic')
    alignment, populations = simulate_structured_alignment()

    FastaWriter(out / 'synthetic.fasta').write_alignment(alignment)
    NexusWriter(out / 'synthetic.nex').write_alignment(alignment)
    PhylipWriter(out / 'synthetic.phy').write_alignment(alignment)
    with open(out / 'meta.csv', 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['id', 'population'])
        writer.writerows(populations.items())

    return {
        'dir': out,
        'fasta': str(out / 'synthetic.fasta'),
        'nexus': str(out / 'synthetic.nex'),
        'phylip': str(out / 'synthetic.phy'),
        'metadata': str(out / 'meta.csv'),
        'n_seqs': len(alignment),
        'length': alignment.length,
        'n_haplotypes': len({s.data for s in alignment}),
    }


@pytest.fixture
def runner():
    """Create a Click test runner."""
    return CliRunner()


def _count(pattern: str, text: str) -> int:
    """Pull the first integer captured by a regex out of CLI output."""
    match = re.search(pattern, text)
    assert match, f'{pattern!r} not found in:\n{text}'
    return int(match.group(1))


class TestLoad:
    """Every format reads back the same alignment."""

    @pytest.mark.parametrize('fmt', ['fasta', 'nexus', 'phylip'])
    def test_load_each_format(self, runner, dataset, fmt):
        """Sequence count and length match the generator."""
        result = runner.invoke(main, ['load', dataset[fmt]])

        assert result.exit_code == 0, result.output
        assert _count(r'Loaded (\d+) sequences', result.output) == dataset['n_seqs']
        assert _count(r'Alignment length: (\d+) bp', result.output) == dataset['length']

    def test_load_with_metadata(self, runner, dataset):
        """Every sample has a population."""
        result = runner.invoke(
            main, ['load', dataset['fasta'], '-m', dataset['metadata']]
        )

        assert result.exit_code == 0, result.output
        assert (
            _count(r'Metadata loaded: (\d+) entries', result.output)
            == (dataset['n_seqs'])
        )

    def test_convert_between_formats(self, runner, dataset, tmp_path):
        """FASTA in, NEXUS out, and the NEXUS loads with the same counts."""
        out = tmp_path / 'converted.nex'
        result = runner.invoke(main, ['load', dataset['fasta'], '-o', str(out)])
        assert result.exit_code == 0, result.output

        result = runner.invoke(main, ['load', str(out)])
        assert result.exit_code == 0, result.output
        assert _count(r'Loaded (\d+) sequences', result.output) == dataset['n_seqs']


class TestNetwork:
    """Every algorithm builds a network over the full haplotype set."""

    @pytest.mark.parametrize('algorithm', ['mst', 'msn', 'tcs', 'mjn', 'tsw', 'pn'])
    def test_each_algorithm(self, runner, dataset, tmp_path, algorithm):
        """Node count is the haplotype count plus any median vectors."""
        out = tmp_path / f'{algorithm}.graphml'
        args = ['network', dataset['nexus'], '-a', algorithm, '-o', str(out)]
        if algorithm == 'pn':
            args += ['--seed', '1']
        result = runner.invoke(main, args)

        assert result.exit_code == 0, result.output
        assert (
            _count(r'Found (\d+) unique haplotypes', result.output)
            == (dataset['n_haplotypes'])
        )
        nodes = _count(r'Nodes: (\d+)', result.output)
        medians = (
            _count(r'Median vectors: (\d+)', result.output)
            if 'Median vectors' in result.output
            else 0
        )
        assert nodes == dataset['n_haplotypes'] + medians
        assert out.exists() and out.stat().st_size > 0

    def test_mjn_reconstructs_median_vectors(self, runner, dataset, tmp_path):
        """The Desert reticulation forces at least one median vector."""
        out = tmp_path / 'mjn.graphml'
        result = runner.invoke(
            main, ['network', dataset['fasta'], '-a', 'mjn', '-o', str(out)]
        )

        assert result.exit_code == 0, result.output
        assert _count(r'Median vectors: (\d+)', result.output) >= 1

    @pytest.mark.parametrize('fmt', ['graphml', 'gml', 'json', 'nexus'])
    def test_output_formats(self, runner, dataset, tmp_path, fmt):
        """Each export format writes a non-empty file."""
        out = tmp_path / f'net.{fmt}'
        result = runner.invoke(
            main,
            ['network', dataset['fasta'], '-a', 'mst', '-o', str(out), '--format', fmt],
        )

        assert result.exit_code == 0, result.output
        assert out.stat().st_size > 0


class TestAnalyze:
    """Analysis on a real-sized network."""

    @pytest.fixture(scope='class')
    def network_file(self, dataset, tmp_path_factory):
        """Build one MJN network for the analysis tests."""
        out = tmp_path_factory.mktemp('net') / 'mjn.graphml'
        result = CliRunner().invoke(
            main, ['network', dataset['fasta'], '-a', 'mjn', '-o', str(out)]
        )
        assert result.exit_code == 0, result.output
        return str(out)

    def test_stats_and_topology(self, runner, network_file):
        """The network is connected and reports its statistics."""
        result = runner.invoke(main, ['analyze', network_file, '--stats', '--topology'])

        assert result.exit_code == 0, result.output
        assert 'Network Statistics' in result.output
        assert 'Connected: True' in result.output

    def test_popgen_writes_json(self, runner, dataset, network_file, tmp_path):
        """Population genetics needs the alignment and lands in JSON."""
        out = tmp_path / 'stats.json'
        result = runner.invoke(
            main,
            [
                'analyze',
                network_file,
                '--popgen',
                '-a',
                dataset['fasta'],
                '-o',
                str(out),
            ],
        )

        assert result.exit_code == 0, result.output
        report = json.loads(out.read_text())
        assert report


class TestVisualize:
    """Figures render from a real-sized network."""

    @pytest.fixture(scope='class')
    def network_file(self, dataset, tmp_path_factory):
        """Build one MJN network for the figure tests."""
        out = tmp_path_factory.mktemp('viz') / 'mjn.graphml'
        result = CliRunner().invoke(
            main, ['network', dataset['fasta'], '-a', 'mjn', '-o', str(out)]
        )
        assert result.exit_code == 0, result.output
        return str(out)

    def test_static_png(self, runner, network_file, tmp_path):
        """A PNG with labels renders and is not empty."""
        pytest.importorskip('matplotlib')
        out = tmp_path / 'net.png'
        result = runner.invoke(
            main, ['visualize', network_file, '-o', str(out), '--show-labels']
        )

        assert result.exit_code == 0, result.output
        assert out.stat().st_size > 1000

    def test_interactive_html(self, runner, network_file, tmp_path):
        """An HTML figure renders."""
        pytest.importorskip('plotly')
        out = tmp_path / 'net.html'
        result = runner.invoke(main, ['visualize', network_file, '-o', str(out)])

        assert result.exit_code == 0, result.output
        assert out.stat().st_size > 1000


def test_populations_constant_matches_metadata(dataset):
    """The metadata file names exactly the documented populations."""
    with open(dataset['metadata']) as handle:
        rows = list(csv.DictReader(handle))

    assert {row['population'] for row in rows} == set(POPULATIONS)
