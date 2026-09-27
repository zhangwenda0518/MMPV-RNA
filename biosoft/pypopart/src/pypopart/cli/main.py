"""
Command-line interface for PyPopART.

This module provides the CLI for constructing and analyzing haplotype networks.
"""

from pathlib import Path
import traceback
from typing import Optional

import click

from pypopart import __version__

#: Exceptions the CLI reports as user-facing errors rather than tracebacks.
_USER_ERRORS = (OSError, ValueError, KeyError)


def _echo(ctx: click.Context, message: str = '') -> None:
    """
    Print an informational message unless --quiet was given.

    Parameters
    ----------
    ctx : click.Context
        Click context carrying the 'quiet' flag.
    message : str
        Message to print.
    """
    if not ctx.obj.get('quiet'):
        click.echo(message)


def _fail(ctx: click.Context, error: BaseException) -> None:
    """
    Report a fatal error and exit with status 1.

    Prints the full traceback first when -v/--verbose was given.

    Parameters
    ----------
    ctx : click.Context
        Click context carrying the 'verbose' count.
    error : BaseException
        The exception to report.

    Raises
    ------
    click.ClickException
        Always raised to terminate the command.
    """
    if ctx.obj.get('verbose', 0) > 0:
        traceback.print_exc()
    raise click.ClickException(str(error))


@click.group(
    help='PyPopART - Pure Python implementation of PopART haplotype network '
    'analysis.\n\nConstruct and visualize haplotype networks from DNA '
    'sequence alignments.'
)
@click.version_option(version=__version__)
@click.option(
    '-v',
    '--verbose',
    count=True,
    help='Increase verbosity (can be repeated: -v, -vv, -vvv)',
)
@click.option('-q', '--quiet', is_flag=True, help='Suppress all output except errors')
@click.pass_context
def main(ctx: click.Context, verbose: int, quiet: bool) -> None:
    """
    Cli: PyPopART - Pure Python implementation of PopART haplotype network analysis.

    Construct and visualize haplotype networks from DNA sequence alignments.

    Parameters
    ----------
    ctx : click.Context
        Click context object.
    verbose : int
        Verbosity level (repeatable flag count).
    quiet : bool
        Suppress all output except errors.
    """
    # Store verbosity in context for subcommands
    ctx.ensure_object(dict)
    ctx.obj['verbose'] = verbose
    ctx.obj['quiet'] = quiet


@main.command(
    help='Load and validate sequence alignment data.\n\nINPUT_FILE: Path to sequence alignment file'
)
@click.argument('input_file', type=click.Path(exists=True, dir_okay=False))
@click.option(
    '-f',
    '--format',
    type=click.Choice(['fasta', 'nexus', 'phylip', 'genbank'], case_sensitive=False),
    help='Input file format (auto-detected if not specified)',
)
@click.option(
    '-m',
    '--metadata',
    type=click.Path(exists=True, dir_okay=False),
    help='Metadata/traits file (CSV format)',
)
@click.option(
    '-o',
    '--output',
    type=click.Path(),
    help='Re-save the alignment (format from extension: .fasta/.nexus/.phy)',
)
@click.pass_context
def load(
    ctx: click.Context,
    input_file: str,
    format: Optional[str],
    metadata: Optional[str],
    output: Optional[str],
) -> None:
    """
    Load and validate sequence alignment data.

    INPUT_FILE: Path to sequence alignment file

    Parameters
    ----------
    ctx : click.Context
        Click context object.
    input_file : str
        Path to the sequence alignment file.
    format : str, optional
        Input format override (fasta, nexus, phylip, genbank).
    metadata : str, optional
        Path to a metadata CSV file.
    output : str, optional
        Path to re-save the alignment to.
    """
    from pypopart.io import load_alignment, save_alignment
    from pypopart.io.metadata import MetadataReader

    _echo(ctx, f'Loading sequences from {input_file}...')

    try:
        alignment = load_alignment(input_file, format=format)
        _echo(ctx, f'✓ Loaded {len(alignment)} sequences')
        _echo(ctx, f'  Alignment length: {alignment.length} bp')

        if metadata:
            _echo(ctx, f'Loading metadata from {metadata}...')
            metadata_dict = MetadataReader(metadata).read_metadata()
            matched = sum(
                1 for seq_id in alignment.sequence_ids if seq_id in metadata_dict
            )
            _echo(
                ctx,
                f'✓ Metadata loaded: {len(metadata_dict)} entries, '
                f'{matched}/{len(alignment)} sequences matched',
            )

        if output:
            suffix = Path(output).suffix.lower()
            out_format = {
                '.fasta': 'fasta',
                '.fa': 'fasta',
                '.fna': 'fasta',
                '.nexus': 'nexus',
                '.nex': 'nexus',
                '.phy': 'phylip',
                '.phylip': 'phylip',
            }.get(suffix, 'fasta')
            save_alignment(alignment, output, format=out_format)
            _echo(ctx, f'✓ Saved alignment to {output} ({out_format})')

        # Display summary statistics
        stats = alignment.calculate_stats()
        _echo(ctx, '\nAlignment Statistics:')
        _echo(ctx, f'  Sequences: {stats.num_sequences}')
        _echo(ctx, f'  Length: {stats.length} bp')
        _echo(ctx, f'  Variable sites: {stats.variable_sites}')
        _echo(ctx, f'  Parsimony informative: {stats.parsimony_informative_sites}')
        _echo(ctx, f'  GC content: {stats.gc_content:.1f}%')

    except _USER_ERRORS as e:
        _fail(ctx, e)


@main.command(
    help='Construct haplotype network from sequence alignment.\n\nINPUT_FILE: Path to sequence alignment file'
)
@click.argument('input_file', type=click.Path(exists=True, dir_okay=False))
@click.option(
    '-a',
    '--algorithm',
    type=click.Choice(['mst', 'msn', 'tcs', 'mjn', 'pn', 'tsw'], case_sensitive=False),
    default='mjn',
    show_default=True,
    help='Network construction algorithm',
)
@click.option(
    '-d',
    '--distance',
    type=click.Choice(
        ['hamming', 'jc', 'k2p', 'tn', 'tamura_nei'], case_sensitive=False
    ),
    default='hamming',
    show_default=True,
    help='Distance metric (tamura_nei is an alias for tn)',
)
@click.option(
    '-e',
    '--epsilon',
    type=float,
    default=0,
    show_default=True,
    help='Epsilon parameter for MSN/MJN',
)
@click.option(
    '-p',
    '--parsimony-limit',
    type=float,
    default=None,
    help='Enable a TCS connection limit at this parsimony confidence '
    '(0-1). PopART default is no limit (fully connected network).',
)
@click.option(
    '--seed',
    type=int,
    default=None,
    help='Random seed for stochastic algorithms (pn)',
)
@click.option('-o', '--output', type=click.Path(), help='Output network file')
@click.option(
    '--format',
    'output_format',
    type=click.Choice(['graphml', 'gml', 'json', 'nexus'], case_sensitive=False),
    default='graphml',
    show_default=True,
    help='Output format',
)
@click.pass_context
def network(
    ctx: click.Context,
    input_file: str,
    algorithm: str,
    distance: str,
    epsilon: float,
    parsimony_limit: float,
    seed: Optional[int],
    output: Optional[str],
    output_format: str,
) -> None:
    """
    Construct haplotype network from sequence alignment.

    INPUT_FILE: Path to sequence alignment file

    Parameters
    ----------
    ctx : click.Context
        Click context object.
    input_file : str
        Path to the sequence alignment file.
    algorithm : str
        Network construction algorithm name.
    distance : str
        Distance metric name or alias.
    epsilon : float
        Epsilon parameter for MSN/MJN.
    parsimony_limit : float
        Parsimony confidence limit for TCS (0-1).
    seed : int, optional
        Random seed for stochastic algorithms.
    output : str, optional
        Output network file path.
    output_format : str
        Output format (graphml, gml, json, nexus).
    """
    from pypopart.algorithms import build
    from pypopart.core.distance import normalize_distance_method
    from pypopart.core.haplotype import identify_haplotypes_from_alignment
    from pypopart.io import load_alignment, save_network

    _echo(ctx, f'Loading sequences from {input_file}...')

    try:
        distance = normalize_distance_method(distance)

        # Load alignment
        alignment = load_alignment(input_file)
        _echo(ctx, f'✓ Loaded {len(alignment)} sequences ({alignment.length} bp)')

        # Identify unique haplotypes for informational purposes
        _echo(ctx, 'Identifying unique haplotypes...')
        haplotypes = identify_haplotypes_from_alignment(alignment)
        _echo(ctx, f'✓ Found {len(haplotypes)} unique haplotypes')

        # Construct network
        _echo(ctx, f'Building {algorithm.upper()} network...')

        algorithm = algorithm.lower()
        algo_params = {}
        if algorithm in ('msn', 'mjn'):
            algo_params['epsilon'] = epsilon
        elif algorithm == 'tcs' and parsimony_limit is not None:
            algo_params['confidence'] = parsimony_limit
            algo_params['connection_limit'] = 'auto'
        elif algorithm == 'pn':
            algo_params['random_seed'] = seed
        algo = build(algorithm, distance_method=distance, **algo_params)

        network = algo.build_network(alignment)
        _echo(ctx, '✓ Network constructed')

        # Display network statistics
        _echo(ctx, '\nNetwork Statistics:')
        _echo(ctx, f'  Nodes: {len(network.graph.nodes)}')
        _echo(ctx, f'  Edges: {len(network.graph.edges)}')

        # Count inferred median/intermediate vectors
        n_medians = sum(
            1
            for _, attrs in network.graph.nodes(data=True)
            if attrs.get('median_vector')
        )
        if n_medians > 0:
            _echo(ctx, f'  Median vectors: {n_medians}')

        # Save network
        if output:
            _echo(ctx, f'\nSaving network to {output}...')
            save_network(network, output, format=output_format)
            _echo(ctx, f'✓ Network saved as {output_format.upper()}')
        else:
            click.echo('\nℹ Use -o/--output to save the network to a file', err=True)

    except _USER_ERRORS as e:
        _fail(ctx, e)


@main.command(
    help='Analyze haplotype network statistics.\n\nNETWORK_FILE: Path to network file (GraphML, GML, or JSON)'
)
@click.argument('network_file', type=click.Path(exists=True, dir_okay=False))
@click.option(
    '--stats',
    is_flag=True,
    help='Calculate and display network statistics',
)
@click.option(
    '--topology',
    is_flag=True,
    help='Analyze network topology',
)
@click.option(
    '--popgen',
    is_flag=True,
    help='Calculate population genetics measures (requires --alignment)',
)
@click.option(
    '-a',
    '--alignment',
    'alignment_file',
    type=click.Path(exists=True, dir_okay=False),
    help='Original sequence alignment (needed for --popgen and diversity)',
)
@click.option(
    '-o', '--output', type=click.Path(), help='Output file for analysis results (JSON)'
)
@click.pass_context
def analyze(
    ctx: click.Context,
    network_file: str,
    stats: bool,
    topology: bool,
    popgen: bool,
    alignment_file: Optional[str],
    output: Optional[str],
) -> None:
    """
    Analyze haplotype network statistics.

    NETWORK_FILE: Path to network file (GraphML, GML, or JSON)

    Parameters
    ----------
    ctx : click.Context
        Click context object.
    network_file : str
        Path to the network file.
    stats : bool
        Whether to print network statistics.
    topology : bool
        Whether to run topology analysis.
    popgen : bool
        Whether to run population genetics analysis.
    alignment_file : str, optional
        Path to the original alignment (required for --popgen).
    output : str, optional
        Path for a JSON results file.
    """
    from pypopart.io import load_alignment, load_network
    from pypopart.stats import (
        calculate_summary_statistics,
        calculate_tajimas_d,
        calculate_topology_summary,
    )

    _echo(ctx, f'Loading network from {network_file}...')

    try:
        network = load_network(network_file)
        _echo(ctx, f'✓ Loaded network with {network.num_nodes} nodes')

        alignment = None
        if alignment_file:
            alignment = load_alignment(alignment_file)
            _echo(ctx, f'✓ Loaded alignment with {len(alignment)} sequences')

        results = {}

        # Network statistics (default when no analysis flag is given)
        if stats or not (stats or topology or popgen):
            _echo(ctx, '\n=== Network Statistics ===')
            summary = calculate_summary_statistics(network, alignment)

            for section, values in summary.items():
                if isinstance(values, dict):
                    _echo(ctx, f'{section}:')
                    for key, value in values.items():
                        if isinstance(value, float):
                            _echo(ctx, f'  {key}: {value:.4f}')
                        else:
                            _echo(ctx, f'  {key}: {value}')
                else:
                    _echo(ctx, f'{section}: {values}')

            results['statistics'] = summary

        # Topology analysis
        if topology:
            _echo(ctx, '\n=== Topology Analysis ===')
            topo_summary = calculate_topology_summary(network)

            _echo(ctx, f'Connected: {topo_summary.get("is_connected")}')
            _echo(ctx, f'Components: {topo_summary.get("num_components")}')
            _echo(
                ctx, f'Star-like patterns: {len(topo_summary.get("star_patterns", []))}'
            )
            ancestral = topo_summary.get('ancestral_candidates', [])
            _echo(ctx, f'Ancestral candidates: {len(ancestral)}')

            results['topology'] = topo_summary

        # Population genetics
        if popgen:
            _echo(ctx, '\n=== Population Genetics ===')
            if alignment is None:
                raise ValueError(
                    'Population genetics analysis requires the original '
                    'alignment; pass it with -a/--alignment'
                )
            tajima = calculate_tajimas_d(alignment)
            _echo(ctx, f"Tajima's D: {tajima.D:.4f}")
            _echo(ctx, f'  Nucleotide diversity (pi): {tajima.pi:.4f}')
            _echo(ctx, f'  Watterson theta: {tajima.theta_w:.4f}')
            _echo(ctx, f'  Segregating sites: {tajima.n_segregating_sites}')

            results['popgen'] = {
                'tajimas_d': tajima.D,
                'pi': tajima.pi,
                'theta_w': tajima.theta_w,
                'segregating_sites': tajima.n_segregating_sites,
                'n_samples': tajima.n_samples,
            }

        # Save results
        if output:
            import json

            with open(output, 'w') as f:
                json.dump(results, f, indent=2, default=str)
            _echo(ctx, f'\n✓ Results saved to {output}')

    except _USER_ERRORS as e:
        _fail(ctx, e)


@main.command(
    help='Visualize haplotype network.\n\nNETWORK_FILE: Path to network file (GraphML, GML, or JSON)'
)
@click.argument('network_file', type=click.Path(exists=True, dir_okay=False))
@click.option(
    '-o',
    '--output',
    type=click.Path(),
    required=True,
    help='Output image file (PNG, PDF, SVG, or HTML)',
)
@click.option(
    '--layout',
    type=click.Choice(
        ['spring', 'circular', 'kamada_kawai'],
        case_sensitive=False,
    ),
    default='spring',
    show_default=True,
    help='Layout algorithm',
)
@click.option(
    '--width',
    type=int,
    default=800,
    show_default=True,
    help='Figure width in pixels',
)
@click.option(
    '--height',
    type=int,
    default=600,
    show_default=True,
    help='Figure height in pixels',
)
@click.option(
    '--interactive',
    is_flag=True,
    help='Create interactive HTML visualization (implied by .html output)',
)
@click.option(
    '--show-labels',
    is_flag=True,
    default=False,
    help='Show node labels',
)
@click.pass_context
def visualize(
    ctx: click.Context,
    network_file: str,
    output: str,
    layout: str,
    width: int,
    height: int,
    interactive: bool,
    show_labels: bool,
) -> None:
    """
    Visualize haplotype network.

    NETWORK_FILE: Path to network file (GraphML, GML, or JSON)

    Parameters
    ----------
    ctx : click.Context
        Click context object.
    network_file : str
        Path to the network file.
    output : str
        Output image file path.
    layout : str
        Layout algorithm name.
    width : int
        Figure width in pixels.
    height : int
        Figure height in pixels.
    interactive : bool
        Whether to force interactive HTML output.
    show_labels : bool
        Whether to show node labels.
    """
    from pypopart.io import load_network

    _echo(ctx, f'Loading network from {network_file}...')

    try:
        network = load_network(network_file)
        _echo(ctx, f'✓ Loaded network with {network.num_nodes} nodes')

        # Determine output format
        output_path = Path(output)
        is_html = output_path.suffix.lower() == '.html'

        if interactive or is_html:
            # Interactive visualization (plotly)
            from pypopart.visualization import InteractiveNetworkPlotter

            _echo(ctx, f'Creating interactive visualization with {layout} layout...')
            plotter = InteractiveNetworkPlotter(network)
            fig = plotter.plot(
                layout_algorithm=layout,
                width=width,
                height=height,
                show_labels=show_labels,
            )
            fig.write_html(str(output_path))
            _echo(ctx, f'✓ Interactive visualization saved to {output}')
            _echo(ctx, f'  Open in browser: file://{output_path.absolute()}')

        else:
            # Static visualization (matplotlib)
            from pypopart.visualization import StaticNetworkPlotter

            _echo(ctx, f'Creating static visualization with {layout} layout...')
            plotter = StaticNetworkPlotter(network)
            fig, _ax = plotter.plot(
                layout_algorithm=layout,
                figsize=(width / 100, height / 100),
                show_labels=show_labels,
            )
            fig.savefig(str(output_path), dpi=150, bbox_inches='tight')
            _echo(ctx, f'✓ Visualization saved to {output}')

    except _USER_ERRORS as e:
        _fail(ctx, e)


@main.command(help='Display information about PyPopART capabilities.')
@click.option(
    '--list-algorithms',
    is_flag=True,
    help='List available network construction algorithms',
)
@click.option(
    '--list-distances',
    is_flag=True,
    help='List available distance metrics',
)
@click.option(
    '--list-formats',
    is_flag=True,
    help='List supported file formats',
)
def info(
    list_algorithms: bool,
    list_distances: bool,
    list_formats: bool,
) -> None:
    """
    Display information about PyPopART capabilities.

    Parameters
    ----------
    list_algorithms : bool
        List available network construction algorithms.
    list_distances : bool
        List available distance metrics.
    list_formats : bool
        List supported file formats.
    """
    if list_algorithms:
        from pypopart.algorithms import list_algorithms as registry_list

        click.echo('Available Network Construction Algorithms:')
        for entry in registry_list():
            click.echo(f'  {entry["name"]:<4}- {entry["description"]}')
        click.echo()

    if list_distances:
        click.echo('Available Distance Metrics:')
        click.echo('  hamming    - Simple Hamming distance (count differences)')
        click.echo('  jc         - Jukes-Cantor correction')
        click.echo('  k2p        - Kimura 2-parameter model')
        click.echo('  tn         - Tamura-Nei model (alias: tamura_nei)')
        click.echo()

    if list_formats:
        click.echo('Supported Input Formats:')
        click.echo('  fasta    - FASTA sequence format')
        click.echo('  nexus    - NEXUS format')
        click.echo('  phylip   - PHYLIP format')
        click.echo('  genbank  - GenBank format')
        click.echo()
        click.echo('Supported Output Formats:')
        click.echo('  graphml  - GraphML (XML-based graph format)')
        click.echo('  gml      - Graph Modelling Language')
        click.echo('  json     - JSON format')
        click.echo('  nexus    - NEXUS format')
        click.echo()

    if not (list_algorithms or list_distances or list_formats):
        click.echo('Use --list-algorithms, --list-distances, or --list-formats')
        click.echo('or run "pypopart --help" for usage information')


if __name__ == '__main__':
    main()
