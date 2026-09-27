"""
Dash-based GUI for PyPopART haplotype network analysis.

This module provides the application shell: the PyPopARTApp class wires
the page layout (gui.layout) to the callback groups (gui.callbacks) and
runs the server. The web interface lets users upload sequence data,
configure network algorithms, visualize results, and export outputs.
"""

import logging
import tempfile

import click
import dash
import dash_bootstrap_components as dbc

from .callbacks import register_all
from .layout import build_layout


def _background_manager(cache_dir: str):
    """
    Build the manager that runs background callbacks.

    Parameters
    ----------
    cache_dir : str
        Directory for the diskcache store that carries job state and
        results between the server and the worker processes.

    Returns
    -------
    dash.DiskcacheManager
        Manager for the app's ``background=True`` callbacks.
    """
    from dash import DiskcacheManager
    import diskcache

    return DiskcacheManager(diskcache.Cache(cache_dir))


class PyPopARTApp:
    """
    Main PyPopART Dash application class.

    Provides a web-based interface for haplotype network analysis.

    Parameters
    ----------
    debug : bool, default=False
        Enable debug mode for development.
    port : int, default=8050
        Port number for the web server.
    """

    def __init__(self, debug: bool = False, port: int = 8050):
        """
        Initialize PyPopART Dash application.

        Parameters
        ----------
        debug : bool, default=False
            Enable debug mode for development.
        port : int, default=8050
            Port number for the web server.
        """
        self.debug = debug
        self.port = port

        # Configure logging
        log_level = logging.DEBUG if debug else logging.INFO
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        )
        self.logger = logging.getLogger(__name__)

        # Network construction and layout run as background callbacks so
        # a stale job can be cancelled when the user changes settings and
        # starts a fresh one. Diskcache hands results between processes.
        self.cache_dir = tempfile.mkdtemp(prefix='pypopart-jobs-')
        self.app = dash.Dash(
            __name__,
            external_stylesheets=[dbc.themes.BOOTSTRAP],
            suppress_callback_exceptions=True,
            background_callback_manager=_background_manager(self.cache_dir),
        )
        self.app.title = 'PyPopART - Haplotype Network Analysis'

        build_layout(self.app)
        register_all(self.app, self.logger)

    def run(self) -> None:
        """Run the Dash application."""
        # Use run() for Dash 2.0+, which replaced run_server()
        self.app.run(debug=self.debug, port=self.port)


def create_app(debug: bool = False, port: int = 8050) -> PyPopARTApp:
    """
    Create a configured PyPopART GUI application.

    Parameters
    ----------
    debug : bool, default=False
        Enable debug mode for development.
    port : int, default=8050
        Port number for the web server.

    Returns
    -------
    PyPopARTApp
        The assembled application.
    """
    return PyPopARTApp(debug=debug, port=port)


@click.command(
    name='pypopart-gui',
    help='Launch the PyPopART GUI application.\n\nOnce started, open your '
    'browser to http://localhost:PORT. Press Ctrl+C to stop the server.',
)
@click.option('--debug', is_flag=True, help='Enable debug mode for development.')
@click.option(
    '--port',
    type=int,
    default=8050,
    show_default=True,
    help='Port number for web server.',
)
def main(debug: bool = False, port: int = 8050) -> None:
    """
    Launch the PyPopART GUI application.

    Once started, open your browser to http://localhost:PORT.
    Press Ctrl+C to stop the server.

    Parameters
    ----------
    debug : bool, default=False
        Enable debug mode.
    port : int, default=8050
        Port number for web server.
    """
    print('=' * 60)
    print('PyPopART GUI - Haplotype Network Analysis')
    print('=' * 60)
    print(f'\n🚀 Starting web server on http://localhost:{port}')
    if debug:
        print('⚠️  Debug mode enabled')
    print('\n📖 Quick Start:')
    print('   1. Upload your sequence alignment (FASTA, NEXUS, or PHYLIP)')
    print('   2. Optionally upload metadata (CSV with population/location data)')
    print('   3. Choose a network algorithm and click "Compute Network"')
    print('   4. Customize the layout and export your results')
    print('\n⚠️  To stop the server, press Ctrl+C')
    print('=' * 60)
    print()

    create_app(debug=debug, port=port).run()


if __name__ == '__main__':
    main()
