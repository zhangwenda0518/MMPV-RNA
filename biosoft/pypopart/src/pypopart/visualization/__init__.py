"""
Visualization module for PyPopART.

Provides functions for creating static and interactive visualizations
of haplotype networks.

The matplotlib- and plotly-based plotters require the optional ``viz``
extra (``pip install 'pypopart[viz]'``) and are imported lazily so the
package remains importable without them.
"""

from .cytoscape_plot import (
    InteractiveCytoscapePlotter,
    create_cytoscape_network,
)

_LAZY_IMPORTS = {
    # Static plotting (matplotlib)
    'StaticNetworkPlotter': ('.static_plot', 'matplotlib'),
    'plot_network': ('.static_plot', 'matplotlib'),
    'create_publication_figure': ('.static_plot', 'matplotlib'),
    # Interactive plotting (plotly)
    'InteractiveNetworkPlotter': ('.interactive_plot', 'plotly'),
    'plot_interactive_network': ('.interactive_plot', 'plotly'),
    'create_interactive_figure': ('.interactive_plot', 'plotly'),
}

__all__ = [
    'StaticNetworkPlotter',
    'plot_network',
    'create_publication_figure',
    'InteractiveNetworkPlotter',
    'plot_interactive_network',
    'create_interactive_figure',
    'InteractiveCytoscapePlotter',
    'create_cytoscape_network',
]


def __getattr__(name: str):
    """
    Lazily import optional plotters, with an install hint when missing.

    Parameters
    ----------
    name : str
        Attribute name being looked up on the package.

    Returns
    -------
    Any
        The requested plotter class or function.
    """
    if name in _LAZY_IMPORTS:
        module_name, backend = _LAZY_IMPORTS[name]
        from importlib import import_module

        try:
            module = import_module(module_name, __name__)
        except ImportError as e:
            raise ImportError(
                f'{name} requires {backend}. '
                "Install the visualization extra with: pip install 'pypopart[viz]'"
            ) from e
        return getattr(module, name)
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
