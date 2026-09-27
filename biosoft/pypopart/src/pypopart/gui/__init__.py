"""PyPopART GUI module using Dash."""

try:
    from pypopart.gui.app import PyPopARTApp, main
except ImportError as e:
    raise ImportError(
        'The PyPopART GUI requires Dash. '
        "Install the GUI extra with: pip install 'pypopart[gui]'"
    ) from e

__all__ = ['PyPopARTApp', 'main']
