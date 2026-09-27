"""Dash callback registration for the PyPopART GUI, grouped by concern."""

from . import display, export, layout, metadata, network, ui, upload


def register_all(app, logger) -> None:
    """
    Register every callback group on the Dash app.

    Parameters
    ----------
    app : dash.Dash
        The Dash application.
    logger : logging.Logger
        Application logger shared by all callbacks.
    """
    upload.register(app, logger)
    network.register(app, logger)
    layout.register(app, logger)
    display.register(app, logger)
    metadata.register(app, logger)
    export.register(app, logger)
    ui.register(app, logger)


__all__ = [
    'register_all',
    'upload',
    'network',
    'layout',
    'display',
    'metadata',
    'export',
    'ui',
]
