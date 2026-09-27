"""
Shared helpers for user-facing status messages.

Transient confirmations go to a single toast in the bottom-right corner;
errors and validation warnings stay inline in the panel they belong to,
where they persist until the user acts on them. A traceback that fades
after four seconds is worse than no traceback at all.
"""

from typing import Tuple

from dash import Output, html

#: How long a success toast stays on screen, in milliseconds.
TOAST_DURATION_MS = 4000


def toast(*message, header: str = 'Success') -> Tuple[object, str, bool]:
    """
    Build the outputs that raise the shared toast.

    Parameters
    ----------
    *message : object
        Children for the toast body: strings, Dash components, or both.
    header : str, default='Success'
        Toast title.

    Returns
    -------
    tuple
        ``(children, header, is_open)`` for the ``app-toast`` outputs.
    """
    children = list(message)
    if len(children) == 1 and isinstance(children[0], (list, tuple)):
        children = list(children[0])
    return html.Div(children), header, True


def no_toast() -> Tuple[object, str, bool]:
    """
    Build the outputs that leave the toast closed.

    Returns
    -------
    tuple
        ``(children, header, is_open)`` with the toast shut.
    """
    return html.Div(), '', False


def task_running(label: str, *extra) -> list:
    """
    Build a ``running`` spec that shows a task in the title bar.

    Parameters
    ----------
    label : str
        What the app is doing, shown beside the spinner while the
        callback runs.
    *extra : tuple
        Further ``(Output, running_value, idle_value)`` entries, such as
        disabling the button that started the task.

    Returns
    -------
    list
        Entries for a callback's ``running=`` argument.
    """
    return [
        (Output('task-indicator', 'style'), {'display': 'flex'}, {'display': 'none'}),
        (Output('task-label', 'children'), label, ''),
        *extra,
    ]
