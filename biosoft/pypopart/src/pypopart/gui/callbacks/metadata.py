"""
Draft-until-commit editing of the metadata table.

Edits typed into the metadata table are held as a draft and only reach
the figure when the user clicks Compute Network, so a half-finished set
of population assignments never repaints the network.
"""

from typing import Dict, List, Optional, Tuple

from dash import ALL, Input, Output, State, callback_context, html
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from ..metadata_edit import (
    apply_population_color,
    build_metadata_rows,
    diff_metadata_rows,
    rows_to_metadata_store,
)


def register(app, logger) -> None:
    """
    Register metadata editing callbacks on the Dash app.

    Parameters
    ----------
    app : dash.Dash
        The Dash application.
    logger : logging.Logger
        Application logger.
    """

    @app.callback(
        Output('metadata-draft-store', 'data'),
        Input('metadata-table', 'data'),
        [State('alignment-store', 'data'), State('metadata-store', 'data')],
        prevent_initial_call=True,
    )
    def capture_metadata_edits(
        edited_rows: Optional[List[Dict]],
        alignment_data: Optional[Dict],
        metadata_data: Optional[Dict],
    ) -> Optional[Dict]:
        """
        Hold table edits as a draft, diffed against the committed rows.

        Parameters
        ----------
        edited_rows : list of dict, optional
            Current contents of the metadata table.
        alignment_data : dict, optional
            Serialized alignment from the alignment store.
        metadata_data : dict, optional
            Committed metadata store contents.

        Returns
        -------
        dict or None
            The draft, or ``None`` when the table matches what is applied.
        """
        if not edited_rows or not alignment_data:
            raise PreventUpdate

        alignment_ids = {seq['id'] for seq in alignment_data['sequences']}
        committed = build_metadata_rows(alignment_ids, metadata_data)
        changes = diff_metadata_rows(committed, edited_rows)

        if not changes:
            return None
        return {'rows': edited_rows, 'edited': [list(change) for change in changes]}

    @app.callback(
        [
            Output('metadata-draft-badge', 'children'),
            Output('metadata-draft-badge-sidebar', 'children'),
            Output('metadata-edit-feedback', 'children'),
        ],
        Input('metadata-draft-store', 'data'),
    )
    def show_draft_badge(
        draft: Optional[Dict],
    ) -> Tuple[object, object, object]:
        """
        Surface pending edits in the Metadata tab and beside the button.

        Parameters
        ----------
        draft : dict, optional
            The uncommitted draft, if any.

        Returns
        -------
        tuple
            Badge for the tab, badge for the sidebar, and a validation
            message for any unusable cell values.
        """
        if not draft or not draft.get('edited'):
            return html.Div(), html.Div(), html.Div()

        count = len(draft['edited'])
        text = f'{count} unapplied edit{"" if count == 1 else "s"}'
        badge = dbc.Badge(text, color='warning', className='me-1')
        sidebar = dbc.Badge(
            f'{text} — recompute to apply',
            color='warning',
            className='d-block text-wrap',
        )

        # Dry run so bad coordinates are flagged before the recompute.
        _, warnings = rows_to_metadata_store(draft['rows'])
        feedback = html.Div()
        if warnings:
            feedback = dbc.Alert(
                [html.Div(warning) for warning in warnings],
                color='warning',
                className='mt-2 mb-0',
            )

        return badge, sidebar, feedback

    @app.callback(
        [
            Output('metadata-store', 'data', allow_duplicate=True),
            Output('metadata-commit-token', 'data'),
            Output('metadata-draft-store', 'data', allow_duplicate=True),
        ],
        Input('compute-button', 'n_clicks'),
        [State('metadata-draft-store', 'data'), State('metadata-store', 'data')],
        prevent_initial_call=True,
    )
    def commit_metadata_edits(
        n_clicks: Optional[int],
        draft: Optional[Dict],
        metadata_data: Optional[Dict],
    ) -> Tuple[object, int, None]:
        """
        Apply any pending edits, then release the network computation.

        This callback is the only producer of ``metadata-commit-token``
        and ``compute_network`` is its only consumer, which serialises the
        two writes: the metadata store is updated *before* the network is
        rebuilt. Firing both off the button directly would race, and the
        network could be built from pre-edit metadata.

        Parameters
        ----------
        n_clicks : int, optional
            Compute button click count.
        draft : dict, optional
            The uncommitted draft, if any.
        metadata_data : dict, optional
            Current metadata store contents.

        Returns
        -------
        tuple
            New metadata store contents, the commit token, and a cleared
            draft.
        """
        from dash import no_update

        if not n_clicks:
            raise PreventUpdate

        if not draft or not draft.get('rows'):
            # Nothing to apply, but the computation still has to run.
            return no_update, n_clicks, None

        store, warnings = rows_to_metadata_store(draft['rows'], metadata_data)
        logger.info(
            f'Applied {len(draft.get("edited", []))} metadata edit(s)'
            + (f' with {len(warnings)} warning(s)' if warnings else '')
        )
        return store, n_clicks, None

    @app.callback(
        Output('metadata-draft-store', 'data', allow_duplicate=True),
        Input('metadata-store', 'data'),
        prevent_initial_call=True,
    )
    def clear_stale_draft(metadata_data: Optional[Dict]) -> None:
        """
        Drop the draft whenever the committed metadata changes.

        Without this, uploading a new CSV mid-draft would leave stale
        values ready to be applied on the next compute.

        Parameters
        ----------
        metadata_data : dict, optional
            New metadata store contents.

        Returns
        -------
        None
            Always clears the draft.
        """
        return None

    @app.callback(
        Output('metadata-table', 'data', allow_duplicate=True),
        Input({'type': 'pop-color', 'pop': ALL}, 'value'),
        State('metadata-table', 'data'),
        prevent_initial_call=True,
    )
    def apply_population_swatch(
        colors: List[Optional[str]],
        rows: Optional[List[Dict]],
    ) -> List[Dict]:
        """
        Recolour every row of a population when its swatch changes.

        Writing back into the table means the change flows through the same
        draft-and-commit path as a typed edit, so the figure still only
        moves on Compute Network.

        Parameters
        ----------
        colors : list of str
            Current value of every population swatch.
        rows : list of dict, optional
            Current contents of the metadata table.

        Returns
        -------
        list of dict
            The table rows with the new colour applied.
        """
        if not rows or not colors:
            raise PreventUpdate

        triggered = callback_context.triggered_id
        if not isinstance(triggered, dict) or triggered.get('type') != 'pop-color':
            raise PreventUpdate

        population = triggered.get('pop')
        # The swatch ids and their values share an order, so the changed
        # population's value can be read back positionally.
        new_color = None
        for control, value in zip(callback_context.inputs_list[0], colors):
            if control.get('id', {}).get('pop') == population:
                new_color = value
                break

        if not new_color:
            raise PreventUpdate

        updated = apply_population_color(rows, population, new_color)
        if updated == rows:
            raise PreventUpdate
        return updated
