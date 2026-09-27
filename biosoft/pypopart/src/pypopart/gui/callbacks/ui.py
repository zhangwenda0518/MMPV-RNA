"""Chrome-level UI behaviour: collapsing and restoring the control panel."""

from dash import Input, Output, State


def register(app, logger) -> None:
    """
    Register UI chrome callbacks on the Dash app.

    Parameters
    ----------
    app : dash.Dash
        The Dash application.
    logger : logging.Logger
        Application logger.
    """
    # Clientside on purpose. The sidebar carries CSS 'resize: horizontal',
    # so a drag writes an *inline* width straight onto the DOM node. A
    # server callback returning a new style dict would make React re-apply
    # the prop and wipe that dragged width on every toggle. Touching
    # el.style here leaves Dash's prop untouched, so the width survives.
    app.clientside_callback(
        """
        function(nClicks, collapsed) {
            if (!nClicks) {
                return window.dash_clientside.no_update;
            }

            const el = document.getElementById('sidebar-panel');
            if (!el) {
                return window.dash_clientside.no_update;
            }

            const nowCollapsed = !collapsed;
            if (nowCollapsed) {
                // Stash whatever width the user dragged to, then hide.
                window.pypopartSidebarWidth =
                    el.style.width || el.getBoundingClientRect().width + 'px';
                el.style.display = 'none';
            } else {
                el.style.display = '';
                if (window.pypopartSidebarWidth) {
                    el.style.width = window.pypopartSidebarWidth;
                }
            }

            // The debounced window-resize handler installed alongside the
            // network already calls cy.resize()/cy.fit(), so reuse it.
            setTimeout(function() {
                window.dispatchEvent(new Event('resize'));
            }, 0);

            // That handler only exists once a network has been drawn, and
            // it is debounced, so also nudge the instance directly.
            try {
                const cy = document.getElementById('network-graph')._cyreg.cy;
                if (cy) {
                    cy.resize();
                    cy.fit(null, 50);
                }
            } catch (e) {
                // No network yet; nothing to re-fit.
            }

            return [
                nowCollapsed,
                nowCollapsed ? '\\u00bb' : '\\u00ab',
                nowCollapsed ? 'Show the control panel' : 'Hide the control panel',
            ];
        }
        """,
        [
            Output('sidebar-collapsed', 'data'),
            Output('sidebar-toggle', 'children'),
            Output('sidebar-toggle', 'title'),
        ],
        Input('sidebar-toggle', 'n_clicks'),
        State('sidebar-collapsed', 'data'),
        prevent_initial_call=True,
    )

    # Box select: with panning off, a plain drag on the canvas draws a
    # selection box instead of moving the view. Cytoscape handles the
    # rest natively (Cmd/Ctrl+click adds, background click clears, and
    # dragging one selected node drags them all).
    app.clientside_callback(
        """
        function(boxSelect) {
            return !boxSelect;
        }
        """,
        Output('network-graph', 'userPanningEnabled'),
        Input('box-select-toggle', 'value'),
    )

    # Publishes the grid settings onto `window` so the Cytoscape event
    # handlers below can read them without a round trip on every drag.
    app.clientside_callback(
        """
        function(enabled, size) {
            window.pypopartGrid = {
                enabled: !!enabled,
                size: size || 0,
            };
            if (window.pypopartDrawGrid) {
                window.pypopartDrawGrid();
            }
            return window.dash_clientside.no_update;
        }
        """,
        Output('grid-size', 'className'),
        [
            Input('snap-to-grid-toggle', 'value'),
            Input('grid-size', 'value'),
        ],
    )

    # Snapping on drop, plus the grid drawn behind the network. The grid is
    # a CSS background on the Cytoscape container, so it has to be redrawn
    # on pan and zoom or it would only line up with the nodes at zoom 1.
    app.clientside_callback(
        """
        function() {
            if (window.pypopartGridSetup) {
                return window.dash_clientside.no_update;
            }
            window.pypopartGridSetup = true;

            setTimeout(function() {
                try {
                    const cy = document.getElementById('network-graph')._cyreg.cy;
                    if (!cy) { return; }
                    const container = cy.container();

                    function settings() {
                        return window.pypopartGrid || {enabled: false, size: 0};
                    }

                    window.pypopartDrawGrid = function() {
                        const grid = settings();
                        if (!grid.enabled || !grid.size) {
                            container.style.backgroundImage = '';
                            return;
                        }
                        const step = grid.size * cy.zoom();
                        const pan = cy.pan();
                        const line = 'rgba(0, 48, 73, 0.12)';
                        container.style.backgroundImage =
                            'linear-gradient(to right, ' + line + ' 1px, transparent 1px),' +
                            'linear-gradient(to bottom, ' + line + ' 1px, transparent 1px)';
                        container.style.backgroundSize = step + 'px ' + step + 'px';
                        container.style.backgroundPosition =
                            pan.x + 'px ' + pan.y + 'px';
                    };

                    cy.on('pan zoom resize', window.pypopartDrawGrid);

                    // Mirrors resolve_grid_collisions in
                    // layout/algorithms.py: same cost model, same
                    // tie-break, so what the user sees on drop matches
                    // what the server stores. Run here as well because
                    // a drag no longer round-trips through a redraw.
                    // resolve_edge_overlaps is deliberately not
                    // mirrored: a node dropped onto a line was put
                    // there on purpose, so a drop only de-collides.
                    const OCCUPIED_PENALTY = 3;
                    const MAX_SEARCH = 12;
                    const STEPS = [[1,0], [-1,0], [0,1], [0,-1]];

                    function freeCell(start, taken, away) {
                        const key = c => c[0] + ',' + c[1];
                        const seen = new Set([key(start)]);
                        let queue = [{cost: 0, cell: start}];
                        while (queue.length) {
                            queue.sort(function(a, b) {
                                if (a.cost !== b.cost) return a.cost - b.cost;
                                return b.align - a.align;
                            });
                            const node = queue.shift();
                            const cell = node.cell;
                            if (!taken.has(key(cell)) &&
                                key(cell) !== key(start)) {
                                return cell;
                            }
                            if (node.cost >= MAX_SEARCH) { continue; }
                            for (const [sx, sy] of STEPS) {
                                const next = [cell[0] + sx, cell[1] + sy];
                                if (seen.has(key(next))) { continue; }
                                seen.add(key(next));
                                const dx = next[0] - start[0];
                                const dy = next[1] - start[1];
                                const len = Math.hypot(dx, dy) || 1;
                                queue.push({
                                    cost: node.cost + 1 +
                                        (taken.has(key(next)) ? OCCUPIED_PENALTY : 0),
                                    align: (dx / len) * away[0] + (dy / len) * away[1],
                                    cell: next,
                                });
                            }
                        }
                        return start;
                    }

                    cy.on('dragfree', 'node', function(evt) {
                        const grid = settings();
                        if (!grid.enabled || !grid.size) { return; }
                        const node = evt.target;
                        const pos = node.position();
                        const cell = [
                            Math.round(pos.x / grid.size),
                            Math.round(pos.y / grid.size),
                        ];

                        // Every other node's cell is spoken for.
                        const taken = new Set();
                        cy.nodes().forEach(function(other) {
                            if (other.id() === node.id()) { return; }
                            const p = other.position();
                            taken.add(
                                Math.round(p.x / grid.size) + ',' +
                                Math.round(p.y / grid.size)
                            );
                        });

                        let target = cell;
                        if (taken.has(cell[0] + ',' + cell[1])) {
                            // Prefer moving away from the closest
                            // connected neighbour, so an edge is not
                            // folded back over itself.
                            let away = [0, 0];
                            let best = Infinity;
                            node.neighborhood('node').forEach(function(nb) {
                                const p = nb.position();
                                const nc = [
                                    Math.round(p.x / grid.size),
                                    Math.round(p.y / grid.size),
                                ];
                                const d = (nc[0] - cell[0]) ** 2 +
                                          (nc[1] - cell[1]) ** 2;
                                if (d < best) {
                                    best = d;
                                    const dx = cell[0] - nc[0];
                                    const dy = cell[1] - nc[1];
                                    const len = Math.hypot(dx, dy);
                                    away = len ? [dx / len, dy / len] : [0, 0];
                                }
                            });
                            target = freeCell(cell, taken, away);
                        }

                        node.position({
                            x: target[0] * grid.size,
                            y: target[1] * grid.size,
                        });
                    });

                    window.pypopartDrawGrid();
                } catch (e) {
                    console.log('Error setting up grid:', e);
                }
            }, 500);

            return window.dash_clientside.no_update;
        }
        """,
        Output('snap-to-grid-toggle', 'className'),
        Input('network-graph', 'elements'),
        prevent_initial_call=True,
    )

    # Refit the mutation tick marks to the edge that is actually visible
    # between two node discs. The server seeds each edge's tick_font from
    # its own idea of node size (fit_tick_font_size in
    # visualization/cytoscape_plot.py); this recomputes it from rendered
    # widths, so the node-size slider and every drag are honoured without
    # a round trip. Same constants and formula as the Python.
    app.clientside_callback(
        """
        function() {
            if (window.pypopartTickSetup) {
                return window.dash_clientside.no_update;
            }
            window.pypopartTickSetup = true;

            setTimeout(function() {
                try {
                    const cy = document.getElementById('network-graph')._cyreg.cy;
                    if (!cy) { return; }

                    const TICK_FONT_SIZE = 14;
                    const MIN_TICK_FONT_SIZE = 7;
                    const TICK_GLYPH_ADVANCE = 0.6;
                    const TICK_FIT_FRACTION = 0.7;

                    function fit(length, radiusA, radiusB, count) {
                        if (count <= 0) { return 0; }
                        const visible = length - radiusA - radiusB;
                        if (visible <= 0) { return 0; }
                        const glyphs = 2 * count - 1;
                        const fitted = visible * TICK_FIT_FRACTION /
                            (glyphs * TICK_GLYPH_ADVANCE);
                        if (fitted < MIN_TICK_FONT_SIZE) { return 0; }
                        return Math.round(Math.min(fitted, TICK_FONT_SIZE) * 100) / 100;
                    }

                    window.pypopartFitTicks = function() {
                        cy.batch(function() {
                            cy.edges().forEach(function(edge) {
                                if (!edge.data('ticks')) { return; }
                                const a = edge.source();
                                const b = edge.target();
                                const pa = a.position();
                                const pb = b.position();
                                const length = Math.hypot(pb.x - pa.x, pb.y - pa.y);
                                const size = fit(
                                    length, a.width() / 2, b.width() / 2,
                                    Number(edge.data('distance')) || 0
                                );
                                if (edge.data('tick_font') !== size) {
                                    edge.data('tick_font', size);
                                }
                            });
                        });
                    };

                    let pending = null;
                    function refit() {
                        if (pending) { clearTimeout(pending); }
                        pending = setTimeout(function() {
                            pending = null;
                            window.pypopartFitTicks();
                        }, 50);
                    }

                    cy.on('dragfree', 'node', refit);
                    cy.on('add', refit);
                    cy.on('style', 'node', refit);
                    refit();
                } catch (e) {
                    console.log('Error setting up tick fitting:', e);
                }
            }, 500);

            return window.dash_clientside.no_update;
        }
        """,
        Output('edge-tick-toggle', 'className'),
        Input('network-graph', 'elements'),
        prevent_initial_call=True,
    )
