/*
 * Sidebar resize grip.
 *
 * Lives in assets/ rather than in a Dash clientside callback because it
 * uses event delegation on `document` and therefore needs no element to
 * exist at load time -- a clientside callback would need an Input that has
 * already fired, and the grip must work before any network is drawn.
 *
 * The dragged width is written straight onto the panel's inline style, the
 * same invariant callbacks/ui.py relies on: returning a `style` dict from a
 * callback makes React re-apply the prop and wipe the user's width.
 */
(function () {
  "use strict";

  var MIN_WIDTH = 250;
  var MAX_WIDTH = 600;

  var dragging = false;
  var startX = 0;
  var startWidth = 0;

  function sidebar() {
    return document.getElementById("sidebar-panel");
  }

  function refitNetwork() {
    try {
      var cy = document.getElementById("network-graph")._cyreg.cy;
      if (cy) {
        cy.resize();
        cy.fit(null, 50);
      }
    } catch (e) {
      // No network drawn yet; nothing to re-fit.
    }
  }

  document.addEventListener("mousedown", function (event) {
    var grip = event.target.closest && event.target.closest("#sidebar-resizer");
    if (!grip) {
      return;
    }
    // Let the collapse button keep its own click.
    if (event.target.closest("#sidebar-toggle")) {
      return;
    }

    var panel = sidebar();
    if (!panel) {
      return;
    }

    dragging = true;
    startX = event.clientX;
    startWidth = panel.getBoundingClientRect().width;
    document.body.style.userSelect = "none";
    document.body.style.cursor = "col-resize";
    event.preventDefault();
  });

  document.addEventListener("mousemove", function (event) {
    if (!dragging) {
      return;
    }
    var panel = sidebar();
    if (!panel) {
      return;
    }

    var width = startWidth + (event.clientX - startX);
    width = Math.min(MAX_WIDTH, Math.max(MIN_WIDTH, width));
    panel.style.width = width + "px";
  });

  document.addEventListener("mouseup", function () {
    if (!dragging) {
      return;
    }
    dragging = false;
    document.body.style.userSelect = "";
    document.body.style.cursor = "";

    // The debounced window-resize handler installed alongside the network
    // never fired for the old native corner drag, so the graph kept its
    // stale width. Nudge it directly.
    refitNetwork();
  });
})();
