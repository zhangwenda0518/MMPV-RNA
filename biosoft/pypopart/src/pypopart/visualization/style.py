"""
Shared node styling for the static, interactive, and Cytoscape plotters.

Single source of truth for the app's palette, for the colour-priority rule
(median colour > explicit node colour > dominant-population colour >
default) and for generating distinct population colours.

The palette constants here are mirrored in ``gui/assets/theme.css``, which
cannot import Python. ``tests/unit/test_style.py`` pins the two together so
they cannot drift apart silently.
"""

import colorsys
from typing import Dict, List, Optional

#: Bauhaus-flavoured pop art palette used across the app chrome, the
#: network canvas and the static exports.
POP_RED = '#d62828'
POP_NAVY = '#003049'
POP_AMBER = '#fcbf49'
POP_BONE = '#eae2b7'
POP_INK = '#1d1d1b'
POP_PAPER = '#ffffff'
POP_TEAL = '#2a9d8f'

#: Bauhaus accent shades, used for the wordmark and app chrome.
POP_ORANGE = '#e8871e'
POP_GREEN = '#2e7d5b'
POP_GOLD = '#d9a521'
POP_PINK = '#d9678a'
POP_BLUE = '#2d6a9f'

#: Off-white ground for the app. Screens are not paper, and a page of
#: pure white behind a figure is glare; exports still get true white.
POP_PARCHMENT = '#f6f2e9'

#: Every palette colour, in the order the CSS declares them.
PALETTE = (
    POP_RED,
    POP_NAVY,
    POP_AMBER,
    POP_BONE,
    POP_INK,
    POP_PAPER,
    POP_TEAL,
    POP_ORANGE,
    POP_GREEN,
    POP_GOLD,
    POP_PINK,
    POP_BLUE,
    POP_PARCHMENT,
)

#: Lowest Rec. 601 luma a node fill may have. Node labels are ink with a
#: thin paper outline, so a dark fill makes them hard to read.
MIN_NODE_LUMA = 120.0

#: Colours assigned to populations, in order. All are light enough to
#: carry an ink label, which is why the navy and ink of the base palette
#: are not among them.
POP_ART_PALETTE = (
    '#f4534d',  # vermilion
    POP_AMBER,
    POP_TEAL,
    '#4cc9f0',  # sky
    '#f77f00',  # orange
    '#90be6d',  # leaf
    '#c77dff',  # lilac
    '#e5989b',  # rose
    '#8ecae6',  # pale blue
    '#f2cc8f',  # sand
    '#52b788',  # emerald
    '#ef476f',  # raspberry
    '#b5e48c',  # lime
    '#d4a5c9',  # orchid
    '#9db4c0',  # slate
    '#ffb703',  # gold
)

#: Default node colour used by every plotter.
DEFAULT_NODE_COLOR = '#8ecae6'

#: Default colour for inferred median/intermediate vertices.
DEFAULT_MEDIAN_COLOR = '#c9c9c4'


def luma(hex_color: str) -> float:
    """
    Measure a colour's perceived brightness.

    Rec. 601 weighting, the usual quick approximation for deciding
    whether a fill is light or dark.

    Parameters
    ----------
    hex_color : str
        Colour as ``#rrggbb``.

    Returns
    -------
    float
        Luma from 0 (black) to 255 (white). Unparseable input returns 0.0,
        which reads as 'dark' and so fails any brightness check.
    """
    value = hex_color.lstrip('#')
    if len(value) != 6:
        return 0.0
    try:
        r, g, b = (int(value[i : i + 2], 16) for i in (0, 2, 4))
    except ValueError:
        return 0.0
    return 0.299 * r + 0.587 * g + 0.114 * b


def node_color(
    node: str,
    haplotype,
    is_median: bool,
    node_color_map: Optional[Dict[str, str]] = None,
    population_colors: Optional[Dict[str, str]] = None,
    median_vector_color: str = DEFAULT_MEDIAN_COLOR,
    default: str = DEFAULT_NODE_COLOR,
) -> str:
    """
    Resolve a node's colour with the shared priority rule.

    Parameters
    ----------
    node : str
        Node ID.
    haplotype : Haplotype or None
        The node's haplotype, if any.
    is_median : bool
        Whether the node is an inferred median vector.
    node_color_map : dict, optional
        Explicit per-node colour overrides.
    population_colors : dict, optional
        Colour per population name.
    median_vector_color : str, default=DEFAULT_MEDIAN_COLOR
        Colour for median vectors.
    default : str, default=DEFAULT_NODE_COLOR
        Fallback colour.

    Returns
    -------
    str
        The resolved colour.
    """
    if is_median:
        return median_vector_color
    if node_color_map and node in node_color_map:
        return node_color_map[node]
    if population_colors and haplotype is not None:
        pop_counts = haplotype.get_frequency_by_population()
        if pop_counts:
            dominant = max(pop_counts.items(), key=lambda item: item[1])[0]
            return population_colors.get(dominant, default)
    return default


def generate_population_colors(populations: List[str]) -> Dict[str, str]:
    """
    Assign a distinct colour to each population.

    Uses the hand-picked pop art palette while it lasts and falls back to
    evenly spaced HSV hues beyond it. The palette is not cycled: repeating
    it would hand two populations the same colour, and distinguishable
    colours matter more than staying on-palette. Every colour it can
    return clears :data:`MIN_NODE_LUMA`.

    Parameters
    ----------
    populations : list of str
        Population names.

    Returns
    -------
    dict
        Mapping of population name to hex colour. Colours are unique for
        any number of populations.
    """
    names = sorted(populations)
    n = len(names)

    if n <= len(POP_ART_PALETTE):
        return {pop: POP_ART_PALETTE[i] for i, pop in enumerate(names)}

    # Low saturation at full value, so even the blue end of the wheel
    # stays above MIN_NODE_LUMA. The old 0.7/0.9 produced fills around
    # luma 86 -- darker than anything in the palette.
    colors = {}
    for i, pop in enumerate(names):
        hue = i / n if n else 0.0
        r, g, b = colorsys.hsv_to_rgb(hue, 0.45, 1.0)
        colors[pop] = '#{:02x}{:02x}{:02x}'.format(
            int(r * 255), int(g * 255), int(b * 255)
        )
    return colors


def apply_pop_art_rcparams() -> Dict[str, object]:
    """
    Apply the app's palette and typography to matplotlib.

    Called by the static plotter so exported figures carry the same look
    as the on-screen network. Fonts fall back to whatever the system has
    when the bundled families are not installed, which matplotlib handles
    on its own.

    Returns
    -------
    dict
        The rcParams that were set, so callers can restore them.
    """
    import matplotlib as mpl
    from matplotlib import font_manager

    params = {
        'figure.facecolor': POP_PAPER,
        'axes.facecolor': POP_PAPER,
        'axes.edgecolor': POP_INK,
        'axes.labelcolor': POP_INK,
        'axes.titlecolor': POP_INK,
        'axes.titleweight': 'bold',
        'text.color': POP_INK,
        'xtick.color': POP_INK,
        'ytick.color': POP_INK,
        'legend.frameon': True,
        'legend.edgecolor': POP_INK,
        'legend.facecolor': POP_PAPER,
    }

    # The app's own faces ship as woff2 for the browser, which matplotlib
    # cannot load. Ask for them only if the system happens to have them
    # installed; otherwise say nothing and let matplotlib use its default,
    # rather than emitting a findfont warning per label.
    available = {font.name for font in font_manager.fontManager.ttflist}
    preferred = [name for name in ('Jost', 'Space Grotesk') if name in available]
    if preferred:
        params['font.family'] = [*preferred, 'DejaVu Sans', 'sans-serif']

    mpl.rcParams.update(params)
    return params
