"""
Interactive network visualization using Dash Cytoscape for PyPopART.

Provides Cytoscape-based interactive plotting for haplotype networks
with manual node repositioning, pie chart nodes, and legend support.
"""

import base64
import math
from typing import Dict, List, Optional, Tuple

import numpy as np

from ..core.graph import HaplotypeNetwork
from .style import (
    DEFAULT_MEDIAN_COLOR,
    DEFAULT_NODE_COLOR,
    POP_AMBER,
    POP_INK,
    POP_NAVY,
    POP_PAPER,
)

#: Type stack for network labels, matching the app chrome. Cytoscape falls
#: back through the list itself when a family is not available.
LABEL_FONT = 'Space Grotesk, Helvetica, Arial, sans-serif'

#: Monospace stack for the mutation tick marks.
TICK_FONT = 'JetBrains Mono, ui-monospace, monospace'

#: Hard ceiling on tick marks drawn for one edge. Beyond this the comb is
#: unreadable at any zoom, so the numeral is shown instead.
MAX_TICK_MARKS = 30

#: Default mutation count above which an edge shows a numeral, not ticks.
DEFAULT_TICK_THRESHOLD = 10

#: Font size, in px, of the tick label on an edge long enough to take it.
TICK_FONT_SIZE = 14.0

#: Smallest tick font, in px, before the strokes merge into a smear. An
#: edge that cannot fit its comb even at this size shows the numeral.
MIN_TICK_FONT_SIZE = 7.0

#: Advance width of one glyph of the tick font as a fraction of the font
#: size. JetBrains Mono, like most monospace faces, is 0.6 em wide.
TICK_GLYPH_ADVANCE = 0.6

#: Share of the visible edge (between the two node discs) the comb may
#: occupy, leaving clear line on either side so it reads as a marker.
TICK_FIT_FRACTION = 0.7


def format_edge_ticks(distance: int, max_ticks: int = MAX_TICK_MARKS) -> str:
    """
    Render a mutation count as tick marks for an edge label.

    One '|' glyph per mutation, space separated. Cytoscape.js has no
    letter-spacing property, so the spacing has to be part of the string.
    Drawn with 'text-rotation: autorotate' the glyph strokes sit
    perpendicular to the edge, which is the PopART convention.

    Parameters
    ----------
    distance : int
        Number of mutations separating the two haplotypes.
    max_ticks : int, default=MAX_TICK_MARKS
        Counts above this return an empty string, so the caller can fall
        back to a numeral.

    Returns
    -------
    str
        Space-separated tick marks, or '' when out of range.
    """
    if distance <= 0 or distance > max_ticks:
        return ''
    return ' '.join('|' * int(distance))


def fit_tick_font_size(
    edge_length: float, radius_a: float, radius_b: float, count: int
) -> float:
    """
    Choose a tick font size so the comb fits between two node discs.

    The ticks are a text label, so the only way to pack them closer on
    a short edge is to shrink the font. The comb is ``count`` pipes with
    a space between each, ``2 * count - 1`` glyphs wide, and may take up
    :data:`TICK_FIT_FRACTION` of the edge that is actually visible
    between the nodes.

    Parameters
    ----------
    edge_length : float
        Centre-to-centre length of the edge, in px.
    radius_a : float
        Radius of the disc at one end, in px.
    radius_b : float
        Radius of the disc at the other end, in px.
    count : int
        Number of tick marks.

    Returns
    -------
    float
        Font size in px, between :data:`MIN_TICK_FONT_SIZE` and
        :data:`TICK_FONT_SIZE`, or ``0.0`` when even the smallest comb
        would not fit, which tells the stylesheet to show the numeral.
    """
    if count <= 0:
        return 0.0
    visible = edge_length - radius_a - radius_b
    if visible <= 0:
        return 0.0
    glyphs = 2 * count - 1
    fitted = visible * TICK_FIT_FRACTION / (glyphs * TICK_GLYPH_ADVANCE)
    if fitted < MIN_TICK_FONT_SIZE:
        return 0.0
    return round(min(fitted, TICK_FONT_SIZE), 2)


def create_edge_tick_stylesheet(
    threshold: int = DEFAULT_TICK_THRESHOLD, show_ticks: bool = True
) -> List[Dict]:
    """
    Build the edge label rules for tick marks and numerals.

    Every selector is prefixed 'edge[distance' so callers can strip and
    replace the whole group in one pass.

    The tick font size comes from each edge's ``tick_font`` datum (see
    :func:`fit_tick_font_size`); an edge whose comb will not fit between
    its nodes carries ``0`` and shows the numeral instead.

    These rules must be appended *after* the base 'edge[label]' rule:
    every edge carries a 'label' key, so that rule matches all edges and
    would otherwise paint a white background pill behind the ticks
    (Cytoscape.js resolves conflicts by stylesheet order, later wins).

    Parameters
    ----------
    threshold : int, default=DEFAULT_TICK_THRESHOLD
        Edges up to this many mutations get ticks; longer ones get the
        numeral.
    show_ticks : bool, default=True
        When False, every edge shows the numeral.

    Returns
    -------
    List[Dict]
        Cytoscape stylesheet rules.
    """
    numeral_style = {
        'label': 'data(label)',
        'text-rotation': 'none',
        'font-family': LABEL_FONT,
        'font-size': '10px',
        'text-background-color': POP_PAPER,
        'text-background-opacity': 0.7,
        'text-background-padding': '3px',
        'color': POP_INK,
    }

    if not show_ticks:
        return [{'selector': 'edge[distance > 0]', 'style': numeral_style}]

    return [
        {
            'selector': f'edge[distance <= {threshold}][tick_font > 0]',
            'style': {
                'label': 'data(ticks)',
                'text-rotation': 'autorotate',
                'font-family': TICK_FONT,
                'font-size': 'data(tick_font)',
                'color': POP_INK,
                # no pill behind tick marks (overrides edge[label])
                'text-background-opacity': 0,
            },
        },
        {
            'selector': f'edge[distance <= {threshold}][tick_font = 0]',
            'style': numeral_style,
        },
        {
            'selector': f'edge[distance > {threshold}]',
            'style': numeral_style,
        },
    ]


def resolve_population_counts(
    hap, population_mapping: Optional[Dict]
) -> Dict[str, int]:
    """
    Count a haplotype's samples per population.

    Prefers the counts the haplotype carries itself and falls back to the
    GUI's sample-to-population mapping, which is the only source when the
    metadata arrived after the alignment was parsed.

    Parameters
    ----------
    hap : Haplotype
        Haplotype whose samples are being counted.
    population_mapping : Dict[str, str], optional
        Mapping of sample_id to population.

    Returns
    -------
    Dict[str, int]
        Population name to sample count. Samples with no population land
        under 'Unassigned'.
    """
    counts = hap.get_frequency_by_population()

    # 'Unassigned'-only counts mean the haplotype never saw the metadata,
    # so the mapping is the better source.
    stale = not counts or (len(counts) == 1 and 'Unassigned' in counts)
    if not (stale and population_mapping and hap.sample_ids):
        return counts

    counts = {}
    for sample_id in hap.sample_ids:
        pop = population_mapping.get(sample_id, 'Unassigned')
        counts[pop] = counts.get(pop, 0) + 1
    return counts


#: Sample IDs listed in a node tooltip before it starts counting the rest.
MAX_TOOLTIP_SAMPLES = 10


def build_node_tooltip(
    label: str,
    hap,
    is_median: bool,
    population_mapping: Optional[Dict] = None,
    max_samples: int = MAX_TOOLTIP_SAMPLES,
) -> Dict:
    """
    Build the hover tooltip payload for one node.

    Returns plain JSON-safe data rather than markup: it travels to the
    browser inside the node's Cytoscape data, and the tooltip is rendered
    there. Sample IDs come from user-supplied FASTA headers, so the
    renderer must insert them as text, never as HTML.

    Parameters
    ----------
    label : str
        Display label for the node, typically its H number.
    hap : Haplotype or None
        Haplotype at this node, if it has one.
    is_median : bool
        Whether the node is an inferred median vector.
    population_mapping : Dict[str, str], optional
        Mapping of sample_id to population.
    max_samples : int, default=MAX_TOOLTIP_SAMPLES
        Sample IDs to list before summarising the remainder.

    Returns
    -------
    Dict
        Keys: ``label``, ``kind`` ('median' or 'haplotype'), ``frequency``,
        ``populations`` (list of ``[name, count]`` pairs), ``samples`` and
        ``extra`` (how many sample IDs were not listed).
    """
    tooltip = {
        'label': label,
        'kind': 'median' if is_median or hap is None else 'haplotype',
        'frequency': 0,
        'populations': [],
        'samples': [],
        'extra': 0,
    }
    if tooltip['kind'] == 'median':
        return tooltip

    sample_ids = list(hap.sample_ids or [])
    counts = resolve_population_counts(hap, population_mapping)

    tooltip['frequency'] = hap.frequency
    tooltip['populations'] = [[pop, counts[pop]] for pop in sorted(counts)]
    tooltip['samples'] = sample_ids[:max_samples]
    tooltip['extra'] = max(0, len(sample_ids) - max_samples)
    return tooltip


class InteractiveCytoscapePlotter:
    """
    Interactive network plotter using Dash Cytoscape.

    Creates interactive visualizations of haplotype networks with
    manual node repositioning, pie chart nodes for population data,
    and customizable legends.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object to visualize.
    """

    def __init__(self, network: HaplotypeNetwork):
        """
        Initialize interactive Cytoscape plotter with a haplotype network.

        Parameters
        ----------
        network : HaplotypeNetwork
            HaplotypeNetwork object to visualize.
        """
        self.network = network
        self.elements = None
        self.stylesheet = None

    @staticmethod
    def generate_pie_chart_svg(pie_data: List[Dict]) -> str:
        """
        Generate SVG pie chart as Data URI for node background.

        Parameters
        ----------
        pie_data : List[Dict]
            List of dicts with 'percent' and 'color' keys for each segment.

        Returns
        -------
        str
            Data URI string with base64-encoded SVG.
        """
        size = 100
        center = size / 2
        radius = size / 2

        svg_parts = [
            f'<svg width="{size}" height="{size}" xmlns="http://www.w3.org/2000/svg">'
        ]

        # Start angle for pie slices (in radians)
        current_angle = -math.pi / 2  # Start at top (12 o'clock position)

        for segment in pie_data:
            percent = segment['percent']
            color = segment['color']

            # Calculate angle for this segment
            angle_size = (percent / 100) * 2 * math.pi

            # Calculate start and end points
            start_x = center + radius * math.cos(current_angle)
            start_y = center + radius * math.sin(current_angle)

            current_angle += angle_size

            end_x = center + radius * math.cos(current_angle)
            end_y = center + radius * math.sin(current_angle)

            # Use large-arc-flag if angle > 180 degrees
            large_arc = 1 if angle_size > math.pi else 0

            # Create pie slice path
            if percent == 100:
                # Full circle
                svg_parts.append(
                    f'<circle cx="{center}" cy="{center}" r="{radius}" fill="{color}"/>'
                )
            else:
                # Pie slice
                path = (
                    f'M {center},{center} '
                    f'L {start_x},{start_y} '
                    f'A {radius},{radius} 0 {large_arc},1 {end_x},{end_y} '
                    f'Z'
                )
                svg_parts.append(f'<path d="{path}" fill="{color}"/>')

        svg_parts.append('</svg>')
        svg_str = ''.join(svg_parts)

        # Encode as Data URI
        encoded = base64.b64encode(svg_str.encode('utf-8')).decode('utf-8')
        return f'data:image/svg+xml;base64,{encoded}'

    def create_elements(
        self,
        layout: Optional[Dict[str, Tuple[float, float]]] = None,
        node_size_scale: float = 20.0,
        population_colors: Optional[Dict[str, str]] = None,
        population_mapping: Optional[Dict[str, str]] = None,
        show_labels: bool = True,
        show_edge_labels: bool = True,
        median_vector_color: str = DEFAULT_MEDIAN_COLOR,
        node_labels: Optional[Dict[str, str]] = None,
        max_tick_marks: int = MAX_TICK_MARKS,
    ) -> List[Dict]:
        """
        Create Cytoscape elements from network data.

        Parameters
        ----------
        layout : Dict[str, Tuple[float, float]], optional
            Pre-computed node positions {node_id: (x, y)}.
        node_size_scale : float, default=20.0
            Scaling factor for node sizes.
        population_colors : Dict[str, str], optional
            Color mapping for populations {pop_name: color}.
        population_mapping : Dict[str, str], optional
            Mapping of sample_id to population {sample_id: population}.
        show_labels : bool, default=True
            Whether to show node labels.
        show_edge_labels : bool, default=True
            Whether to show edge labels with mutation counts.
        median_vector_color : str, default=DEFAULT_MEDIAN_COLOR
            Color for median vector nodes.
        node_labels : Dict[str, str], optional
            Custom labels for nodes {node_id: label}.
        max_tick_marks : int, default=MAX_TICK_MARKS
            Cap on tick marks emitted per edge.

        Returns
        -------
        List[Dict]
            List of Cytoscape element dictionaries.
        """
        elements = []
        graph = self.network._graph

        # Create nodes
        for node in graph.nodes():
            hap = self.network.get_haplotype(node)
            is_median = self.network.is_median_vector(node)

            # Get position
            pos = layout.get(node, (0, 0)) if layout else (0, 0)

            # Calculate size
            if is_median:
                size = node_size_scale * 0.8
            elif hap:
                # Use sqrt of frequency, but ensure minimum size
                size = max(
                    node_size_scale * 0.5, node_size_scale * np.sqrt(hap.frequency)
                )
            else:
                size = node_size_scale * 0.5

            # Determine label
            if show_labels:
                if node_labels and node in node_labels:
                    label = node_labels[node]
                else:
                    label = node
            else:
                label = ''

            # Build node data
            node_data = {
                'id': node,
                'label': label,
                'size': size,
                'is_median': is_median,
            }

            # Add population pie chart data if available
            if not is_median and hap and population_colors:
                pop_counts = resolve_population_counts(hap, population_mapping)

                # Pie charts only make sense for mixed-population nodes;
                # single-population nodes get that population's solid colour.
                if len(pop_counts) > 1:
                    # Prepare pie chart display for all nodes with populations
                    total = sum(pop_counts.values())
                    pie_data = []
                    pie_colors = []
                    pie_sizes = []

                    for pop, count in sorted(pop_counts.items()):
                        if count > 0:
                            percent = (count / total) * 100
                            # Unassigned samples take the median grey, not a population colour
                            if pop == 'Unassigned':
                                color = DEFAULT_MEDIAN_COLOR  # unassigned samples
                            else:
                                color = population_colors.get(pop, '#cccccc')

                            pie_data.append(
                                {
                                    'population': pop,
                                    'value': count,
                                    'percent': percent,
                                    'color': color,
                                }
                            )
                            pie_colors.append(color)
                            pie_sizes.append(percent)

                    # Store pie chart data for custom rendering
                    node_data['pie_data'] = pie_data
                    node_data['has_pie'] = True
                    node_data['pie_colors'] = pie_colors
                    node_data['pie_sizes'] = pie_sizes

                    # Generate SVG pie chart as Data URI
                    node_data['pie_svg'] = self.generate_pie_chart_svg(pie_data)

                    # Use transparent background to show pie chart
                    node_data['color'] = 'transparent'
                elif pop_counts:
                    # Exactly one population: solid colour for it
                    node_data['has_pie'] = False
                    (pop,) = pop_counts
                    if pop == 'Unassigned':
                        node_data['color'] = DEFAULT_MEDIAN_COLOR
                    else:
                        node_data['color'] = population_colors.get(
                            pop, DEFAULT_NODE_COLOR
                        )
                else:
                    node_data['has_pie'] = False
                    node_data['color'] = DEFAULT_NODE_COLOR
            else:
                node_data['has_pie'] = False
                if is_median:
                    node_data['color'] = median_vector_color
                else:
                    node_data['color'] = DEFAULT_NODE_COLOR  # lightblue

            # Everything the hover tooltip needs, rendered entirely in the
            # browser from this payload -- no server round trip per hover.
            node_data['tooltip'] = build_node_tooltip(
                label or node, hap, is_median, population_mapping
            )

            # Create element with position
            element = {
                'data': node_data,
                'position': {
                    'x': float(pos[0] * 100),
                    'y': float(pos[1] * 100),
                },  # Scale for visibility
                'grabbable': True,
            }

            elements.append(element)

        # Node centres and radii in canvas px, for fitting the tick comb.
        centres = {
            el['data']['id']: (el['position']['x'], el['position']['y'])
            for el in elements
        }
        radii = {el['data']['id']: el['data']['size'] / 2 for el in elements}

        # Create edges
        for u, v in graph.edges():
            # Get mutation count from distance attribute
            distance = graph[u][v].get('distance', 1)
            # Get weight for layout purposes (should be uniform)
            weight = graph[u][v].get('weight', 1.0)

            ticks = (
                format_edge_ticks(int(distance), max_tick_marks)
                if show_edge_labels
                else ''
            )
            # Initial fit from the server's idea of node size; the browser
            # refits from rendered widths (see callbacks/ui.py).
            tick_font = 0.0
            if ticks:
                (ux, uy), (vx, vy) = centres[u], centres[v]
                tick_font = fit_tick_font_size(
                    math.hypot(vx - ux, vy - uy), radii[u], radii[v], int(distance)
                )

            edge_data = {
                'id': f'{u}-{v}',
                'source': u,
                'target': v,
                'distance': distance,  # Mutation count for labels and proportional layout
                'weight': weight,  # Uniform weight for standard layouts
                'label': str(int(distance))
                if show_edge_labels and distance > 0
                else '',
                'ticks': ticks,
                'tick_font': tick_font,
            }

            elements.append({'data': edge_data})

        self.elements = elements
        return elements

    def create_stylesheet(
        self,
        population_colors: Optional[Dict[str, str]] = None,
        median_vector_color: str = DEFAULT_MEDIAN_COLOR,
    ) -> List[Dict]:
        """
        Create Cytoscape stylesheet for network visualization.

        Parameters
        ----------
        population_colors : Dict[str, str], optional
            Color mapping for populations.
        median_vector_color : str, default=DEFAULT_MEDIAN_COLOR
            Color for median vector nodes.

        Returns
        -------
        List[Dict]
            List of stylesheet dictionaries.
        """
        stylesheet = [
            # Default node style - circular markers
            {
                'selector': 'node',
                'style': {
                    'content': 'data(label)',
                    'text-valign': 'center',
                    'text-halign': 'center',
                    'background-color': 'data(color)',
                    'width': 'data(size)',
                    'height': 'data(size)',
                    'shape': 'ellipse',
                    'border-width': 3,
                    'border-color': POP_INK,
                    'font-family': LABEL_FONT,
                    'font-size': '10px',
                    'font-weight': 'bold',
                    'color': POP_INK,
                    'text-outline-width': 2,
                    'text-outline-color': POP_PAPER,
                },
            },
            # Median vector style - also circular but distinguished by color.
            # '[?field]' is the Cytoscape.js truthy test; '[field = true]'
            # is not valid selector syntax and silently matched every node,
            # painting the whole network median-grey.
            {
                'selector': 'node[?is_median]',
                'style': {
                    'shape': 'ellipse',
                    'background-color': median_vector_color,
                },
            },
            # Default edge style
            {
                'selector': 'edge',
                'style': {
                    'width': 'mapData(weight, 1, 10, 3, 1)',
                    'line-color': POP_NAVY,
                    'target-arrow-color': POP_NAVY,
                    'curve-style': 'bezier',
                    'opacity': 0.75,
                },
            },
            # Edge label style
            {
                'selector': 'edge[label]',
                'style': {
                    'label': 'data(label)',
                    'font-family': LABEL_FONT,
                    'font-size': '10px',
                    'text-background-color': POP_PAPER,
                    'text-background-opacity': 0.7,
                    'text-background-padding': '3px',
                    'color': POP_INK,
                },
            },
            # Highlighted/selected node style
            {
                'selector': 'node:selected',
                'style': {
                    'border-width': 6,
                    'border-color': POP_AMBER,
                    'z-index': 999,  # Bring to front
                },
            },
        ]

        self.stylesheet = stylesheet
        return stylesheet

    def create_pie_stylesheet(self, population_colors: Dict[str, str]) -> List[Dict]:
        """
        Create stylesheet with pie chart support for nodes.

        Parameters
        ----------
        population_colors : Dict[str, str]
            Color mapping for populations.

        Returns
        -------
        List[Dict]
            List of stylesheet rules for pie chart nodes.
        """
        pie_styles = []

        # Style for pie chart nodes - use SVG background image
        pie_styles.append(
            {
                'selector': 'node[pie_svg]',
                'style': {
                    'background-color': 'transparent',
                    'background-opacity': 0.1,
                    'background-image': 'data(pie_svg)',
                    'background-fit': 'contain',
                    'background-clip': 'node',
                    'border-width': 2,
                    'border-color': POP_INK,
                },
            }
        )

        # Override border for selected pie chart nodes
        # This must come AFTER the node[pie_svg] style to have higher specificity
        pie_styles.append(
            {
                'selector': 'node[pie_svg]:selected',
                'style': {
                    'border-width': 4,
                    'border-color': POP_AMBER,
                    'z-index': 999,
                },
            }
        )

        return pie_styles

    def generate_population_colors(self, populations: List[str]) -> Dict[str, str]:
        """
        Generate distinct colors for populations using HSV color space.

        Parameters
        ----------
        populations : List[str]
            List of population names.

        Returns
        -------
        Dict[str, str]
            Dictionary mapping population names to hex colors.
        """
        from .style import generate_population_colors as shared

        return shared(populations)


def create_cytoscape_network(
    network: HaplotypeNetwork,
    layout: Optional[Dict[str, Tuple[float, float]]] = None,
    population_colors: Optional[Dict[str, str]] = None,
    population_mapping: Optional[Dict[str, str]] = None,
    node_size_scale: float = 20.0,
    show_labels: bool = True,
    show_edge_labels: bool = True,
    median_vector_color: str = DEFAULT_MEDIAN_COLOR,
    node_labels: Optional[Dict[str, str]] = None,
    show_edge_ticks: bool = True,
    edge_tick_threshold: int = DEFAULT_TICK_THRESHOLD,
    max_tick_marks: int = MAX_TICK_MARKS,
) -> Tuple[List[Dict], List[Dict]]:
    """
    Create Cytoscape elements and stylesheet for a haplotype network.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object to visualize.
    layout : Dict[str, Tuple[float, float]], optional
        Pre-computed node positions {node_id: (x, y)}.
    population_colors : Dict[str, str], optional
        Color mapping for populations.
    population_mapping : Dict[str, str], optional
        Mapping of sample_id to population {sample_id: population}.
    node_size_scale : float, default=20.0
        Scaling factor for node sizes.
    show_labels : bool, default=True
        Whether to show node labels.
    show_edge_labels : bool, default=True
        Whether to show edge labels with mutation counts.
    median_vector_color : str, default=DEFAULT_MEDIAN_COLOR
        Color for median vector nodes.
    node_labels : Dict[str, str], optional
        Custom labels for nodes {node_id: label}.
    show_edge_ticks : bool, default=True
        Draw mutation counts as tick marks instead of numerals.
    edge_tick_threshold : int, default=DEFAULT_TICK_THRESHOLD
        Edges with more mutations than this fall back to a numeral.
    max_tick_marks : int, default=MAX_TICK_MARKS
        Cap on tick marks emitted per edge.

    Returns
    -------
    Tuple[List[Dict], List[Dict]]
        Tuple of (elements, stylesheet) for Cytoscape component.
    """
    plotter = InteractiveCytoscapePlotter(network)

    # Generate population colors if needed and populations exist
    if population_colors is None:
        # Check if any haplotypes have population data
        populations = set()
        for node in network._graph.nodes():
            if not network.is_median_vector(node):
                hap = network.get_haplotype(node)
                if hap:
                    pop_counts = hap.get_frequency_by_population()
                    if pop_counts:
                        populations.update(pop_counts.keys())

        if populations:
            population_colors = plotter.generate_population_colors(list(populations))

    elements = plotter.create_elements(
        layout=layout,
        node_size_scale=node_size_scale,
        population_colors=population_colors,
        population_mapping=population_mapping,
        show_labels=show_labels,
        show_edge_labels=show_edge_labels,
        median_vector_color=median_vector_color,
        node_labels=node_labels,
        max_tick_marks=max_tick_marks,
    )

    stylesheet = plotter.create_stylesheet(
        population_colors=population_colors,
        median_vector_color=median_vector_color,
    )

    # Add pie chart styles if we have population colors
    if population_colors:
        stylesheet.extend(plotter.create_pie_stylesheet(population_colors))

    # Must come last: these override the catch-all edge[label] rule.
    if show_edge_labels:
        stylesheet.extend(
            create_edge_tick_stylesheet(edge_tick_threshold, show_edge_ticks)
        )

    return elements, stylesheet
