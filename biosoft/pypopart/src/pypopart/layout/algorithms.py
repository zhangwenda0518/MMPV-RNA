"""
Layout algorithms for network visualization in PyPopART.

Provides various layout algorithms for positioning nodes in haplotype networks,
including force-directed, hierarchical, spectral, and custom layouts.

Algorithm selection guide:

For small networks (<50 nodes):
- KamadaKawaiLayout: Best quality, slow
- ForceDirectedLayout: Good quality, moderate speed

For medium networks (50-500 nodes):
- ForceDirectedLayout: Default choice, good balance
- SpectralLayout: Faster alternative, good quality
- HierarchicalLayout: Very fast, tree-like structure

For large networks (>500 nodes):
- SpectralLayout: Fast, maintains structure
- HierarchicalLayout: Fastest option
- CircularLayout: Simple, very fast

Special purposes:
- RadialLayout: Emphasize central node
- CircularLayout: Show connectivity patterns
"""

import json
import math
from typing import Dict, Iterable, List, Optional, Tuple

import networkx as nx
import numpy as np

from ..core.graph import HaplotypeNetwork


class LayoutAlgorithm:
    """
    Base class for layout algorithms.

    Provides interface for computing node positions in network visualizations.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.
    """

    def __init__(self, network: HaplotypeNetwork):
        """
        Initialize layout algorithm with a network.

        Parameters
        ----------
        network : HaplotypeNetwork
            HaplotypeNetwork object.
        """
        self.network = network
        self.graph = network._graph

    def compute(self, **kwargs) -> Dict[str, Tuple[float, float]]:
        """
        Compute node positions.

        Parameters
        ----------
        **kwargs : dict
            Algorithm-specific parameters.

        Returns
        -------
        Dict[str, Tuple[float, float]]
            Dictionary mapping node IDs to (x, y) positions.
        """
        raise NotImplementedError('Subclasses must implement compute()')

    def save_layout(
        self, layout: Dict[str, Tuple[float, float]], filename: str
    ) -> None:
        """
        Save layout to a JSON file.

        Parameters
        ----------
        layout : Dict[str, Tuple[float, float]]
            Node positions dictionary.
        filename : str
            Output filename.
        """
        # Convert tuples to lists for JSON serialization
        layout_serializable = {node: list(pos) for node, pos in layout.items()}

        with open(filename, 'w') as f:
            json.dump(layout_serializable, f, indent=2)

    @staticmethod
    def load_layout(filename: str) -> Dict[str, Tuple[float, float]]:
        """
        Load layout from a JSON file.

        Parameters
        ----------
        filename : str
            Input filename.

        Returns
        -------
        Dict[str, Tuple[float, float]]
            Dictionary mapping node IDs to (x, y) positions.
        """
        with open(filename, 'r') as f:
            layout_data = json.load(f)

        # Convert lists back to tuples
        return {node: tuple(pos) for node, pos in layout_data.items()}


class ForceDirectedLayout(LayoutAlgorithm):
    """
    Force-directed layout using spring algorithm.

    Simulates physical spring forces between connected nodes to create
    aesthetically pleasing layouts. Uses the Fruchterman-Reingold algorithm
    implemented in NetworkX's spring_layout.

    Notes
    -----
    **Performance**

    - Time complexity: O(iterations * N^2) where N is number of nodes
    - Typical runtime: ~25ms for 100 nodes, 50 iterations
    - Best for: Networks with 10-500 nodes
    - Quality: Good balance between speed and aesthetic quality

    For very large networks (>500 nodes), consider using:
    - HierarchicalLayout (fastest, ~0.1ms for 100 nodes)
    - CircularLayout (very fast, ~0.2ms for 100 nodes)
    - SpectralLayout (faster alternative to force-directed)
    """

    def compute(
        self,
        k: Optional[float] = None,
        iterations: int = 50,
        seed: Optional[int] = None,
        **kwargs,
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compute force-directed layout.

        Parameters
        ----------
        k : float, optional
            Optimal distance between nodes (None for auto).
            Smaller values bring nodes closer together.
        iterations : int, default=50
            Number of iterations for optimization.
            More iterations = better quality but slower.
            Default 50 is good for most networks.
        seed : int, optional
            Random seed for reproducibility.
        **kwargs : dict
            Additional parameters passed to spring_layout.

        Returns
        -------
        Dict[str, Tuple[float, float]]
            Node positions dictionary.
        """
        layout = nx.spring_layout(
            self.graph, k=k, iterations=iterations, seed=seed, **kwargs
        )
        # Convert numpy arrays to tuples
        return {node: tuple(pos) for node, pos in layout.items()}


class CircularLayout(LayoutAlgorithm):
    """
    Circular layout arranging nodes in a circle.

    Places nodes evenly spaced around a circle, useful for showing
    connectivity patterns.
    """

    def compute(
        self, scale: float = 1.0, center: Optional[Tuple[float, float]] = None, **kwargs
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compute circular layout.

        Parameters
        ----------
        scale : float, default=1.0
            Scale factor for the layout.
        center : Tuple[float, float], optional
            Center position (x, y).
        **kwargs : dict
            Additional parameters passed to circular_layout.

        Returns
        -------
        Dict[str, Tuple[float, float]]
            Node positions dictionary.
        """
        layout = nx.circular_layout(self.graph, scale=scale, center=center, **kwargs)
        # Convert numpy arrays to tuples
        return {node: tuple(pos) for node, pos in layout.items()}


class RadialLayout(LayoutAlgorithm):
    """
    Radial layout with center node and concentric rings.

    Places a central node at the origin and arranges other nodes
    in concentric circles based on distance from center.
    """

    def compute(
        self, center_node: Optional[str] = None, scale: float = 1.0, **kwargs
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compute radial layout.

        Parameters
        ----------
        center_node : str, optional
            Node to place at center (most connected if None).
        scale : float, default=1.0
            Scale factor for the layout.
        **kwargs : dict
            Additional parameters.

        Returns
        -------
        Dict[str, Tuple[float, float]]
            Node positions dictionary.
        """
        if not self.graph.nodes():
            return {}

        # Find center node if not specified
        if center_node is None:
            # Use node with highest degree
            center_node = max(self.graph.nodes(), key=lambda n: self.graph.degree(n))

        # Calculate distances from center using shortest path
        try:
            distances = nx.single_source_shortest_path_length(self.graph, center_node)
            # Add disconnected nodes at max distance + 1
            max_dist = max(distances.values()) if distances else 0
            for node in self.graph.nodes():
                if node not in distances:
                    distances[node] = max_dist + 1
        except nx.NetworkXError:
            # If graph is disconnected, use all nodes at max distance
            distances = dict.fromkeys(self.graph.nodes(), 1)
            distances[center_node] = 0

        # Group nodes by distance
        max_distance = max(distances.values()) if distances else 0
        rings: Dict[int, List[str]] = {i: [] for i in range(max_distance + 1)}

        for node, dist in distances.items():
            rings[dist].append(node)

        # Position nodes
        positions = {}
        positions[center_node] = (0.0, 0.0)

        for ring_idx, nodes in rings.items():
            if ring_idx == 0:
                continue

            radius = scale * ring_idx
            n_nodes = len(nodes)

            for i, node in enumerate(nodes):
                angle = 2 * np.pi * i / n_nodes
                x = radius * np.cos(angle)
                y = radius * np.sin(angle)
                positions[node] = (x, y)

        return positions


class HierarchicalLayout(LayoutAlgorithm):
    """
    Hierarchical layout arranging nodes in levels.

    Creates a tree-like structure with nodes arranged in horizontal
    levels based on distance from a root node.
    """

    def compute(
        self,
        root_node: Optional[str] = None,
        vertical: bool = True,
        width: float = 2.0,
        height: float = 2.0,
        **kwargs,
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compute hierarchical layout.

        Parameters
        ----------
        root_node : str, optional
            Root node for hierarchy (most connected if None).
        vertical : bool, default=True
            If True, levels are horizontal; if False, levels are vertical.
        width : float, default=2.0
            Total width of the layout.
        height : float, default=2.0
            Total height of the layout.
        **kwargs : dict
            Additional parameters.

        Returns
        -------
        Dict[str, Tuple[float, float]]
            Node positions dictionary.
        """
        if not self.graph.nodes():
            return {}

        # Find root node if not specified
        if root_node is None:
            root_node = max(self.graph.nodes(), key=lambda n: self.graph.degree(n))

        # Calculate distances from root using BFS
        try:
            distances = nx.single_source_shortest_path_length(self.graph, root_node)
            # Add disconnected nodes at max level + 1
            max_dist = max(distances.values()) if distances else 0
            for node in self.graph.nodes():
                if node not in distances:
                    distances[node] = max_dist + 1
        except nx.NetworkXError:
            # If graph is disconnected, use default distances
            distances = dict.fromkeys(self.graph.nodes(), 1)
            distances[root_node] = 0

        # Group nodes by level
        max_level = max(distances.values()) if distances else 0
        levels: Dict[int, List[str]] = {i: [] for i in range(max_level + 1)}

        for node, level in distances.items():
            levels[level].append(node)

        # Position nodes
        positions = {}

        for level_idx, nodes in levels.items():
            n_nodes = len(nodes)

            if vertical:
                # Horizontal levels (top to bottom)
                y = height * (1 - level_idx / max(max_level, 1))

                for i, node in enumerate(nodes):
                    if n_nodes == 1:
                        x = width / 2
                    else:
                        x = width * i / (n_nodes - 1)
                    positions[node] = (x, y)
            else:
                # Vertical levels (left to right)
                x = width * level_idx / max(max_level, 1)

                for i, node in enumerate(nodes):
                    if n_nodes == 1:
                        y = height / 2
                    else:
                        y = height * (1 - i / (n_nodes - 1))
                    positions[node] = (x, y)

        return positions


class KamadaKawaiLayout(LayoutAlgorithm):
    """
    Kamada-Kawai layout algorithm.

    Uses energy minimization to position nodes based on graph-theoretic
    distances. Produces high-quality layouts but is computationally expensive
    for large networks.

    Notes
    -----
    **Performance**

    - Time complexity: O(N^3) where N is number of nodes
    - Typical runtime: ~190ms for 100 nodes
    - Best for: Small networks (<50 nodes) where layout quality is critical
    - Quality: Excellent, minimizes stress based on graph distances

    For large networks, use ForceDirectedLayout or SpectralLayout instead.
    Kamada-Kawai can be very slow for networks with >100 nodes.
    """

    def compute(
        self, scale: float = 1.0, center: Optional[Tuple[float, float]] = None, **kwargs
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compute Kamada-Kawai layout.

        Parameters
        ----------
        scale : float, default=1.0
            Scale factor for the layout.
        center : Tuple[float, float], optional
            Center position (x, y).
        **kwargs : dict
            Additional parameters passed to kamada_kawai_layout.

        Returns
        -------
        Dict[str, Tuple[float, float]]
            Node positions dictionary.

        Warnings
        --------
        This algorithm can be very slow for large networks (>100 nodes).
        Consider using ForceDirectedLayout or SpectralLayout as faster alternatives.
        """
        layout = nx.kamada_kawai_layout(
            self.graph, scale=scale, center=center, **kwargs
        )
        # Convert numpy arrays to tuples
        return {node: tuple(pos) for node, pos in layout.items()}


class SpectralLayout(LayoutAlgorithm):
    """
    Spectral layout using graph Laplacian eigenvectors.

    Uses the eigenvectors of the graph Laplacian matrix to position nodes.
    This is a fast alternative to force-directed layouts that works well
    for large networks.

    Notes
    -----
    **Performance**

    - Time complexity: O(N^2) where N is number of nodes
    - Typical runtime: ~5-10ms for 100 nodes
    - Best for: Large networks (100-1000+ nodes)
    - Quality: Good, respects graph structure efficiently

    Spectral layout is much faster than Kamada-Kawai and comparable
    to force-directed layouts while maintaining good quality.
    Particularly effective for networks with clear clustering structure.
    """

    def compute(
        self,
        scale: float = 1.0,
        center: Optional[Tuple[float, float]] = None,
        dim: int = 2,
        **kwargs,
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compute spectral layout using graph Laplacian.

        Parameters
        ----------
        scale : float, default=1.0
            Scale factor for the layout.
        center : Tuple[float, float], optional
            Center position (x, y).
        dim : int, default=2
            Dimensionality of layout (default 2 for 2D visualization).
        **kwargs : dict
            Additional parameters passed to spectral_layout.

        Returns
        -------
        Dict[str, Tuple[float, float]]
            Node positions dictionary.
        """
        layout = nx.spectral_layout(
            self.graph, scale=scale, center=center, dim=dim, **kwargs
        )
        # Convert numpy arrays to tuples
        return {node: tuple(pos) for node, pos in layout.items()}


class ManualLayout(LayoutAlgorithm):
    """
    Manual layout with user-specified positions.

    Allows manual positioning of nodes or adjustment of existing layouts.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.
    initial_positions : Dict[str, Tuple[float, float]], optional
        Starting positions for nodes.
    """

    def __init__(
        self,
        network: HaplotypeNetwork,
        initial_positions: Optional[Dict[str, Tuple[float, float]]] = None,
    ):
        """
        Initialize manual layout.

        Parameters
        ----------
        network : HaplotypeNetwork
            HaplotypeNetwork object.
        initial_positions : Dict[str, Tuple[float, float]], optional
            Starting positions for nodes.
        """
        super().__init__(network)
        self.positions = initial_positions or {}

    def compute(self, **kwargs) -> Dict[str, Tuple[float, float]]:
        """
        Return current manual positions.

        Parameters
        ----------
        **kwargs : dict
            Ignored.

        Returns
        -------
        Dict[str, Tuple[float, float]]
            Node positions dictionary.
        """
        # Fill in missing nodes with default layout
        if len(self.positions) < len(self.graph.nodes()):
            default_layout = nx.spring_layout(self.graph)
            for node in self.graph.nodes():
                if node not in self.positions:
                    self.positions[node] = default_layout[node]

        return self.positions

    def set_position(self, node: str, position: Tuple[float, float]) -> None:
        """
        Set position for a specific node.

        Parameters
        ----------
        node : str
            Node ID.
        position : Tuple[float, float]
            (x, y) coordinates.
        """
        if node not in self.graph.nodes():
            raise ValueError(f"Node '{node}' not in network")

        self.positions[node] = position

    def move_node(self, node: str, dx: float, dy: float) -> None:
        """
        Move a node by a relative offset.

        Parameters
        ----------
        node : str
            Node ID.
        dx : float
            X offset.
        dy : float
            Y offset.
        """
        if node not in self.positions:
            raise ValueError(f"Node '{node}' has no position set")

        x, y = self.positions[node]
        self.positions[node] = (x + dx, y + dy)


class LayoutManager:
    """
    Manager for network layout computation and persistence.

    Provides high-level interface for computing, saving, and loading
    network layouts with optional caching.

    Parameters
    ----------
    network : HaplotypeNetwork
        The haplotype network to compute layouts for.
    enable_cache : bool, optional
        If True, cache layout results for repeated calls with same parameters.
        Default is True. Disable if network structure changes between calls.

    Examples
    --------
    >>> manager = LayoutManager(network)
    >>> layout = manager.compute_layout('spring', iterations=50, seed=42)
    >>> # Second call with same parameters uses cached result
    >>> layout2 = manager.compute_layout('spring', iterations=50, seed=42)
    """

    def __init__(self, network: HaplotypeNetwork, enable_cache: bool = True):
        """
        Initialize layout manager.

        Parameters
        ----------
        network : HaplotypeNetwork
            HaplotypeNetwork object.
        enable_cache : bool, default=True
            Enable caching of layout computations. Default True.
        """
        self.network = network
        self._enable_cache = enable_cache
        self._cache: Dict[str, Dict[str, Tuple[float, float]]] = {}
        self._algorithms = {
            'force_directed': ForceDirectedLayout,
            'spring': ForceDirectedLayout,  # Alias
            'circular': CircularLayout,
            'radial': RadialLayout,
            'hierarchical': HierarchicalLayout,
            'kamada_kawai': KamadaKawaiLayout,
            'spectral': SpectralLayout,
            'shell': CircularLayout,  # single-shell layout is circular
            'manual': ManualLayout,
        }

    def compute_layout(
        self, algorithm: str = 'force_directed', use_cache: bool = True, **kwargs
    ) -> Dict[str, Tuple[float, float]]:
        """
        Compute network layout using specified algorithm.

        Parameters
        ----------
        algorithm : str, default='force_directed'
            Layout algorithm name.
        use_cache : bool, default=True
            If True and caching is enabled, return cached result if available.
            Default True.
        **kwargs : dict
            Algorithm-specific parameters.

        Returns
        -------
        Dict[str, Tuple[float, float]]
            Node positions dictionary.

        Raises
        ------
        ValueError :
            If algorithm not recognized.

        Notes
        -----
        Results are cached based on algorithm name and parameters. To force
        recomputation, set use_cache=False or clear_cache().
        """
        if algorithm not in self._algorithms:
            available = ', '.join(self._algorithms.keys())
            raise ValueError(
                f"Unknown layout algorithm '{algorithm}'. Available: {available}"
            )

        # Check cache if enabled
        if self._enable_cache and use_cache:
            cache_key = self._make_cache_key(algorithm, kwargs)
            if cache_key in self._cache:
                return self._cache[cache_key]

        # Compute layout
        layout_class = self._algorithms[algorithm]
        layout_algo = layout_class(self.network)
        result = layout_algo.compute(**kwargs)

        # Store in cache if enabled
        if self._enable_cache:
            cache_key = self._make_cache_key(algorithm, kwargs)
            self._cache[cache_key] = result

        return result

    def _make_cache_key(self, algorithm: str, kwargs: Dict) -> str:
        """
        Create a cache key from algorithm name and parameters.

        Parameters
        ----------
        algorithm : str
            Algorithm name.
        kwargs : dict
            Algorithm parameters.

        Returns
        -------
        str
            Cache key string.
        """
        import json

        # Sort kwargs for consistent keys
        sorted_kwargs = json.dumps(kwargs, sort_keys=True, default=str)
        return f'{algorithm}:{sorted_kwargs}'

    def clear_cache(self) -> None:
        """
        Clear the layout cache.

        Use this after the network structure has changed to ensure
        fresh layout computations.
        """
        self._cache.clear()

    def save_layout(
        self, layout: Dict[str, Tuple[float, float]], filename: str
    ) -> None:
        """
        Save layout to file.

        Parameters
        ----------
        layout : Dict[str, Tuple[float, float]]
            Node positions dictionary.
        filename : str
            Output filename (JSON format).
        """
        algo = LayoutAlgorithm(self.network)
        algo.save_layout(layout, filename)

    def load_layout(self, filename: str) -> Dict[str, Tuple[float, float]]:
        """
        Load layout from file.

        Parameters
        ----------
        filename : str
            Input filename (JSON format).

        Returns
        -------
        Dict[str, Tuple[float, float]]
            Node positions dictionary.
        """
        return LayoutAlgorithm.load_layout(filename)

    def get_available_algorithms(self) -> List[str]:
        """
        Get list of available layout algorithms.

        Returns
        -------
        List[str]
            List of algorithm names.
        """
        return list(self._algorithms.keys())


def snap_to_grid(
    positions: Dict[str, Tuple[float, float]], grid_size: float
) -> Dict[str, Tuple[float, float]]:
    """
    Quantise node positions onto a regular grid.

    Applied both to computed layouts and to positions persisted after a
    manual drag, so a snapped network stays snapped either way.

    Parameters
    ----------
    positions : Dict[str, Tuple[float, float]]
        Node positions to quantise.
    grid_size : float
        Grid spacing in the same units as ``positions``. Values of zero or
        less return the positions unchanged, which is how the feature is
        switched off.

    Returns
    -------
    Dict[str, Tuple[float, float]]
        New mapping with every coordinate on the nearest grid intersection.
    """
    if grid_size <= 0:
        return dict(positions)

    return {
        node: (
            round(float(pos[0]) / grid_size) * grid_size,
            round(float(pos[1]) / grid_size) * grid_size,
        )
        for node, pos in positions.items()
    }


#: Extra cost charged for routing a displaced node through a cell that is
#: already taken, so a detour around a cluster wins over a path straight
#: across it when both take the same number of moves.
OCCUPIED_STEP_PENALTY = 3.0

#: How far, in grid cells, to search for a free intersection before giving
#: up and leaving a node where it landed. Bounds the search on a network
#: dense enough to have no free cell nearby.
MAX_COLLISION_SEARCH = 12

#: The four grid directions a node may step in.
_GRID_STEPS = ((1, 0), (-1, 0), (0, 1), (0, -1))


def resolve_grid_collisions(
    positions: Dict[str, Tuple[float, float]],
    grid_size: float,
    neighbours: Optional[Dict[str, List[str]]] = None,
    movable: Optional[Iterable[str]] = None,
) -> Dict[str, Tuple[float, float]]:
    """
    Spread nodes that snapped onto the same grid intersection.

    Snapping quantises positions, so nodes that were merely close end up
    exactly coincident and one hides the other. Each surplus node is
    walked outwards to the nearest free intersection, preferring to move
    away from its closest connected neighbour so an edge is not folded
    back over itself.

    Ties on distance are broken by that away-direction, and a route that
    crosses occupied cells is charged :data:`OCCUPIED_STEP_PENALTY` per
    crossing, so a node steps around a cluster rather than through it.

    Parameters
    ----------
    positions : Dict[str, Tuple[float, float]]
        Node positions, already snapped to the grid.
    grid_size : float
        Grid spacing, in the same units as ``positions``. Zero or less
        returns the positions unchanged.
    neighbours : Dict[str, List[str]], optional
        Adjacency, used to decide which way is "away". Nodes with no
        neighbours simply take the nearest free cell.
    movable : iterable of str, optional
        Nodes allowed to move. Everything else holds its cell, which is
        how a single dragged node is displaced rather than the network
        rearranging itself around it. Defaults to every node.

    Returns
    -------
    Dict[str, Tuple[float, float]]
        New positions with at most one node per intersection.
    """
    if grid_size <= 0 or not positions:
        return dict(positions)

    cells = {
        node: (round(pos[0] / grid_size), round(pos[1] / grid_size))
        for node, pos in positions.items()
    }
    movable_set = (
        set(cells) if movable is None else {node for node in movable if node in cells}
    )

    # Pinned nodes claim their cell first; the rest are placed in a
    # stable order so the same input always gives the same output.
    taken = {}
    for node in sorted(cells):
        if node not in movable_set:
            taken.setdefault(cells[node], node)

    resolved = {node: cells[node] for node in cells if node not in movable_set}

    for node in sorted(movable_set):
        start = cells[node]
        if start not in taken:
            taken[start] = node
            resolved[node] = start
            continue

        target = _nearest_free_cell(start, taken, _away_vector(node, cells, neighbours))
        taken[target] = node
        resolved[node] = target

    return {
        node: (grid_cell[0] * grid_size, grid_cell[1] * grid_size)
        for node, grid_cell in resolved.items()
    }


def _away_vector(
    node: str,
    cells: Dict[str, Tuple[int, int]],
    neighbours: Optional[Dict[str, List[str]]],
) -> Tuple[float, float]:
    """
    Point away from a node's closest connected neighbour.

    Parameters
    ----------
    node : str
        Node being displaced.
    cells : Dict[str, Tuple[int, int]]
        Grid cell of every node.
    neighbours : Dict[str, List[str]], optional
        Adjacency.

    Returns
    -------
    Tuple[float, float]
        Unit-ish direction to prefer, or ``(0.0, 0.0)`` when the node has
        no neighbours to move away from.
    """
    linked = [n for n in (neighbours or {}).get(node, []) if n in cells]
    if not linked:
        return (0.0, 0.0)

    x, y = cells[node]
    nearest = min(linked, key=lambda n: (cells[n][0] - x) ** 2 + (cells[n][1] - y) ** 2)
    dx, dy = x - cells[nearest][0], y - cells[nearest][1]
    length = math.hypot(dx, dy)
    if length == 0:
        return (0.0, 0.0)
    return (dx / length, dy / length)


def _nearest_free_cell(
    start: Tuple[int, int],
    taken: Dict[Tuple[int, int], str],
    away: Tuple[float, float],
) -> Tuple[int, int]:
    """
    Find the cheapest unoccupied cell reachable from a start cell.

    A uniform-cost search over the grid: one unit per step, plus a
    penalty for stepping through a cell that is already taken.

    Parameters
    ----------
    start : Tuple[int, int]
        Occupied cell the node landed on.
    taken : Dict[Tuple[int, int], str]
        Cells already claimed.
    away : Tuple[float, float]
        Preferred direction, used only to break ties.

    Returns
    -------
    Tuple[int, int]
        A free cell, or ``start`` if none was found within
        :data:`MAX_COLLISION_SEARCH` cells.
    """
    import heapq

    seen = {start}
    queue = [(0.0, 0.0, start, start)]
    while queue:
        cost, _, _, current = heapq.heappop(queue)
        if current not in taken and current != start:
            return current
        if cost >= MAX_COLLISION_SEARCH:
            continue

        for step_x, step_y in _GRID_STEPS:
            nxt = (current[0] + step_x, current[1] + step_y)
            if nxt in seen:
                continue
            seen.add(nxt)
            # Charge for crossing an occupied cell, not for landing free.
            step_cost = 1.0 + (OCCUPIED_STEP_PENALTY if nxt in taken else 0.0)

            # Ties on cost go to the cell best aligned with the away
            # direction, then to the cell itself so the result never
            # depends on iteration order.
            dx, dy = nxt[0] - start[0], nxt[1] - start[1]
            length = math.hypot(dx, dy) or 1.0
            alignment = (dx / length) * away[0] + (dy / length) * away[1]
            heapq.heappush(queue, (cost + step_cost, -alignment, nxt, nxt))

    return start


#: How many rounds of overlap fixing to run. Moving one node can expose a
#: new overlap elsewhere, so the pass repeats, but only this many times so
#: a pathological layout cannot loop.
MAX_OVERLAP_PASSES = 3

#: Furthest ring, in grid cells, searched around the anchor when rotating a
#: node off an overlapped edge. Bounds the candidate set on a big grid.
MAX_OVERLAP_RADIUS = 6

#: Smallest angle, in degrees, wanted between two edges at one node once
#: a node has been rotated. Two edges a few degrees apart hide each
#: other's tick marks nearly as well as exactly collinear ones, so a
#: placement that keeps this much clearance is tried first.
MIN_EDGE_SEPARATION_DEG = 15.0


def resolve_edge_overlaps(
    positions: Dict[str, Tuple[float, float]],
    grid_size: float,
    neighbours: Optional[Dict[str, List[str]]],
    movable: Optional[Iterable[str]] = None,
) -> Dict[str, Tuple[float, float]]:
    """
    Rotate nodes so that no two edges share a lattice vector.

    On a grid two edges leaving one node along the same direction lie on
    top of each other: a parent at the origin with children one and two
    cells to the right draws the second edge straight through the first
    child. The same happens when an unrelated node happens to sit on the
    segment between two connected nodes. Either way one edge is hidden.

    For each overlapped edge the lower-degree endpoint is swung around
    the other one to the nearest free cell whose direction differs from
    every other edge at that anchor, preferring the smallest change of
    angle and radius so the layout keeps its shape.

    Meant to run after :func:`resolve_grid_collisions` on a computed
    layout. Manual drags are left alone: a user who drops a node on a
    line did so on purpose.

    Parameters
    ----------
    positions : Dict[str, Tuple[float, float]]
        Node positions, already snapped and de-collided.
    grid_size : float
        Grid spacing, in the same units as ``positions``. Zero or less
        returns the positions unchanged.
    neighbours : Dict[str, List[str]], optional
        Adjacency. Without it there are no edges to overlap.
    movable : iterable of str, optional
        Nodes allowed to move. Defaults to every node.

    Returns
    -------
    Dict[str, Tuple[float, float]]
        New positions with overlapped edges rotated apart where a free
        direction could be found.
    """
    if grid_size <= 0 or not positions or not neighbours:
        return dict(positions)

    cells = {
        node: (round(pos[0] / grid_size), round(pos[1] / grid_size))
        for node, pos in positions.items()
    }
    movable_set = (
        set(cells) if movable is None else {node for node in movable if node in cells}
    )
    adjacency = {
        node: sorted(n for n in neighbours.get(node, []) if n in cells and n != node)
        for node in cells
    }
    edges = sorted({tuple(sorted((u, v))) for u in adjacency for v in adjacency[u]})
    if not edges:
        return dict(positions)

    taken = {}
    for node in sorted(cells):
        taken.setdefault(cells[node], node)

    for _ in range(MAX_OVERLAP_PASSES):
        changed = False
        # Longest edge first: of two collinear edges from one node, the
        # outer child is the one to swing, and once it has moved the
        # inner edge is clear and is skipped when its turn comes.
        for u, v in sorted(edges, key=lambda e: (-_lattice_length(cells, e), e)):
            if not _edge_is_overlapped(u, v, cells, adjacency):
                continue

            mover, anchor = _pick_mover(u, v, adjacency, movable_set)
            if mover is None:
                continue

            origin = cells[mover]
            offset = (origin[0] - cells[anchor][0], origin[1] - cells[anchor][1])
            candidates = _overlap_candidates(cells[anchor], offset, taken)
            # Well-separated placements first; if none exists settle for
            # any cell that at least ends the exact overlap.
            for separation in (MIN_EDGE_SEPARATION_DEG, 0.0):
                target = next(
                    (
                        cell
                        for cell in candidates
                        if _placement_is_clear(
                            mover, {**cells, mover: cell}, adjacency, edges, separation
                        )
                    ),
                    None,
                )
                if target is not None:
                    break
            if target is None:
                continue

            # A cell can be shared when collision resolution gave up on a
            # dense pile-up, so only release it if this node holds it.
            if taken.get(origin) == mover:
                del taken[origin]
            taken[target] = mover
            cells = {**cells, mover: target}
            changed = True

        if not changed:
            break

    return {
        node: (grid_cell[0] * grid_size, grid_cell[1] * grid_size)
        for node, grid_cell in cells.items()
    }


def _lattice_length(cells: Dict[str, Tuple[int, int]], edge: Tuple[str, str]) -> float:
    """
    Length of an edge in grid cells.

    Parameters
    ----------
    cells : Dict[str, Tuple[int, int]]
        Grid cell of every node.
    edge : Tuple[str, str]
        The edge's endpoints.

    Returns
    -------
    float
        Euclidean distance between the two cells.
    """
    (ax, ay), (bx, by) = cells[edge[0]], cells[edge[1]]
    return math.hypot(bx - ax, by - ay)


def _reduced_direction(start: Tuple[int, int], end: Tuple[int, int]) -> Tuple[int, int]:
    """
    Reduce a lattice vector to its primitive direction.

    Parameters
    ----------
    start : Tuple[int, int]
        Cell the vector leaves from.
    end : Tuple[int, int]
        Cell the vector points at.

    Returns
    -------
    Tuple[int, int]
        The vector divided by the gcd of its components, so any two
        collinear same-sense vectors compare equal. ``(0, 0)`` when the
        cells coincide.
    """
    dx, dy = end[0] - start[0], end[1] - start[1]
    divisor = math.gcd(abs(dx), abs(dy))
    if divisor == 0:
        return (0, 0)
    return (dx // divisor, dy // divisor)


def _cell_on_segment(
    cell: Tuple[int, int], start: Tuple[int, int], end: Tuple[int, int]
) -> bool:
    """
    Test whether a cell lies strictly inside a lattice segment.

    Parameters
    ----------
    cell : Tuple[int, int]
        Cell to test.
    start : Tuple[int, int]
        One end of the segment.
    end : Tuple[int, int]
        The other end.

    Returns
    -------
    bool
        True when the cell is on the segment and is neither endpoint.
    """
    if cell == start or cell == end:
        return False
    sx, sy = end[0] - start[0], end[1] - start[1]
    cx, cy = cell[0] - start[0], cell[1] - start[1]
    if sx * cy - sy * cx != 0:
        return False
    dot = sx * cx + sy * cy
    return 0 < dot < sx * sx + sy * sy


def _edge_is_overlapped(
    u: str,
    v: str,
    cells: Dict[str, Tuple[int, int]],
    adjacency: Dict[str, List[str]],
) -> bool:
    """
    Decide whether an edge is hidden under another edge or a node.

    Parameters
    ----------
    u : str
        One endpoint.
    v : str
        The other endpoint.
    cells : Dict[str, Tuple[int, int]]
        Grid cell of every node.
    adjacency : Dict[str, List[str]]
        Adjacency restricted to placed nodes.

    Returns
    -------
    bool
        True when another edge leaves either endpoint along the same
        direction, or when a third node sits on the segment.
    """
    for here, there in ((u, v), (v, u)):
        direction = _reduced_direction(cells[here], cells[there])
        for other in adjacency[here]:
            if other != there and (
                _reduced_direction(cells[here], cells[other]) == direction
            ):
                return True

    start, end = cells[u], cells[v]
    return any(
        _cell_on_segment(cell, start, end)
        for node, cell in cells.items()
        if node != u and node != v
    )


def _pick_mover(
    u: str,
    v: str,
    adjacency: Dict[str, List[str]],
    movable: set,
) -> Tuple[Optional[str], Optional[str]]:
    """
    Choose which endpoint of an overlapped edge to swing.

    Parameters
    ----------
    u : str
        One endpoint.
    v : str
        The other endpoint.
    adjacency : Dict[str, List[str]]
        Adjacency restricted to placed nodes.
    movable : set
        Nodes allowed to move.

    Returns
    -------
    Tuple[Optional[str], Optional[str]]
        ``(mover, anchor)``: the lower-degree movable endpoint and the
        one it rotates around, or ``(None, None)`` when neither may
        move.
    """
    candidates = sorted(
        (node for node in (u, v) if node in movable),
        key=lambda node: (len(adjacency[node]), node),
    )
    if not candidates:
        return None, None
    mover = candidates[0]
    return mover, v if mover == u else u


def _overlap_candidates(
    anchor: Tuple[int, int],
    offset: Tuple[int, int],
    taken: Dict[Tuple[int, int], str],
) -> List[Tuple[int, int]]:
    """
    List free cells around an anchor, closest in angle first.

    Parameters
    ----------
    anchor : Tuple[int, int]
        Cell the node rotates around.
    offset : Tuple[int, int]
        The node's current offset from the anchor.
    taken : Dict[Tuple[int, int], str]
        Cells already claimed.

    Returns
    -------
    List[Tuple[int, int]]
        Unoccupied cells within :data:`MAX_OVERLAP_RADIUS`, smallest
        move first. A move is scored as the turn from ``offset`` in
        eighth-turns plus the relative change in radius, so a small
        swing at the same distance beats a long slide outwards; ties go
        to the cell itself so the order is deterministic.
    """
    radius = math.hypot(*offset)
    reach = min(int(math.ceil(radius)) + 1, MAX_OVERLAP_RADIUS)

    def rank(cell: Tuple[int, int]) -> Tuple[float, Tuple[int, int]]:
        """
        Score one candidate cell; lower is a smaller move.

        Parameters
        ----------
        cell : Tuple[int, int]
            Candidate cell.

        Returns
        -------
        Tuple[float, Tuple[int, int]]
            The move score and the cell itself as a tie-break.
        """
        dx, dy = cell[0] - anchor[0], cell[1] - anchor[1]
        length = math.hypot(dx, dy)
        if radius == 0 or length == 0:
            angle = 0.0
        else:
            cosine = (dx * offset[0] + dy * offset[1]) / (length * radius)
            angle = math.acos(max(-1.0, min(1.0, cosine)))
        score = angle / (math.pi / 4) + abs(length - radius) / max(radius, 1.0)
        return (round(score, 9), cell)

    ring = [
        (anchor[0] + dx, anchor[1] + dy)
        for dx in range(-reach, reach + 1)
        for dy in range(-reach, reach + 1)
        if (dx, dy) != (0, 0)
    ]
    return sorted((cell for cell in ring if cell not in taken), key=rank)


def _placement_is_clear(
    node: str,
    cells: Dict[str, Tuple[int, int]],
    adjacency: Dict[str, List[str]],
    edges: List[Tuple[str, str]],
    min_separation_deg: float = 0.0,
) -> bool:
    """
    Check that a trial position for a node creates no overlap.

    Parameters
    ----------
    node : str
        Node that was moved.
    cells : Dict[str, Tuple[int, int]]
        Grid cells with the trial position applied.
    adjacency : Dict[str, List[str]]
        Adjacency restricted to placed nodes.
    edges : List[Tuple[str, str]]
        Every edge, so edges not touching ``node`` can be checked for
        now passing under it.
    min_separation_deg : float, default=0.0
        Additionally require this much angle between each of the node's
        edges and every other edge sharing an endpoint with it.

    Returns
    -------
    bool
        True when none of the node's edges overlap, no other edge runs
        through its new cell, and the angular clearance is met.
    """
    if any(
        _edge_is_overlapped(node, other, cells, adjacency) for other in adjacency[node]
    ):
        return False
    here = cells[node]
    if any(
        _cell_on_segment(here, cells[a], cells[b])
        for a, b in edges
        if node != a and node != b
    ):
        return False
    if min_separation_deg <= 0:
        return True

    threshold = math.radians(min_separation_deg)
    for neighbour in adjacency[node]:
        # Clearance at the node's own end and at the neighbour's end.
        for pivot, far in ((node, neighbour), (neighbour, node)):
            for other in adjacency[pivot]:
                if other != far and (
                    _angle_between(cells[pivot], cells[far], cells[other]) < threshold
                ):
                    return False
    return True


def _angle_between(
    pivot: Tuple[int, int], a: Tuple[int, int], b: Tuple[int, int]
) -> float:
    """
    Angle at a pivot cell between the rays to two other cells.

    Parameters
    ----------
    pivot : Tuple[int, int]
        Cell the two rays leave from.
    a : Tuple[int, int]
        End of the first ray.
    b : Tuple[int, int]
        End of the second ray.

    Returns
    -------
    float
        Angle in radians, in ``[0, pi]``; zero when either ray has no
        length.
    """
    ax, ay = a[0] - pivot[0], a[1] - pivot[1]
    bx, by = b[0] - pivot[0], b[1] - pivot[1]
    norm = math.hypot(ax, ay) * math.hypot(bx, by)
    if norm == 0:
        return 0.0
    return math.acos(max(-1.0, min(1.0, (ax * bx + ay * by) / norm)))
