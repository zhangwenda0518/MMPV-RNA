"""Unit tests for layout algorithms."""

import json
import os
import tempfile

import pytest

from pypopart.core.graph import HaplotypeNetwork
from pypopart.core.haplotype import Haplotype
from pypopart.core.sequence import Sequence
from pypopart.layout.algorithms import (
    CircularLayout,
    ForceDirectedLayout,
    HierarchicalLayout,
    KamadaKawaiLayout,
    LayoutAlgorithm,
    LayoutManager,
    ManualLayout,
    RadialLayout,
    SpectralLayout,
)


@pytest.fixture
def simple_network():
    """Create a simple haplotype network for testing."""
    network = HaplotypeNetwork('TestNetwork')

    # Create haplotypes
    seq1 = Sequence('hap1', 'ATCG')
    seq2 = Sequence('hap2', 'ATCC')
    seq3 = Sequence('hap3', 'GTCG')

    hap1 = Haplotype(seq1, sample_ids=['s1', 's2'])
    hap2 = Haplotype(seq2, sample_ids=['s3'])
    hap3 = Haplotype(seq3, sample_ids=['s4'])

    network.add_haplotype(hap1)
    network.add_haplotype(hap2)
    network.add_haplotype(hap3)

    # Create a simple star topology
    network.add_edge('hap1', 'hap2', distance=1)
    network.add_edge('hap1', 'hap3', distance=1)

    return network


@pytest.fixture
def disconnected_network():
    """Create a disconnected network for testing."""
    network = HaplotypeNetwork('Disconnected')

    seq1 = Sequence('hap1', 'ATCG')
    seq2 = Sequence('hap2', 'GTCC')
    seq3 = Sequence('hap3', 'CCGG')

    hap1 = Haplotype(seq1, sample_ids=['s1'])
    hap2 = Haplotype(seq2, sample_ids=['s2'])
    hap3 = Haplotype(seq3, sample_ids=['s3'])

    network.add_haplotype(hap1)
    network.add_haplotype(hap2)
    network.add_haplotype(hap3)

    # Only connect hap1 and hap2, leaving hap3 disconnected
    network.add_edge('hap1', 'hap2', distance=1)

    return network


class TestLayoutAlgorithm:
    """Test base LayoutAlgorithm class."""

    def test_init(self, simple_network):
        """Test initialization."""
        algo = LayoutAlgorithm(simple_network)
        assert algo.network == simple_network
        assert algo.graph == simple_network._graph

    def test_compute_not_implemented(self, simple_network):
        """Test that compute raises NotImplementedError."""
        algo = LayoutAlgorithm(simple_network)

        with pytest.raises(NotImplementedError):
            algo.compute()

    def test_save_load_layout(self, simple_network):
        """Test saving and loading layouts."""
        algo = LayoutAlgorithm(simple_network)

        layout = {'hap1': (0.0, 0.0), 'hap2': (1.0, 0.0), 'hap3': (0.5, 1.0)}

        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(tmpdir, 'layout.json')

            # Save layout
            algo.save_layout(layout, filename)
            assert os.path.exists(filename)

            # Load layout
            loaded = LayoutAlgorithm.load_layout(filename)
            assert loaded == layout

            # Verify JSON structure
            with open(filename, 'r') as f:
                data = json.load(f)
                assert 'hap1' in data
                assert data['hap1'] == [0.0, 0.0]


class TestForceDirectedLayout:
    """Test ForceDirectedLayout class."""

    def test_compute(self, simple_network):
        """Test force-directed layout computation."""
        algo = ForceDirectedLayout(simple_network)
        layout = algo.compute()

        assert len(layout) == 3
        assert 'hap1' in layout
        assert 'hap2' in layout
        assert 'hap3' in layout

        # Check that positions are tuples of floats
        for _node, pos in layout.items():
            assert isinstance(pos, tuple)
            assert len(pos) == 2
            assert all(isinstance(x, (int, float)) for x in pos)

    def test_compute_with_parameters(self, simple_network):
        """Test force-directed layout with parameters."""
        algo = ForceDirectedLayout(simple_network)

        layout1 = algo.compute(k=0.5, iterations=100, seed=42)
        layout2 = algo.compute(k=0.5, iterations=100, seed=42)

        # Same seed should give same layout
        for node in layout1:
            assert layout1[node] == layout2[node]

    def test_empty_network(self):
        """Test with empty network."""
        network = HaplotypeNetwork('Empty')
        algo = ForceDirectedLayout(network)

        layout = algo.compute()
        assert len(layout) == 0


class TestCircularLayout:
    """Test CircularLayout class."""

    def test_compute(self, simple_network):
        """Test circular layout computation."""
        algo = CircularLayout(simple_network)
        layout = algo.compute()

        assert len(layout) == 3

        # Check that nodes are arranged in a circle
        # All nodes should be approximately same distance from center
        distances = []
        for pos in layout.values():
            dist = (pos[0] ** 2 + pos[1] ** 2) ** 0.5
            distances.append(dist)

        # All distances should be similar (within tolerance)
        assert max(distances) - min(distances) < 0.01

    def test_compute_with_scale(self, simple_network):
        """Test circular layout with custom scale."""
        algo = CircularLayout(simple_network)

        layout1 = algo.compute(scale=1.0)
        layout2 = algo.compute(scale=2.0)

        # Larger scale should give larger distances
        dist1 = sum((pos[0] ** 2 + pos[1] ** 2) ** 0.5 for pos in layout1.values())
        dist2 = sum((pos[0] ** 2 + pos[1] ** 2) ** 0.5 for pos in layout2.values())

        assert dist2 > dist1

    def test_compute_with_center(self, simple_network):
        """Test circular layout with custom center."""
        algo = CircularLayout(simple_network)

        center = (5.0, 5.0)
        layout = algo.compute(center=center)

        # Average position should be near the center
        avg_x = sum(pos[0] for pos in layout.values()) / len(layout)
        avg_y = sum(pos[1] for pos in layout.values()) / len(layout)

        assert abs(avg_x - center[0]) < 0.1
        assert abs(avg_y - center[1]) < 0.1


class TestRadialLayout:
    """Test RadialLayout class."""

    def test_compute(self, simple_network):
        """Test radial layout computation."""
        algo = RadialLayout(simple_network)
        layout = algo.compute()

        assert len(layout) == 3

        # hap1 should be at center (most connected)
        assert layout['hap1'] == (0.0, 0.0)

        # Other nodes should be at distance 1
        for node in ['hap2', 'hap3']:
            x, y = layout[node]
            dist = (x**2 + y**2) ** 0.5
            assert abs(dist - 1.0) < 0.01

    def test_compute_with_center_node(self, simple_network):
        """Test radial layout with specified center."""
        algo = RadialLayout(simple_network)

        layout = algo.compute(center_node='hap2')

        # hap2 should be at center
        assert layout['hap2'] == (0.0, 0.0)

    def test_compute_with_scale(self, simple_network):
        """Test radial layout with custom scale."""
        algo = RadialLayout(simple_network)

        layout = algo.compute(scale=2.0)

        # Nodes at distance 1 from center should be at radius 2.0
        for node in ['hap2', 'hap3']:
            x, y = layout[node]
            dist = (x**2 + y**2) ** 0.5
            assert abs(dist - 2.0) < 0.01

    def test_disconnected_network(self, disconnected_network):
        """Test radial layout with disconnected network."""
        algo = RadialLayout(disconnected_network)
        layout = algo.compute()

        assert len(layout) == 3
        # Should handle disconnected nodes gracefully
        # All nodes should have positions
        for node in ['hap1', 'hap2', 'hap3']:
            assert node in layout


class TestHierarchicalLayout:
    """Test HierarchicalLayout class."""

    def test_compute_vertical(self, simple_network):
        """Test vertical hierarchical layout."""
        algo = HierarchicalLayout(simple_network)
        layout = algo.compute(vertical=True)

        assert len(layout) == 3

        # Root (hap1) should be at top
        assert layout['hap1'][1] > layout['hap2'][1]
        assert layout['hap1'][1] > layout['hap3'][1]

        # Leaf nodes should be at same level
        assert abs(layout['hap2'][1] - layout['hap3'][1]) < 0.01

    def test_compute_horizontal(self, simple_network):
        """Test horizontal hierarchical layout."""
        algo = HierarchicalLayout(simple_network)
        layout = algo.compute(vertical=False)

        assert len(layout) == 3

        # Root (hap1) should be at left
        assert layout['hap1'][0] < layout['hap2'][0]
        assert layout['hap1'][0] < layout['hap3'][0]

        # Leaf nodes should be at same level
        assert abs(layout['hap2'][0] - layout['hap3'][0]) < 0.01

    def test_compute_with_root(self, simple_network):
        """Test hierarchical layout with specified root."""
        algo = HierarchicalLayout(simple_network)

        layout = algo.compute(root_node='hap2', vertical=True)

        # hap2 should be at top
        assert layout['hap2'][1] >= layout['hap1'][1]

    def test_compute_with_dimensions(self, simple_network):
        """Test hierarchical layout with custom dimensions."""
        algo = HierarchicalLayout(simple_network)

        layout = algo.compute(width=4.0, height=3.0)

        # Check that layout fits within specified dimensions
        x_coords = [pos[0] for pos in layout.values()]
        y_coords = [pos[1] for pos in layout.values()]

        assert min(x_coords) >= 0
        assert max(x_coords) <= 4.0
        assert min(y_coords) >= 0
        assert max(y_coords) <= 3.0

    def test_disconnected_network(self, disconnected_network):
        """Test hierarchical layout with disconnected network."""
        algo = HierarchicalLayout(disconnected_network)
        layout = algo.compute()

        # Should compute positions for all nodes including disconnected ones
        assert len(layout) == 3
        for node in ['hap1', 'hap2', 'hap3']:
            assert node in layout


class TestKamadaKawaiLayout:
    """Test KamadaKawaiLayout class."""

    def test_compute(self, simple_network):
        """Test Kamada-Kawai layout computation."""
        algo = KamadaKawaiLayout(simple_network)
        layout = algo.compute()

        assert len(layout) == 3
        assert all(node in layout for node in ['hap1', 'hap2', 'hap3'])

    def test_compute_with_parameters(self, simple_network):
        """Test Kamada-Kawai layout with parameters."""
        algo = KamadaKawaiLayout(simple_network)

        layout = algo.compute(scale=2.0, center=(5.0, 5.0))

        assert len(layout) == 3


class TestSpectralLayout:
    """Test SpectralLayout class."""

    def test_compute(self, simple_network):
        """Test spectral layout computation."""
        algo = SpectralLayout(simple_network)
        layout = algo.compute()

        assert len(layout) == 3
        assert all(node in layout for node in ['hap1', 'hap2', 'hap3'])

        # Check that positions are tuples of floats
        for _node, pos in layout.items():
            assert isinstance(pos, tuple)
            assert len(pos) == 2
            assert all(isinstance(x, (int, float)) for x in pos)

    def test_compute_with_parameters(self, simple_network):
        """Test spectral layout with parameters."""
        algo = SpectralLayout(simple_network)

        layout = algo.compute(scale=2.0, center=(5.0, 5.0))

        assert len(layout) == 3

        # Check center is approximately correct
        avg_x = sum(pos[0] for pos in layout.values()) / len(layout)
        avg_y = sum(pos[1] for pos in layout.values()) / len(layout)

        assert abs(avg_x - 5.0) < 1.0
        assert abs(avg_y - 5.0) < 1.0


class TestManualLayout:
    """Test ManualLayout class."""

    def test_init_empty(self, simple_network):
        """Test initialization without positions."""
        algo = ManualLayout(simple_network)

        assert algo.positions == {}

    def test_init_with_positions(self, simple_network):
        """Test initialization with positions."""
        initial = {'hap1': (0.0, 0.0), 'hap2': (1.0, 0.0)}

        algo = ManualLayout(simple_network, initial_positions=initial)

        assert algo.positions == initial

    def test_compute_fills_missing(self, simple_network):
        """Test that compute fills in missing nodes."""
        initial = {'hap1': (0.0, 0.0)}
        algo = ManualLayout(simple_network, initial_positions=initial)

        layout = algo.compute()

        assert len(layout) == 3
        assert layout['hap1'] == (0.0, 0.0)
        assert 'hap2' in layout
        assert 'hap3' in layout

    def test_set_position(self, simple_network):
        """Test setting node position."""
        algo = ManualLayout(simple_network)

        algo.set_position('hap1', (1.5, 2.5))

        assert algo.positions['hap1'] == (1.5, 2.5)

    def test_set_position_invalid_node(self, simple_network):
        """Test setting position for non-existent node."""
        algo = ManualLayout(simple_network)

        with pytest.raises(ValueError, match='not in network'):
            algo.set_position('invalid', (0.0, 0.0))

    def test_move_node(self, simple_network):
        """Test moving node by offset."""
        algo = ManualLayout(simple_network)
        algo.set_position('hap1', (1.0, 1.0))

        algo.move_node('hap1', 0.5, -0.3)

        assert algo.positions['hap1'] == (1.5, 0.7)

    def test_move_node_without_position(self, simple_network):
        """Test moving node without initial position."""
        algo = ManualLayout(simple_network)

        with pytest.raises(ValueError, match='has no position set'):
            algo.move_node('hap1', 1.0, 1.0)


class TestLayoutManager:
    """Test LayoutManager class."""

    def test_init(self, simple_network):
        """Test initialization."""
        manager = LayoutManager(simple_network)

        assert manager.network == simple_network
        assert len(manager._algorithms) > 0

    def test_get_available_algorithms(self, simple_network):
        """Test getting available algorithms."""
        manager = LayoutManager(simple_network)

        algos = manager.get_available_algorithms()

        assert 'force_directed' in algos
        assert 'circular' in algos
        assert 'radial' in algos
        assert 'hierarchical' in algos
        assert 'kamada_kawai' in algos
        assert 'spectral' in algos
        assert 'manual' in algos

    def test_compute_layout_force_directed(self, simple_network):
        """Test computing force-directed layout."""
        manager = LayoutManager(simple_network)

        layout = manager.compute_layout('force_directed')

        assert len(layout) == 3
        assert all(node in layout for node in ['hap1', 'hap2', 'hap3'])

    def test_compute_layout_circular(self, simple_network):
        """Test computing circular layout."""
        manager = LayoutManager(simple_network)

        layout = manager.compute_layout('circular')

        assert len(layout) == 3

    def test_compute_layout_radial(self, simple_network):
        """Test computing radial layout."""
        manager = LayoutManager(simple_network)

        layout = manager.compute_layout('radial')

        assert len(layout) == 3

    def test_compute_layout_hierarchical(self, simple_network):
        """Test computing hierarchical layout."""
        manager = LayoutManager(simple_network)

        layout = manager.compute_layout('hierarchical')

        assert len(layout) == 3

    def test_compute_layout_kamada_kawai(self, simple_network):
        """Test computing Kamada-Kawai layout."""
        manager = LayoutManager(simple_network)

        layout = manager.compute_layout('kamada_kawai')

        assert len(layout) == 3

    def test_compute_layout_spectral(self, simple_network):
        """Test computing spectral layout."""
        manager = LayoutManager(simple_network)

        layout = manager.compute_layout('spectral')

        assert len(layout) == 3

    def test_compute_layout_with_params(self, simple_network):
        """Test computing layout with parameters."""
        manager = LayoutManager(simple_network)

        layout = manager.compute_layout('force_directed', iterations=100, seed=42)

        assert len(layout) == 3

    def test_compute_layout_invalid_algorithm(self, simple_network):
        """Test computing layout with invalid algorithm."""
        manager = LayoutManager(simple_network)

        with pytest.raises(ValueError, match='Unknown layout algorithm'):
            manager.compute_layout('invalid_algo')

    def test_save_load_layout(self, simple_network):
        """Test saving and loading layout."""
        manager = LayoutManager(simple_network)

        # Compute a layout
        layout = manager.compute_layout('circular')

        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(tmpdir, 'layout.json')

            # Save
            manager.save_layout(layout, filename)
            assert os.path.exists(filename)

            # Load
            loaded = manager.load_layout(filename)

            assert loaded == layout

    def test_caching_enabled(self, simple_network):
        """Test that caching works when enabled."""
        manager = LayoutManager(simple_network, enable_cache=True)

        # First computation
        layout1 = manager.compute_layout('spring', iterations=50, seed=42)

        # Second computation with same parameters should use cache
        layout2 = manager.compute_layout('spring', iterations=50, seed=42)

        # Results should be identical (same object from cache)
        assert layout1 == layout2

        # Different parameters should compute new layout
        layout3 = manager.compute_layout('spring', iterations=100, seed=42)
        assert layout3 != layout1  # Different iterations

    def test_caching_disabled(self, simple_network):
        """Test that caching can be disabled."""
        manager = LayoutManager(simple_network, enable_cache=False)

        # Compute layouts
        layout1 = manager.compute_layout('spring', iterations=50, seed=42)
        layout2 = manager.compute_layout('spring', iterations=50, seed=42)

        # Should be equal but potentially recomputed
        assert layout1 == layout2

    def test_clear_cache(self, simple_network):
        """Test clearing the layout cache."""
        manager = LayoutManager(simple_network, enable_cache=True)

        # Compute and cache
        layout1 = manager.compute_layout('circular')

        # Clear cache
        manager.clear_cache()

        # Should recompute (but result should be identical for deterministic algorithms)
        layout2 = manager.compute_layout('circular')
        assert layout1 == layout2


class TestEmptyNetwork:
    """Test layout algorithms with empty network."""

    def test_force_directed_empty(self):
        """Test force-directed layout with empty network."""
        network = HaplotypeNetwork('Empty')
        algo = ForceDirectedLayout(network)

        layout = algo.compute()
        assert len(layout) == 0

    def test_circular_empty(self):
        """Test circular layout with empty network."""
        network = HaplotypeNetwork('Empty')
        algo = CircularLayout(network)

        layout = algo.compute()
        assert len(layout) == 0

    def test_radial_empty(self):
        """Test radial layout with empty network."""
        network = HaplotypeNetwork('Empty')
        algo = RadialLayout(network)

        layout = algo.compute()
        assert len(layout) == 0

    def test_hierarchical_empty(self):
        """Test hierarchical layout with empty network."""
        network = HaplotypeNetwork('Empty')
        algo = HierarchicalLayout(network)

        layout = algo.compute()
        assert len(layout) == 0


class TestSingleNodeNetwork:
    """Test layout algorithms with single node."""

    @pytest.fixture
    def single_node_network(self):
        """Create network with single node."""
        network = HaplotypeNetwork('Single')
        seq = Sequence('hap1', 'ATCG')
        hap = Haplotype(seq, sample_ids=['s1'])
        network.add_haplotype(hap)
        return network

    def test_force_directed_single(self, single_node_network):
        """Test force-directed layout with single node."""
        algo = ForceDirectedLayout(single_node_network)
        layout = algo.compute()

        assert len(layout) == 1
        assert 'hap1' in layout

    def test_radial_single(self, single_node_network):
        """Test radial layout with single node."""
        algo = RadialLayout(single_node_network)
        layout = algo.compute()

        assert len(layout) == 1
        assert layout['hap1'] == (0.0, 0.0)

    def test_hierarchical_single(self, single_node_network):
        """Test hierarchical layout with single node."""
        algo = HierarchicalLayout(single_node_network)
        layout = algo.compute()

        assert len(layout) == 1


class TestSnapToGrid:
    """Quantising node positions onto a regular grid."""

    def test_rounds_to_nearest_intersection(self):
        """Coordinates move to the closest multiple, up or down."""
        from pypopart.layout.algorithms import snap_to_grid

        snapped = snap_to_grid({'a': (0.24, 0.26), 'b': (1.4, 1.6)}, 0.5)

        assert snapped['a'] == (0.0, 0.5)
        assert snapped['b'] == (1.5, 1.5)

    def test_handles_negative_coordinates(self):
        """Layouts are centred on the origin, so negatives are the norm."""
        from pypopart.layout.algorithms import snap_to_grid

        snapped = snap_to_grid({'a': (-0.24, -0.76)}, 0.5)

        assert snapped['a'] == (-0.0, -1.0)

    def test_is_idempotent(self):
        """It runs on every elements change, so re-running must be a no-op."""
        from pypopart.layout.algorithms import snap_to_grid

        once = snap_to_grid({'a': (0.31, -1.2), 'b': (2.7, 0.0)}, 0.25)
        twice = snap_to_grid(once, 0.25)

        assert once == twice

    @pytest.mark.parametrize('grid', [0, -1, 0.0])
    def test_zero_grid_passes_through(self, grid):
        """A grid of zero is how the feature is switched off."""
        from pypopart.layout.algorithms import snap_to_grid

        positions = {'a': (0.123, 4.567)}
        assert snap_to_grid(positions, grid) == positions

    def test_does_not_mutate_input(self):
        """Callers keep their own copy of the layout."""
        from pypopart.layout.algorithms import snap_to_grid

        positions = {'a': (0.31, 0.62)}
        snap_to_grid(positions, 0.5)

        assert positions == {'a': (0.31, 0.62)}


class TestPositionRoundTrip:
    """Cytoscape renders stored coordinates at 100x and reports them back."""

    @pytest.mark.parametrize('pixel', [137.0, -212.5, 0.0, 1e-3, 4321.25])
    def test_scaling_is_lossless_at_six_places(self, pixel):
        """Storing a dragged pixel position must not move the node."""
        from pypopart.gui.callbacks.layout import CYTOSCAPE_POSITION_SCALE

        stored = round(pixel / CYTOSCAPE_POSITION_SCALE, 6)
        assert stored * CYTOSCAPE_POSITION_SCALE == pytest.approx(pixel, abs=1e-4)


class TestMergeNodePositions:
    """Manual drags overlay the computed layout."""

    def test_drag_wins_over_layout(self):
        """The whole point of the separate store."""
        from pypopart.gui.serialization import merge_node_positions

        merged = merge_node_positions({'a': [1.0, 2.0]}, {'a': [9.0, 9.0]})

        assert merged['a'] == (9.0, 9.0)

    def test_layout_used_where_no_drag(self):
        """Untouched nodes keep their computed positions."""
        from pypopart.gui.serialization import merge_node_positions

        merged = merge_node_positions({'a': [1.0, 2.0], 'b': [3.0, 4.0]}, {'a': [9, 9]})

        assert merged['b'] == (3.0, 4.0)

    def test_stale_drags_are_dropped(self):
        """A drag for a node the new layout lacks must not resurrect it."""
        from pypopart.gui.serialization import merge_node_positions

        merged = merge_node_positions({'a': [1.0, 2.0]}, {'gone': [9.0, 9.0]})

        assert set(merged) == {'a'}

    @pytest.mark.parametrize('dragged', [None, {}])
    def test_no_drags_returns_the_layout(self, dragged):
        """An empty drag store is the normal state."""
        from pypopart.gui.serialization import merge_node_positions

        assert merge_node_positions({'a': [1.0, 2.0]}, dragged) == {'a': (1.0, 2.0)}

    def test_no_layout_returns_empty(self):
        """Nothing to render before a layout exists."""
        from pypopart.gui.serialization import merge_node_positions

        assert merge_node_positions(None, {'a': [1.0, 2.0]}) == {}


class TestResolveGridCollisions:
    """Snapping can drop two nodes on one intersection; spread them."""

    def test_coincident_nodes_are_separated(self):
        """The whole point: one node per intersection."""
        from pypopart.layout.algorithms import resolve_grid_collisions

        out = resolve_grid_collisions(
            {'a': (0.0, 0.0), 'b': (0.0, 0.0), 'c': (0.0, 0.0)}, 0.5
        )

        assert len(set(out.values())) == 3

    def test_result_stays_on_the_grid(self):
        """A displaced node lands on an intersection, not between them."""
        from pypopart.layout.algorithms import resolve_grid_collisions

        out = resolve_grid_collisions({'a': (0.0, 0.0), 'b': (0.0, 0.0)}, 0.25)

        for x, y in out.values():
            assert abs(x / 0.25 - round(x / 0.25)) < 1e-9
            assert abs(y / 0.25 - round(y / 0.25)) < 1e-9

    def test_moves_away_from_the_nearest_neighbour(self):
        """Displacing toward a neighbour would fold the edge back."""
        from pypopart.layout.algorithms import resolve_grid_collisions

        out = resolve_grid_collisions(
            {'pinned': (0.0, 0.0), 'moved': (0.0, 0.0), 'left': (-1.0, 0.0)},
            1.0,
            {'moved': ['left'], 'left': ['moved']},
            movable=['moved'],
        )

        assert out['moved'] == (1.0, 0.0)

    def test_takes_the_fewest_moves(self):
        """A free cell one step away beats one two steps away."""
        from pypopart.layout.algorithms import resolve_grid_collisions

        out = resolve_grid_collisions(
            {'pinned': (0.0, 0.0), 'moved': (0.0, 0.0)}, 1.0, movable=['moved']
        )
        x, y = out['moved']

        assert abs(x) + abs(y) == 1.0

    def test_steps_around_an_occupied_run(self):
        """Crossing occupied cells is penalised, so it goes around."""
        from pypopart.layout.algorithms import resolve_grid_collisions

        # A wall of pinned nodes to the right of the collision.
        positions = {f'wall{i}': (float(i), 0.0) for i in range(1, 4)}
        positions.update({'pinned': (0.0, 0.0), 'moved': (0.0, 0.0)})

        out = resolve_grid_collisions(
            positions, 1.0, {'moved': ['far']}, movable=['moved']
        )

        # Anywhere but further along the blocked row.
        assert out['moved'] not in {(1.0, 0.0), (2.0, 0.0), (3.0, 0.0)}
        assert len(set(out.values())) == len(positions)

    def test_pinned_nodes_do_not_move(self):
        """Dropping a node onto another displaces the one being dragged."""
        from pypopart.layout.algorithms import resolve_grid_collisions

        out = resolve_grid_collisions(
            {'kept': (2.0, 2.0), 'dragged': (2.0, 2.0)},
            0.5,
            {'dragged': ['kept'], 'kept': ['dragged']},
            movable=['dragged'],
        )

        assert out['kept'] == (2.0, 2.0)
        assert out['dragged'] != (2.0, 2.0)

    def test_no_collisions_is_a_no_op(self):
        """Already-distinct positions come back untouched."""
        from pypopart.layout.algorithms import resolve_grid_collisions

        positions = {'a': (0.0, 0.0), 'b': (1.0, 1.0), 'c': (-2.0, 0.5)}

        assert resolve_grid_collisions(positions, 0.5) == positions

    @pytest.mark.parametrize('grid', [0, -1])
    def test_disabled_grid_passes_through(self, grid):
        """With snapping off there is no lattice to de-conflict on."""
        from pypopart.layout.algorithms import resolve_grid_collisions

        positions = {'a': (0.0, 0.0), 'b': (0.0, 0.0)}

        assert resolve_grid_collisions(positions, grid) == positions

    def test_is_deterministic(self):
        """The same input must always give the same layout."""
        from pypopart.layout.algorithms import resolve_grid_collisions

        positions = dict.fromkeys('abcdefgh', (0.0, 0.0))
        first = resolve_grid_collisions(positions, 1.0)
        second = resolve_grid_collisions(positions, 1.0)

        assert first == second

    def test_handles_a_dense_pile_up(self):
        """Twenty nodes on one cell still each get their own."""
        from pypopart.layout.algorithms import resolve_grid_collisions

        positions = {f'n{i:02d}': (0.0, 0.0) for i in range(20)}
        out = resolve_grid_collisions(positions, 1.0)

        assert len(set(out.values())) == 20

    def test_does_not_mutate_input(self):
        """Callers keep their own copy."""
        from pypopart.layout.algorithms import resolve_grid_collisions

        positions = {'a': (0.0, 0.0), 'b': (0.0, 0.0)}
        resolve_grid_collisions(positions, 1.0)

        assert positions == {'a': (0.0, 0.0), 'b': (0.0, 0.0)}


class TestResolveEdgeOverlaps:
    """Two edges on one lattice vector hide each other; rotate one off."""

    @staticmethod
    def _overlapped(positions, adjacency):
        """Return the edges that are still overlapped after resolution."""
        from pypopart.layout.algorithms import _edge_is_overlapped

        cells = {n: (round(x), round(y)) for n, (x, y) in positions.items()}
        edges = {tuple(sorted((u, v))) for u in adjacency for v in adjacency[u]}
        return sorted(e for e in edges if _edge_is_overlapped(*e, cells, adjacency))

    def test_far_child_swings_off_the_shared_vector(self):
        """Parent at 0, children at 1 and 2: the outer child rotates."""
        from pypopart.layout.algorithms import resolve_edge_overlaps

        adjacency = {'p': ['a', 'b'], 'a': ['p'], 'b': ['p']}
        out = resolve_edge_overlaps(
            {'p': (0.0, 0.0), 'a': (1.0, 0.0), 'b': (2.0, 0.0)}, 1.0, adjacency
        )

        assert out['p'] == (0.0, 0.0)
        assert out['a'] == (1.0, 0.0)
        assert out['b'][1] != 0.0
        assert len(set(out.values())) == 3
        assert self._overlapped(out, adjacency) == []

    def test_result_stays_on_the_grid(self):
        """A rotated node lands on an intersection."""
        from pypopart.layout.algorithms import resolve_edge_overlaps

        adjacency = {'p': ['a', 'b'], 'a': ['p'], 'b': ['p']}
        out = resolve_edge_overlaps(
            {'p': (0.0, 0.0), 'a': (0.5, 0.0), 'b': (1.0, 0.0)}, 0.5, adjacency
        )

        for x, y in out.values():
            assert abs(x / 0.5 - round(x / 0.5)) < 1e-9
            assert abs(y / 0.5 - round(y / 0.5)) < 1e-9

    def test_edge_through_an_unrelated_node_is_cleared(self):
        """An edge must not pass under a node it does not connect."""
        from pypopart.layout.algorithms import resolve_edge_overlaps

        adjacency = {'p': ['q'], 'q': ['p'], 'x': []}
        out = resolve_edge_overlaps(
            {'p': (0.0, 0.0), 'q': (2.0, 0.0), 'x': (1.0, 0.0)}, 1.0, adjacency
        )

        assert out['x'] == (1.0, 0.0)
        assert self._overlapped(out, adjacency) == []

    def test_collinear_fan_gets_distinct_directions(self):
        """Five children in a row from one parent end up on five vectors."""
        from pypopart.layout.algorithms import resolve_edge_overlaps

        positions = {'p': (0.0, 0.0)}
        positions.update({f'c{i}': (float(i), 0.0) for i in range(1, 6)})
        adjacency = {'p': [f'c{i}' for i in range(1, 6)]}
        adjacency.update({f'c{i}': ['p'] for i in range(1, 6)})

        out = resolve_edge_overlaps(positions, 1.0, adjacency)

        assert len(set(out.values())) == 6
        assert self._overlapped(out, adjacency) == []

    def test_no_overlaps_is_a_no_op(self):
        """A clean layout comes back untouched."""
        from pypopart.layout.algorithms import resolve_edge_overlaps

        positions = {'p': (0.0, 0.0), 'a': (1.0, 0.0), 'b': (0.0, 1.0)}
        adjacency = {'p': ['a', 'b'], 'a': ['p'], 'b': ['p']}

        assert resolve_edge_overlaps(positions, 1.0, adjacency) == positions

    @pytest.mark.parametrize('grid', [0, -1])
    def test_disabled_grid_passes_through(self, grid):
        """With snapping off there is no lattice and nothing to do."""
        from pypopart.layout.algorithms import resolve_edge_overlaps

        positions = {'p': (0.0, 0.0), 'a': (1.0, 0.0), 'b': (2.0, 0.0)}
        adjacency = {'p': ['a', 'b'], 'a': ['p'], 'b': ['p']}

        assert resolve_edge_overlaps(positions, grid, adjacency) == positions

    def test_without_adjacency_passes_through(self):
        """No edges, no overlaps."""
        from pypopart.layout.algorithms import resolve_edge_overlaps

        positions = {'p': (0.0, 0.0), 'a': (1.0, 0.0), 'b': (2.0, 0.0)}

        assert resolve_edge_overlaps(positions, 1.0, None) == positions

    def test_pinned_nodes_do_not_move(self):
        """Only movable nodes rotate; a fully pinned overlap is left alone."""
        from pypopart.layout.algorithms import resolve_edge_overlaps

        positions = {'p': (0.0, 0.0), 'a': (1.0, 0.0), 'b': (2.0, 0.0)}
        adjacency = {'p': ['a', 'b'], 'a': ['p'], 'b': ['p']}

        out = resolve_edge_overlaps(positions, 1.0, adjacency, movable=['b'])
        assert out['p'] == (0.0, 0.0) and out['a'] == (1.0, 0.0)
        assert out['b'] != (2.0, 0.0)

        assert resolve_edge_overlaps(positions, 1.0, adjacency, movable=[]) == (
            positions
        )

    def test_is_deterministic(self):
        """The same input must always give the same layout."""
        from pypopart.layout.algorithms import resolve_edge_overlaps

        positions = {'p': (0.0, 0.0)}
        positions.update({f'c{i}': (float(i), 0.0) for i in range(1, 8)})
        adjacency = {'p': [f'c{i}' for i in range(1, 8)]}
        adjacency.update({f'c{i}': ['p'] for i in range(1, 8)})

        assert resolve_edge_overlaps(positions, 1.0, adjacency) == (
            resolve_edge_overlaps(positions, 1.0, adjacency)
        )

    def test_large_star_terminates_clean(self):
        """Thirty leaves packed in a block all get their own direction."""
        from pypopart.layout.algorithms import resolve_edge_overlaps

        positions = {'hub': (0.0, 0.0)}
        positions.update(
            {f'l{i:02d}': (float(i % 6 + 1), float(i // 6)) for i in range(30)}
        )
        adjacency = {'hub': [n for n in positions if n != 'hub']}
        adjacency.update({n: ['hub'] for n in positions if n != 'hub'})

        out = resolve_edge_overlaps(positions, 1.0, adjacency)

        assert len(set(out.values())) == 31
        assert self._overlapped(out, adjacency) == []

    def test_does_not_mutate_input(self):
        """Callers keep their own copy."""
        from pypopart.layout.algorithms import resolve_edge_overlaps

        positions = {'p': (0.0, 0.0), 'a': (1.0, 0.0), 'b': (2.0, 0.0)}
        adjacency = {'p': ['a', 'b'], 'a': ['p'], 'b': ['p']}
        resolve_edge_overlaps(positions, 1.0, adjacency)

        assert positions == {'p': (0.0, 0.0), 'a': (1.0, 0.0), 'b': (2.0, 0.0)}

    def test_shared_cell_after_failed_collision_search_does_not_crash(self):
        """A node moving off a cell it does not own must not raise."""
        from pypopart.layout.algorithms import resolve_edge_overlaps

        # 'b' and 'c' share (2, 0): only one of them holds the cell.
        positions = {'p': (0.0, 0.0), 'a': (1.0, 0.0), 'b': (2.0, 0.0), 'c': (2.0, 0.0)}
        adjacency = {'p': ['a', 'b', 'c'], 'a': ['p'], 'b': ['p'], 'c': ['p']}

        out = resolve_edge_overlaps(positions, 1.0, adjacency)

        assert self._overlapped(out, adjacency) == []
