"""
Haplotype network representation and analysis for PyPopART.

Provides the HaplotypeNetwork class for building and analyzing
haplotype networks from DNA sequence data.
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Set, Tuple

import networkx as nx

from .haplotype import Haplotype
from .sequence import Sequence


@dataclass
class NetworkStats:
    """Statistics about a haplotype network."""

    num_nodes: int
    num_edges: int
    num_haplotypes: int
    num_median_vectors: int
    total_samples: int
    diameter: int
    avg_degree: float
    num_components: int


class HaplotypeNetwork:
    """
    Represents a haplotype network using NetworkX.

    A haplotype network is a graph where nodes represent unique haplotypes
    (or inferred median vectors) and edges represent mutational relationships.
    Node sizes typically reflect haplotype frequencies, and edge weights
    represent genetic distances.

    Parameters
    ----------
    name : str, optional
        Optional name for the network.
    """

    def __init__(self, name: Optional[str] = None):
        """
        Initialize an empty haplotype network.

        Parameters
        ----------
        name : str, optional
            Optional name for the network.
        """
        self.name = name or 'HaplotypeNetwork'
        self._graph = nx.Graph()
        self._haplotype_map: Dict[str, Haplotype] = {}
        self._median_vectors: Set[str] = set()
        self.metadata: Dict[str, Any] = {}
        # Mutation counter for invalidating cached derived metrics
        self._version = 0
        self._metric_cache: Dict[str, Any] = {}

    @classmethod
    def from_serialized(cls, network_data: Dict) -> 'HaplotypeNetwork':
        """
        Reconstruct a HaplotypeNetwork from serialized data.

        Parameters
        ----------
        network_data : Dict
            Dictionary containing serialized network data with 'nodes' and 'edges'.

        Returns
        -------
        HaplotypeNetwork
            Reconstructed HaplotypeNetwork object.
        """
        network = cls()

        # Reconstruct haplotypes from nodes
        for node in network_data.get('nodes', []):
            node_id = node['id']
            sequence_data = node.get('sequence', '')
            is_median = node.get('is_median', False)
            sample_ids = node.get('sample_ids', node.get('samples', []))

            # Create a Sequence object
            seq = Sequence(id=node_id, data=sequence_data)

            # Create a Haplotype object
            haplotype = Haplotype(sequence=seq, sample_ids=sample_ids)

            # Add to network
            network.add_haplotype(haplotype, median_vector=is_median)

        # Add edges; older serializations carry only 'weight', which holds
        # the genetic distance, so fall back to it when 'distance' is absent.
        for edge in network_data.get('edges', []):
            source = edge['source']
            target = edge['target']
            weight = edge.get('weight', 1.0)
            distance = edge.get('distance', weight if 'weight' in edge else 0)
            network.add_edge(source, target, distance=distance, weight=weight)

        return network

    @classmethod
    def from_networkx(cls, graph: nx.Graph) -> 'HaplotypeNetwork':
        """
        Build a HaplotypeNetwork from a plain NetworkX graph.

        Intended for graphs loaded from exported network files (GraphML,
        GML, JSON), where node attributes such as ``sequence``, ``frequency``,
        ``median_vector`` and ``sample_ids`` may be present, possibly as
        strings after serialization round-trips.

        Parameters
        ----------
        graph : nx.Graph
            NetworkX graph with optional haplotype node/edge attributes.

        Returns
        -------
        HaplotypeNetwork
            Reconstructed HaplotypeNetwork object.
        """

        def _as_bool(value) -> bool:
            """
            Coerce a possibly-stringified flag to a bool.

            Parameters
            ----------
            value : Any
                Attribute value from the graph.

            Returns
            -------
            bool
                True when the value represents a true flag.
            """
            if isinstance(value, str):
                return value.strip().lower() in ('true', '1', 'yes')
            return bool(value)

        def _as_list(value) -> list:
            """
            Coerce a possibly-serialised sequence to a list.

            Parameters
            ----------
            value : Any
                Attribute value from the graph.

            Returns
            -------
            list
                The value as a list, empty when absent.
            """
            if value is None:
                return []
            if isinstance(value, str):
                return [item for item in value.split(',') if item]
            return list(value)

        nodes = [
            {
                'id': str(node_id),
                'sequence': attrs.get('sequence', ''),
                'is_median': _as_bool(
                    attrs.get('median_vector', attrs.get('is_median', False))
                ),
                'sample_ids': _as_list(attrs.get('sample_ids')),
            }
            for node_id, attrs in graph.nodes(data=True)
        ]
        edges = [
            {
                'source': str(source),
                'target': str(target),
                **{
                    key: attrs[key]
                    for key in ('weight', 'distance')
                    if key in attrs and attrs[key] is not None
                },
            }
            for source, target, attrs in graph.edges(data=True)
        ]

        network = cls.from_serialized({'nodes': nodes, 'edges': edges})
        network.metadata.update(graph.graph)
        return network

    @property
    def graph(self) -> nx.Graph:
        """
        Get the underlying NetworkX graph.

        Returns
        -------
        nx.Graph
            The underlying NetworkX graph.
        """
        return self._graph

    def add_haplotype(self, haplotype: Haplotype, median_vector: bool = False) -> None:
        """
        Add a haplotype as a node to the network.

        Parameters
        ----------
        haplotype : Haplotype
            Haplotype object to add.
        median_vector : bool, default=False
            Whether this is an inferred median vector.
        """
        node_id = haplotype.id

        if node_id in self._graph:
            raise ValueError(f"Node '{node_id}' already exists in network")

        # Store haplotype reference
        self._haplotype_map[node_id] = haplotype

        # Track median vectors
        if median_vector:
            self._median_vectors.add(node_id)

        self._touch()

        # Add node with attributes
        self._graph.add_node(
            node_id,
            haplotype=haplotype,
            frequency=haplotype.frequency,
            sequence=haplotype.data,
            median_vector=median_vector,
            sample_ids=haplotype.sample_ids,
            populations=list(haplotype.get_populations()),
        )

    def remove_haplotype(self, haplotype_id: str) -> None:
        """
        Remove a haplotype node from the network.

        Parameters
        ----------
        haplotype_id : str
            ID of haplotype to remove.

        Raises
        ------
        KeyError
            If haplotype not found.
        """
        if haplotype_id not in self._graph:
            raise KeyError(f"Haplotype '{haplotype_id}' not found in network")

        self._touch()
        self._graph.remove_node(haplotype_id)
        self._haplotype_map.pop(haplotype_id, None)
        self._median_vectors.discard(haplotype_id)

    def add_edge(
        self,
        source: str,
        target: str,
        distance: int = 0,
        weight: float = 1.0,
        **attributes,
    ) -> None:
        """
        Add an edge between two haplotypes.

        Parameters
        ----------
        source : str
            Source haplotype ID.
        target : str
            Target haplotype ID.
        distance : int, default=0
            Genetic distance between the two haplotypes.
        weight : float, default=1.0
            Edge weight used by layout and drawing code.
        **attributes : dict
            Additional edge attributes.
        """
        if source not in self._graph:
            raise KeyError(f"Source node '{source}' not found in network")
        if target not in self._graph:
            raise KeyError(f"Target node '{target}' not found in network")

        self._touch()
        self._graph.add_edge(
            source, target, weight=weight, distance=distance, **attributes
        )

    def remove_edge(self, source: str, target: str) -> None:
        """
        Remove an edge from the network.

        Parameters
        ----------
        source : str
            Source haplotype ID.
        target : str
            Target haplotype ID.

        Raises
        ------
        KeyError
            If edge not found.
        """
        if not self._graph.has_edge(source, target):
            raise KeyError(f'Edge ({source}, {target}) not found in network')

        self._touch()
        self._graph.remove_edge(source, target)

    def get_haplotype(self, haplotype_id: str) -> Haplotype:
        """
        Get haplotype by ID.

        Parameters
        ----------
        haplotype_id : str
            Haplotype identifier.

        Returns
        -------
        Haplotype
            Haplotype object.
            Raises :
            KeyError :
            If haplotype not found.
        """
        if haplotype_id not in self._haplotype_map:
            raise KeyError(f"Haplotype '{haplotype_id}' not found in network")

        return self._haplotype_map[haplotype_id]

    def has_node(self, haplotype_id: str) -> bool:
        """
        Check if node exists in network.

        Parameters
        ----------
        haplotype_id : str
            Haplotype identifier.

        Returns
        -------
        bool
            True if node exists.
        """
        return haplotype_id in self._graph

    def has_edge(self, source: str, target: str) -> bool:
        """
        Check if edge exists in network.

        Parameters
        ----------
        source : str
            Source haplotype ID.
        target : str
            Target haplotype ID.

        Returns
        -------
        bool
            True if edge exists.
        """
        return self._graph.has_edge(source, target)

    def get_edge_distance(self, source: str, target: str) -> float:
        """
        Get distance for an edge.

        Parameters
        ----------
        source : str
            Source haplotype ID.
        target : str
            Target haplotype ID.

        Returns
        -------
        float
            Edge distance.
            Raises :
            KeyError :
            If edge not found.
        """
        if not self.has_edge(source, target):
            raise KeyError(f'Edge ({source}, {target}) not found')

        return self._graph[source][target]['distance']

    def get_neighbors(self, haplotype_id: str) -> List[str]:
        """
        Get neighboring haplotype IDs.

        Parameters
        ----------
        haplotype_id : str
            Haplotype identifier.

        Returns
        -------
        List[str]
            List of neighbor IDs.
            Raises :
            KeyError :
            If haplotype not found.
        """
        if not self.has_node(haplotype_id):
            raise KeyError(f"Haplotype '{haplotype_id}' not found")

        return list(self._graph.neighbors(haplotype_id))

    def get_degree(self, haplotype_id: str) -> int:
        """
        Get degree (number of connections) for a node.

        Parameters
        ----------
        haplotype_id : str
            Haplotype identifier.

        Returns
        -------
        int
            Node degree.
        """
        if not self.has_node(haplotype_id):
            raise KeyError(f"Haplotype '{haplotype_id}' not found")

        return self._graph.degree[haplotype_id]

    @property
    def num_nodes(self) -> int:
        """
        Get number of nodes in network.

        Returns
        -------
        int
            Number of nodes in the network.
        """
        return self._graph.number_of_nodes()

    @property
    def num_edges(self) -> int:
        """
        Get number of edges in network.

        Returns
        -------
        int
            Number of edges in the network.
        """
        return self._graph.number_of_edges()

    @property
    def nodes(self) -> List[str]:
        """
        Get list of node IDs.

        Returns
        -------
        List[str]
            List of node IDs.
        """
        return list(self._graph.nodes())

    @property
    def edges(self) -> List[Tuple[str, str]]:
        """
        Get list of edges as (source, target) tuples.

        Returns
        -------
        List[Tuple[str, str]]
            List of edges as (source, target) tuples.
        """
        return list(self._graph.edges())

    @property
    def haplotypes(self) -> List[Haplotype]:
        """
        Get list of all haplotypes (excluding median vectors).

        Returns
        -------
        List[Haplotype]
            List of all haplotypes (excluding median vectors).
        """
        return [
            hap
            for hap_id, hap in self._haplotype_map.items()
            if hap_id not in self._median_vectors
        ]

    @property
    def median_vector_ids(self) -> List[str]:
        """
        Get list of median vector node IDs.

        Returns
        -------
        List[str]
            List of median vector node IDs.
        """
        return sorted(self._median_vectors)

    def is_median_vector(self, node_id: str) -> bool:
        """
        Check if a node is a median vector.

        Parameters
        ----------
        node_id : str
            Node identifier.

        Returns
        -------
        bool
            True if node is a median vector.
        """
        return node_id in self._median_vectors

    def is_connected(self) -> bool:
        """
        Check if network is fully connected.

        Returns
        -------
        bool
            True if all nodes are in one connected component.
        """
        return nx.is_connected(self._graph)

    def get_connected_components(self) -> List[Set[str]]:
        """
        Get connected components of the network.

        Returns
        -------
        List[Set[str]]
            List of sets, each containing node IDs in a component.
        """
        return [set(component) for component in nx.connected_components(self._graph)]

    def _touch(self) -> None:
        """Invalidate cached metrics after a structural change."""
        self._version += 1
        self._metric_cache.clear()

    def _cached(self, key: str, compute):
        """
        Return a cached derived metric, recomputing after mutations.

        Parameters
        ----------
        key : str
            Cache key naming the metric.
        compute : callable
            Zero-argument function producing the value.

        Returns
        -------
        Any
            The cached or freshly computed value.
        """
        if key not in self._metric_cache:
            self._metric_cache[key] = compute()
        return self._metric_cache[key]

    def calculate_diameter(self) -> int:
        """
        Calculate network diameter (longest shortest path).

        Cached until the network is structurally modified.

        Returns
        -------
        int
            Network diameter, or -1 if not connected.
        """

        def compute() -> int:
            """
            Compute the diameter, ignoring the cache.

            Returns
            -------
            int
                Network diameter, or -1 when disconnected.
            """
            if not self.is_connected():
                return -1
            return nx.diameter(self._graph)

        return self._cached('diameter', compute)

    def get_shortest_path(self, source: str, target: str) -> List[str]:
        """
        Find shortest path between two nodes.

        Parameters
        ----------
        source : str
            Source node ID.
        target : str
            Target node ID.

        Returns
        -------
        List[str]
            List of node IDs in the shortest path.

        Raises
        ------
        networkx.NetworkXNoPath
            If no path exists between the nodes.
        """
        return nx.shortest_path(self._graph, source, target)

    def get_shortest_path_length(self, source: str, target: str) -> float:
        """
        Get the genetic-distance-weighted shortest path between two nodes.

        Uses each edge's 'distance' attribute as its weight, matching
        PopART's weighted path lengths (Graph::pathLength); an edge
        without a distance counts as one step.

        Parameters
        ----------
        source : str
            Source node ID.
        target : str
            Target node ID.

        Returns
        -------
        float
            Total distance along the shortest weighted path.

        Raises
        ------
        networkx.NetworkXNoPath
            If no path exists between the nodes.
        """
        return nx.shortest_path_length(
            self._graph,
            source,
            target,
            weight=lambda u, v, attrs: attrs.get('distance') or 1,
        )

    def calculate_centrality(self) -> Dict[str, float]:
        """
        Calculate betweenness centrality for all nodes.

        Cached until the network is structurally modified.

        Returns
        -------
        dict
            Dictionary mapping node ID to centrality score.
        """
        return self._cached(
            'centrality', lambda: nx.betweenness_centrality(self._graph)
        )

    def get_total_samples(self) -> int:
        """
        Get total number of samples represented in the network.

        Returns
        -------
        int
            Total sample count.
        """
        return sum(hap.frequency for hap in self._haplotype_map.values())

    def calculate_stats(self) -> NetworkStats:
        """
        Calculate comprehensive network statistics.

        Cached until the network is structurally modified.

        Returns
        -------
        NetworkStats
            NetworkStats object with network metrics.
        """
        return self._cached('stats', self._compute_stats)

    def _compute_stats(self) -> NetworkStats:
        """
        Compute network statistics (uncached).

        Returns
        -------
        NetworkStats
            NetworkStats object with network metrics.
        """
        num_nodes = self.num_nodes
        num_edges = self.num_edges
        num_haplotypes = len(self.haplotypes)
        num_median = len(self._median_vectors)
        total_samples = self.get_total_samples()

        # Calculate diameter (handle disconnected networks)
        try:
            diameter = self.calculate_diameter()
        except Exception:
            diameter = -1

        # Calculate average degree
        if num_nodes > 0:
            avg_degree = (2 * num_edges) / num_nodes
        else:
            avg_degree = 0.0

        # Count connected components
        num_components = nx.number_connected_components(self._graph)

        return NetworkStats(
            num_nodes=num_nodes,
            num_edges=num_edges,
            num_haplotypes=num_haplotypes,
            num_median_vectors=num_median,
            total_samples=total_samples,
            diameter=diameter,
            avg_degree=avg_degree,
            num_components=num_components,
        )

    def validate(self) -> None:
        """
        Validate network structure.

        Raises
        ------
        ValueError: If network is invalid
        """
        # Check all nodes have haplotypes
        for node_id in self._graph.nodes():
            if node_id not in self._haplotype_map:
                raise ValueError(f"Node '{node_id}' missing haplotype reference")

        # Check all edges have valid endpoints
        for source, target in self._graph.edges():
            if source not in self._graph:
                raise ValueError(f"Edge references non-existent source '{source}'")
            if target not in self._graph:
                raise ValueError(f"Edge references non-existent target '{target}'")

        # Check edge distances are non-negative
        for source, target, data in self._graph.edges(data=True):
            distance = data.get('distance', 0)
            if distance < 0:
                raise ValueError(
                    f'Edge ({source}, {target}) has negative distance {distance}'
                )

    def to_networkx(self) -> nx.Graph:
        """
        Get the underlying NetworkX graph.

        Returns
        -------
        nx.Graph
            NetworkX Graph object.
        """
        return self._graph.copy()

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert network to dictionary representation.

        Returns
        -------
        Dict[str, Any]
            Dictionary with network data.
        """
        nodes_data = []
        for node_id in self._graph.nodes():
            hap = self._haplotype_map[node_id]
            nodes_data.append(
                {
                    'id': node_id,
                    'frequency': hap.frequency,
                    'sequence': hap.data,
                    'median_vector': self.is_median_vector(node_id),
                    'sample_ids': hap.sample_ids,
                    'populations': list(hap.get_populations()),
                }
            )

        edges_data = []
        for source, target, data in self._graph.edges(data=True):
            edges_data.append(
                {'source': source, 'target': target, 'distance': data['distance']}
            )

        return {
            'name': self.name,
            'nodes': nodes_data,
            'edges': edges_data,
            'metadata': self.metadata,
        }

    def __len__(self) -> int:
        """Return number of nodes in network."""
        return self.num_nodes

    def __str__(self) -> str:
        """Return string representation."""
        stats = self.calculate_stats()
        return (
            f'{self.name}: {stats.num_haplotypes} haplotypes, '
            f'{stats.num_median_vectors} median vectors, '
            f'{stats.num_edges} edges'
        )

    def __repr__(self) -> str:
        """Detailed representation."""
        return f'HaplotypeNetwork(nodes={self.num_nodes}, edges={self.num_edges})'
