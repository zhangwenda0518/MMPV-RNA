"""
File I/O module for PyPopART.

Provides readers and writers for various sequence and network formats.
"""

from pathlib import Path
from typing import Optional, Union

import networkx as nx

from .fasta import FastaReader, FastaWriter
from .genbank import GenBankReader
from .metadata import MetadataReader, MetadataWriter
from .network_export import (
    CSVExporter,
    CytoscapeExporter,
    GMLExporter,
    GraphMLExporter,
    JSONExporter,
)
from .nexus import NexusReader, NexusWriter
from .phylip import PhylipReader, PhylipWriter


def load_alignment(filepath: Union[str, Path], format: Optional[str] = None):
    """
    Load sequence alignment from file.

    Auto-detects format if not specified based on file extension.

    Parameters
    ----------
    filepath : str or Path
        Path to sequence file.
    format : str, optional
        File format: 'fasta', 'nexus', 'phylip', 'genbank'
        If None, auto-detect from extension.

    Returns
    -------
    Alignment
        Alignment object with the file's sequences.
    """
    filepath = Path(filepath)

    # Auto-detect format
    if format is None:
        suffix = filepath.suffix.lower()
        if suffix in ['.fasta', '.fa', '.fna']:
            format = 'fasta'
        elif suffix in ['.nexus', '.nex']:
            format = 'nexus'
        elif suffix in ['.phy', '.phylip']:
            format = 'phylip'
        elif suffix in ['.gb', '.gbk']:
            format = 'genbank'
        else:
            # Try FASTA as default
            format = 'fasta'

    # Load with appropriate reader
    if format.lower() == 'fasta':
        reader = FastaReader(filepath)
    elif format.lower() == 'nexus':
        reader = NexusReader(filepath)
    elif format.lower() == 'phylip':
        reader = PhylipReader(filepath)
    elif format.lower() == 'genbank':
        reader = GenBankReader(filepath)
    else:
        raise ValueError(f'Unknown format: {format}')

    return reader.read_alignment()


def save_alignment(alignment, filepath: Union[str, Path], format: str = 'fasta'):
    """
    Save alignment to file.

    Parameters
    ----------
    alignment : Alignment
        Alignment to save.
    filepath : str or Path
        Output file path.
    format : str, default='fasta'
        Output format: 'fasta', 'nexus', 'phylip'.
    """
    filepath = Path(filepath)

    if format.lower() == 'fasta':
        writer = FastaWriter(filepath)
    elif format.lower() == 'nexus':
        writer = NexusWriter(filepath)
    elif format.lower() == 'phylip':
        writer = PhylipWriter(filepath)
    else:
        raise ValueError(f'Unknown format: {format}')

    writer.write_alignment(alignment)


def load_network(filepath: Union[str, Path], format: Optional[str] = None):
    """
    Load network from file.

    Parameters
    ----------
    filepath : str or Path
        Path to network file.
    format : str, optional
        File format: 'graphml', 'gml', 'json'
        If None, auto-detect from extension.

    Returns
    -------
    HaplotypeNetwork
        Reconstructed network object.
    """
    from pypopart.core.graph import HaplotypeNetwork

    filepath = Path(filepath)

    # Auto-detect format
    if format is None:
        suffix = filepath.suffix.lower()
        if suffix == '.graphml':
            format = 'graphml'
        elif suffix == '.gml':
            format = 'gml'
        elif suffix == '.json':
            format = 'json'
        else:
            format = 'graphml'

    # Load with appropriate method
    if format.lower() == 'graphml':
        graph = nx.read_graphml(str(filepath))
    elif format.lower() == 'gml':
        graph = nx.read_gml(str(filepath))
    elif format.lower() == 'json':
        import json

        with open(filepath) as f:
            data = json.load(f)
        if 'links' in data:
            # NetworkX node-link format
            graph = nx.node_link_graph(data)
        else:
            # PyPopART JSONExporter format: nodes/edges with 'attributes'
            graph = nx.Graph()
            graph.graph.update(data.get('metadata', {}))
            for node in data.get('nodes', []):
                graph.add_node(node['id'], **node.get('attributes', {}))
            for edge in data.get('edges', []):
                graph.add_edge(
                    edge['source'], edge['target'], **edge.get('attributes', {})
                )
    else:
        raise ValueError(f'Unknown format: {format}')

    return HaplotypeNetwork.from_networkx(graph)


def save_network(
    network: nx.Graph, filepath: Union[str, Path], format: str = 'graphml'
):
    """
    Save network to file.

    Parameters
    ----------
    network : nx.Graph
        Network to save.
    filepath : str or Path
        Output file path.
    format : str, default='graphml'
        Output format: 'graphml', 'gml', 'json', 'nexus', 'cytoscape', 'csv'.
    """
    filepath = Path(filepath)

    if format.lower() == 'graphml':
        exporter = GraphMLExporter(filepath)
        exporter.export(network)
    elif format.lower() == 'gml':
        exporter = GMLExporter(filepath)
        exporter.export(network)
    elif format.lower() == 'json':
        exporter = JSONExporter(filepath)
        exporter.export(network)
    elif format.lower() == 'nexus':
        writer = NexusWriter(filepath)
        writer.write_network(network)
    elif format.lower() == 'cytoscape':
        exporter = CytoscapeExporter(filepath)
        exporter.export(network)
    elif format.lower() == 'csv':
        exporter = CSVExporter(filepath)
        exporter.export(network)
    else:
        raise ValueError(f'Unknown format: {format}')


__all__ = [
    'FastaReader',
    'FastaWriter',
    'NexusReader',
    'NexusWriter',
    'PhylipReader',
    'PhylipWriter',
    'GenBankReader',
    'MetadataReader',
    'MetadataWriter',
    'GraphMLExporter',
    'GMLExporter',
    'CytoscapeExporter',
    'JSONExporter',
    'CSVExporter',
    'load_alignment',
    'save_alignment',
    'load_network',
    'save_network',
]
