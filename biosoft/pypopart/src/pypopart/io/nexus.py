"""
NEXUS file format reader and writer for PyPopART.

Supports PopART-style NEXUS files with traits blocks.
"""

import gzip
import io
from pathlib import Path
import re
from typing import TYPE_CHECKING, Dict, Optional, TextIO, Tuple, Union

from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence

if TYPE_CHECKING:
    from pypopart.core.graph import HaplotypeNetwork


class NexusReader:
    """
    Reader for NEXUS format files.

    Supports PopART-style NEXUS with traits blocks.

    Parameters
    ----------
    filepath : str or Path
        Path to NEXUS file.
    validate : bool, default=True
        Whether to validate sequences and alignment.
    """

    def __init__(self, filepath: Union[str, Path], validate: bool = True):
        """
        Initialize NEXUS reader.

        Parameters
        ----------
        filepath : str or Path
            Path to NEXUS file.
        validate : bool, default=True
            Whether to validate sequences and alignment.
        """
        self.filepath = Path(filepath)
        self.validate = validate

        if not self.filepath.exists():
            raise FileNotFoundError(f'File not found: {filepath}')

        self.traits = {}

    @classmethod
    def from_string(cls, text: str, validate: bool = True) -> 'NexusReader':
        """
        Create a reader over in-memory text instead of a file.

        Lets callers (e.g. the GUI handling uploads) parse content without
        writing a temporary file to disk.

        Parameters
        ----------
        text : str
            Complete file content in this reader's format.
        validate : bool, default=True
            Whether to validate sequences.

        Returns
        -------
        NexusReader
            Reader that parses the given text.
        """
        reader = cls.__new__(cls)
        reader.filepath = None
        reader.validate = validate
        reader.traits = {}
        reader._text = text
        return reader

    def _open_file(self) -> TextIO:
        """
        Open file handling gzip compression.

        Returns
        -------
        TextIO
            Open file handle, transparently decompressed.
        """
        if getattr(self, '_text', None) is not None:
            return io.StringIO(self._text)
        if self.filepath.suffix == '.gz':
            return gzip.open(self.filepath, 'rt')
        else:
            return open(self.filepath, 'r')

    def _parse_dimensions(self, content: str) -> Tuple[int, int]:
        """
        Parse DIMENSIONS block.

        Parameters
        ----------
        content : str
            NEXUS file content.

        Returns
        -------
        Tuple[int, int]
            Tuple of (ntax, nchar).
        """
        dimensions_match = re.search(
            r'DIMENSIONS\s+NTAX=(\d+)\s+NCHAR=(\d+)', content, re.IGNORECASE
        )

        if dimensions_match:
            ntax = int(dimensions_match.group(1))
            nchar = int(dimensions_match.group(2))
            return ntax, nchar

        return 0, 0

    def _parse_matrix(self, content: str) -> Dict[str, str]:
        """
        Parse MATRIX block.

        Parameters
        ----------
        content : str
            NEXUS file content.

        Returns
        -------
        Dict[str, str]
            Dictionary mapping sequence IDs to sequence data.
        """
        sequences = {}

        # Find MATRIX block
        matrix_match = re.search(
            r'MATRIX\s+(.*?);\s*END;', content, re.IGNORECASE | re.DOTALL
        )

        if not matrix_match:
            raise ValueError('No MATRIX block found in NEXUS file')

        matrix_content = matrix_match.group(1)

        # Parse sequences (handle both interleaved and sequential formats)
        current_id = None

        for line in matrix_content.split('\n'):
            line = line.strip()
            if not line:
                continue

            # Check if line starts with sequence ID
            parts = line.split(None, 1)
            if len(parts) == 2:
                seq_id, seq_data = parts
                seq_data = seq_data.replace(' ', '').replace('\t', '')

                if seq_id in sequences:
                    # Interleaved format - append to existing sequence
                    sequences[seq_id] += seq_data
                else:
                    # New sequence
                    sequences[seq_id] = seq_data
            elif len(parts) == 1:
                # Continuation of previous sequence
                seq_data = parts[0].replace(' ', '').replace('\t', '')
                if current_id:
                    sequences[current_id] += seq_data

        return sequences

    def _parse_traits(self, content: str) -> Dict[str, Dict[str, str]]:
        """
        Parse TRAITS block (PopART extension).

        Parameters
        ----------
        content : str
            NEXUS file content.

        Returns
        -------
        Dict[str, Dict[str, str]]
            Dictionary mapping sequence IDs to trait dictionaries.
        """
        traits = {}

        # Find TRAITS block
        traits_match = re.search(
            r'BEGIN TRAITS;(.*?)END;', content, re.IGNORECASE | re.DOTALL
        )

        if not traits_match:
            return traits

        traits_content = traits_match.group(1)

        # Parse TRAITLABELS
        labels_match = re.search(r'TRAITLABELS\s+(.*?);', traits_content, re.IGNORECASE)

        trait_labels = []
        if labels_match:
            trait_labels = labels_match.group(1).split()

        # Parse MATRIX
        matrix_match = re.search(
            r'MATRIX\s+(.*?);', traits_content, re.IGNORECASE | re.DOTALL
        )

        if matrix_match:
            matrix_content = matrix_match.group(1)
            for line in matrix_content.split('\n'):
                line = line.strip()
                if not line:
                    continue

                parts = line.split()
                if len(parts) >= 2:
                    seq_id = parts[0]
                    trait_values = parts[1:]

                    traits[seq_id] = {}
                    for i, value in enumerate(trait_values):
                        label = (
                            trait_labels[i] if i < len(trait_labels) else f'trait_{i}'
                        )
                        traits[seq_id][label] = value

        return traits

    def read_alignment(self, progress_callback=None) -> Alignment:
        """
        Read alignment from NEXUS file.

        Parameters
        ----------
        progress_callback : callable, optional
            Optional callback function(current, total).

        Returns
        -------
        Alignment
            Alignment object with metadata.
        """
        with self._open_file() as handle:
            content = handle.read()

        # Parse dimensions
        ntax, nchar = self._parse_dimensions(content)

        # Parse sequences
        seq_dict = self._parse_matrix(content)

        # Parse traits
        self.traits = self._parse_traits(content)

        # Create Sequence objects
        sequences = []
        count = 0

        for seq_id, seq_data in seq_dict.items():
            metadata = self.traits.get(seq_id, {})

            seq = Sequence(id=seq_id, data=seq_data.upper(), metadata=metadata)

            # Validation happens automatically in Sequence.__init__

            sequences.append(seq)

            count += 1
            if progress_callback:
                progress_callback(count, ntax)

        alignment = Alignment(sequences)

        if self.validate:
            alignment.validate()

        return alignment

    def get_traits(self) -> Dict[str, Dict[str, str]]:
        """
        Get parsed traits/metadata.

        Returns
        -------
        Dict[str, Dict[str, str]]
            Dictionary mapping sequence IDs to trait dictionaries.
        """
        return self.traits


class NexusWriter:
    """
    Writer for NEXUS format files.

    Supports PopART-style NEXUS with traits blocks.

    Parameters
    ----------
    filepath : str or Path
        Output file path.
    interleaved : bool, default=False
        Whether to write in interleaved format.
    compress : str, optional
        Compression format ('gzip' or None).
    """

    def __init__(
        self,
        filepath: Union[str, Path],
        interleaved: bool = False,
        compress: Optional[str] = None,
    ):
        """
        Initialize NEXUS writer.

        Parameters
        ----------
        filepath : str or Path
            Output file path.
        interleaved : bool, default=False
            Whether to write in interleaved format.
        compress : str, optional
            Compression format ('gzip' or None).
        """
        self.filepath = Path(filepath)
        self.interleaved = interleaved
        self.compress = compress

        if compress == 'gzip' and not str(self.filepath).endswith('.gz'):
            self.filepath = Path(str(self.filepath) + '.gz')

    def _open_file(self) -> TextIO:
        """
        Open file for writing with optional compression.

        Returns
        -------
        TextIO
            Open file handle for writing, optionally compressed.
        """
        if self.compress == 'gzip':
            return gzip.open(self.filepath, 'wt')
        else:
            return open(self.filepath, 'w')

    def write_alignment(
        self, alignment: Alignment, include_traits: bool = True, progress_callback=None
    ) -> None:
        """
        Write alignment to NEXUS file.

        Parameters
        ----------
        alignment : Alignment
            Alignment object.
        include_traits : bool, default=True
            Whether to include traits block.
        progress_callback : callable, optional
            Optional callback function(current, total).
        """
        with self._open_file() as handle:
            # Write header
            handle.write('#NEXUS\n\n')

            # Write DATA block
            handle.write('BEGIN DATA;\n')
            handle.write(
                f'  DIMENSIONS NTAX={len(alignment)} NCHAR={alignment.length};\n'
            )
            handle.write('  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n')
            handle.write('  MATRIX\n')

            # Write sequences
            count = 0
            for seq in alignment:
                # Pad ID to align sequences
                padded_id = seq.id.ljust(20)
                handle.write(f'    {padded_id} {seq.data}\n')

                count += 1
                if progress_callback:
                    progress_callback(count, len(alignment))

            handle.write('  ;\n')
            handle.write('END;\n')

            # Write TRAITS block if requested and metadata exists
            if include_traits:
                # Collect all trait keys
                all_traits = set()
                for seq in alignment:
                    all_traits.update(seq.metadata.keys())

                if all_traits:
                    handle.write('\n')
                    handle.write('BEGIN TRAITS;\n')
                    handle.write('  DIMENSIONS NTRAITS={};\n'.format(len(all_traits)))

                    # Write trait labels
                    trait_list = sorted(all_traits)
                    handle.write('  TRAITLABELS {};\n'.format(' '.join(trait_list)))

                    # Write trait matrix
                    handle.write('  MATRIX\n')
                    for seq in alignment:
                        padded_id = seq.id.ljust(20)
                        trait_values = [
                            str(seq.metadata.get(trait, '?')) for trait in trait_list
                        ]
                        handle.write(f'    {padded_id} {" ".join(trait_values)}\n')

                    handle.write('  ;\n')
                    handle.write('END;\n')

    def write_network(
        self,
        network: 'HaplotypeNetwork',
        layout: Optional[Dict[str, Tuple[float, float]]] = None,
    ) -> None:
        """
        Write a haplotype network as a PopART-style NEXUS file.

        Produces a DATA block holding the sampled haplotype sequences and
        a NETWORK block in the layout PopART itself saves: a TRANSLATE
        table of vertex labels, VERTICES with coordinates, VLABELS with
        label offsets and EDGES with mutation counts. Median vectors have
        no sequence and appear only as vertices.

        Parameters
        ----------
        network : HaplotypeNetwork
            Network to write.
        layout : Dict[str, Tuple[float, float]], optional
            Node positions. Computed with a seeded spring layout when
            omitted, so the same network always writes the same file.
        """
        import networkx as nx

        graph = network.graph
        nodes = list(graph.nodes())
        index = {node: i for i, node in enumerate(nodes)}
        if layout is None:
            layout = nx.spring_layout(graph, seed=42, scale=200.0)

        sampled = [n for n in nodes if not network.is_median_vector(n)]
        sequences = {
            n: str(graph.nodes[n].get('sequence') or '')
            for n in sampled
            if graph.nodes[n].get('sequence')
        }

        with self._open_file() as handle:
            handle.write('#NEXUS\n\n')

            if sequences:
                length = max(len(seq) for seq in sequences.values())
                handle.write('BEGIN DATA;\n')
                handle.write(f'  DIMENSIONS NTAX={len(sequences)} NCHAR={length};\n')
                handle.write('  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n')
                handle.write('  MATRIX\n')
                for node, seq in sequences.items():
                    handle.write(f'    {node.ljust(20)} {seq}\n')
                handle.write('  ;\nEND;\n\n')

            handle.write('BEGIN NETWORK;\n')
            handle.write(
                f'  DIMENSIONS ntax={len(sampled)} nvertices={len(nodes)} '
                f'nedges={graph.number_of_edges()};\n'
            )
            handle.write('  FORMAT VSize=10 EView=Numbers;\n')

            handle.write('  TRANSLATE\n')
            for i, node in enumerate(sampled, start=1):
                handle.write(f'    {i} {node},\n')
            handle.write('  ;\n')

            handle.write('  VERTICES\n')
            for i, node in enumerate(nodes, start=1):
                x, y = layout.get(node, (0.0, 0.0))
                handle.write(f'    {i} {float(x):.4f} {float(y):.4f},\n')
            handle.write('  ;\n')

            handle.write('  VLABELS\n')
            for i, node in enumerate(nodes, start=1):
                x, y = layout.get(node, (0.0, 0.0))
                handle.write(f'    {i} {float(x) + 5.0:.4f} {float(y) - 5.0:.4f},\n')
            handle.write('  ;\n')

            # PopART lists vertices 1-based but refers to them 0-based in
            # the edge table; match that so PopART can read the file.
            handle.write('  EDGES\n')
            for i, (u, v) in enumerate(graph.edges(), start=1):
                weight = graph[u][v].get('distance', graph[u][v].get('weight', 1))
                handle.write(f'    {i} {index[u]} {index[v]} {int(weight)},\n')
            handle.write('  ;\n')

            handle.write('END;\n')
