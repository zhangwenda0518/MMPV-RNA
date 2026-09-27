"""FASTA file format reader and writer for PyPopART."""

import gzip
import io
from pathlib import Path
from typing import Iterator, Optional, TextIO, Union
import zipfile

from Bio import SeqIO

from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence


class FastaReader:
    """
    Reader for FASTA format sequence files.

    Supports plain text, gzip, and zip compressed files.

    Parameters
    ----------
    filepath : str or Path
        Path to FASTA file.
    validate : bool, default=True
        Whether to validate sequences.
    """

    def __init__(self, filepath: Union[str, Path], validate: bool = True):
        """
        Initialize FASTA reader.

        Parameters
        ----------
        filepath : str or Path
            Path to FASTA file.
        validate : bool, default=True
            Whether to validate sequences.
        """
        self.filepath = Path(filepath)
        self.validate = validate

        if not self.filepath.exists():
            raise FileNotFoundError(f'File not found: {filepath}')

    @classmethod
    def from_string(cls, text: str, validate: bool = True) -> 'FastaReader':
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
        FastaReader
            Reader that parses the given text.
        """
        reader = cls.__new__(cls)
        reader.filepath = None
        reader.validate = validate
        reader._text = text
        return reader

    def _open_file(self) -> TextIO:
        """
        Open file handling compression automatically.

        Returns
        -------
        TextIO
            File handle.
        """
        if getattr(self, '_text', None) is not None:
            return io.StringIO(self._text)
        if self.filepath.suffix == '.gz':
            return gzip.open(self.filepath, 'rt')
        elif self.filepath.suffix == '.zip':
            with zipfile.ZipFile(self.filepath) as zf:
                # Assume first file in zip is the FASTA
                names = zf.namelist()
                if not names:
                    raise ValueError('Empty zip file')
                return zf.open(names[0], 'r')
        else:
            return open(self.filepath, 'r')

    def read_sequences(self, progress_callback=None) -> Iterator[Sequence]:
        """
        Read sequences from FASTA file.

        Parameters
        ----------
        progress_callback : callable, optional
            Optional callback function(current, total).

        Yields
        ------
        Sequence
            Sequence objects.
        """
        count = 0
        with self._open_file() as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                # Extract metadata from description
                metadata = {}
                if '|' in record.description:
                    # Parse pipe-separated metadata
                    parts = record.description.split('|')
                    if len(parts) > 1:
                        for part in parts[1:]:
                            if '=' in part:
                                key, value = part.split('=', 1)
                                metadata[key.strip()] = value.strip()

                seq = Sequence(
                    id=record.id,
                    data=str(record.seq).upper(),
                    description=record.description,
                    metadata=metadata,
                )

                # Validation happens automatically in Sequence.__init__

                count += 1
                if progress_callback:
                    progress_callback(count, None)

                yield seq

    def read_alignment(self, progress_callback=None) -> Alignment:
        """
        Read alignment from FASTA file.

        Parameters
        ----------
        progress_callback : callable, optional
            Optional callback function(current, total).

        Returns
        -------
        Alignment
            Alignment object.
        """
        sequences = list(self.read_sequences(progress_callback))
        alignment = Alignment(sequences)

        if self.validate:
            alignment.validate()

        return alignment


class FastaWriter:
    """
    Writer for FASTA format sequence files.

    Parameters
    ----------
    filepath : str or Path
        Output file path.
    line_length : int, default=80
        Maximum line length for sequences (0 for no wrapping).
    compress : str, optional
        Compression format ('gzip' or None).
    """

    def __init__(
        self,
        filepath: Union[str, Path],
        line_length: int = 80,
        compress: Optional[str] = None,
    ):
        """
        Initialize FASTA writer.

        Parameters
        ----------
        filepath : str or Path
            Output file path.
        line_length : int, default=80
            Maximum line length for sequences (0 for no wrapping).
        compress : str, optional
            Compression format ('gzip' or None).
        """
        self.filepath = Path(filepath)
        self.line_length = line_length
        self.compress = compress

        # Add compression extension if needed
        if compress == 'gzip' and not str(self.filepath).endswith('.gz'):
            self.filepath = Path(str(self.filepath) + '.gz')

    def _open_file(self) -> TextIO:
        """
        Open file for writing with optional compression.

        Returns
        -------
        TextIO
            File handle.
        """
        if self.compress == 'gzip':
            return gzip.open(self.filepath, 'wt')
        else:
            return open(self.filepath, 'w')

    def write_sequences(
        self, sequences: Iterator[Sequence], progress_callback=None
    ) -> int:
        """
        Write sequences to FASTA file.

        Parameters
        ----------
        sequences : Iterator[Sequence]
            Iterable of Sequence objects.
        progress_callback : callable, optional
            Optional callback function(current, total).

        Returns
        -------
        int
            Number of sequences written.
        """
        count = 0

        with self._open_file() as handle:
            for seq in sequences:
                # Build header
                header = f'>{seq.id}'
                if seq.description and seq.description != seq.id:
                    header = f'>{seq.description}'

                # Add metadata to header
                if seq.metadata:
                    metadata_str = '|'.join(f'{k}={v}' for k, v in seq.metadata.items())
                    header += f'|{metadata_str}'

                handle.write(header + '\n')

                # Write sequence data with line wrapping
                if self.line_length > 0:
                    for i in range(0, len(seq.data), self.line_length):
                        handle.write(seq.data[i : i + self.line_length] + '\n')
                else:
                    handle.write(seq.data + '\n')

                count += 1
                if progress_callback:
                    progress_callback(count, None)

        return count

    def write_alignment(self, alignment: Alignment, progress_callback=None) -> int:
        """
        Write alignment to FASTA file.

        Parameters
        ----------
        alignment : Alignment
            Alignment object.
        progress_callback : callable, optional
            Optional callback function(current, total).

        Returns
        -------
        int
            Number of sequences written.
        """
        return self.write_sequences(iter(alignment), progress_callback)
