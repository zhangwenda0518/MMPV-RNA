"""GenBank file format reader for PyPopART."""

import gzip
from pathlib import Path
from typing import Iterator, TextIO, Union

from Bio import SeqIO

from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence


class GenBankReader:
    """
    Reader for GenBank format sequence files.

    Parameters
    ----------
    filepath : str or Path
        Path to GenBank file.
    validate : bool, default=True
        Whether to validate sequences.
    """

    def __init__(self, filepath: Union[str, Path], validate: bool = True):
        """
        Initialize GenBank reader.

        Parameters
        ----------
        filepath : str or Path
            Path to GenBank file.
        validate : bool, default=True
            Whether to validate sequences.
        """
        self.filepath = Path(filepath)
        self.validate = validate

        if not self.filepath.exists():
            raise FileNotFoundError(f'File not found: {filepath}')

    def _open_file(self) -> TextIO:
        """
        Open file handling gzip compression.

        Returns
        -------
        TextIO
            Open file handle, transparently decompressed.
        """
        if self.filepath.suffix == '.gz':
            return gzip.open(self.filepath, 'rt')
        else:
            return open(self.filepath, 'r')

    def read_sequences(self, progress_callback=None) -> Iterator[Sequence]:
        """
        Read sequences from GenBank file.

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
            for record in SeqIO.parse(handle, 'genbank'):
                # Extract metadata from annotations
                metadata = {}

                if record.annotations:
                    # Add relevant annotations as metadata
                    for key in ['organism', 'source', 'taxonomy', 'date']:
                        if key in record.annotations:
                            value = record.annotations[key]
                            if isinstance(value, list):
                                value = '; '.join(str(v) for v in value)
                            metadata[key] = str(value)

                # Add features as metadata if present
                if record.features:
                    feature_types = {f.type for f in record.features}
                    metadata['features'] = ', '.join(feature_types)

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
        Read all records into an Alignment.

        Parameters
        ----------
        progress_callback : callable, optional
            Optional callback function(current, total).

        Returns
        -------
        Alignment
            Alignment of all sequences in the file.
        """
        alignment = Alignment(list(self.read_sequences(progress_callback)))
        if self.validate:
            alignment.validate()
        return alignment
