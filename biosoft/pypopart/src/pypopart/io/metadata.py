"""
Metadata file reader and writer for PyPopART.

Handles CSV-based metadata/traits files that can be linked to sequences.
Includes support for geographic coordinates (latitude/longitude).
"""

import csv
import gzip
import io
from pathlib import Path
from typing import Dict, List, Optional, TextIO, Tuple, Union

from pypopart.core.alignment import Alignment


def parse_coordinate(value: str) -> float:
    """
    Parse a coordinate string to a float.

    Handles various formats:
    - Decimal degrees: "45.5", "-123.4"
    - With degree symbol: "45.5°", "-123.4°"

    Parameters
    ----------
    value : str
        Coordinate string.

    Returns
    -------
    float
        Coordinate as float.

    Raises
    ------
    ValueError: If coordinate cannot be parsed
    """
    try:
        # Remove degree symbol and whitespace
        cleaned = value.strip().replace('°', '').replace(' ', '')
        return float(cleaned)
    except (ValueError, AttributeError) as e:
        raise ValueError(f'Invalid coordinate value: {value}') from e


def validate_latitude(lat: float) -> None:
    """
    Validate latitude value.

    Parameters
    ----------
    lat : float
        Latitude value.

    Raises
    ------
    ValueError: If latitude is out of range [-90, 90]
    """
    if not -90 <= lat <= 90:
        raise ValueError(f'Latitude must be between -90 and 90, got {lat}')


def validate_longitude(lon: float) -> None:
    """
    Validate longitude value.

    Parameters
    ----------
    lon : float
        Longitude value.

    Raises
    ------
    ValueError: If longitude is out of range [-180, 180]
    """
    if not -180 <= lon <= 180:
        raise ValueError(f'Longitude must be between -180 and 180, got {lon}')


def extract_coordinates(
    metadata: Dict[str, str],
    lat_column: str = 'latitude',
    lon_column: str = 'longitude',
    validate: bool = True,
) -> Optional[Tuple[float, float]]:
    """
    Extract and validate geographic coordinates from metadata.

    Parameters
    ----------
    metadata : Dict[str, str]
        Metadata dictionary.
    lat_column : str, default='latitude'
        Name of latitude column.
    lon_column : str, default='longitude'
        Name of longitude column.
    validate : bool, default=True
        Whether to validate coordinate ranges.

    Returns
    -------
    Tuple[float, float], optional
        Tuple of (latitude, longitude) or None if coordinates not present.

    Raises
    ------
    ValueError: If coordinates are invalid
    """
    if lat_column not in metadata or lon_column not in metadata:
        return None

    lat = parse_coordinate(metadata[lat_column])
    lon = parse_coordinate(metadata[lon_column])

    if validate:
        validate_latitude(lat)
        validate_longitude(lon)

    return (lat, lon)


class MetadataReader:
    """
    Reader for CSV-based metadata files.

    Parameters
    ----------
    filepath : str or Path
        Path to metadata CSV file.
    id_column : str, default='id'
        Name of column containing sequence IDs.
    delimiter : str, default=','
        CSV delimiter character.
    validate : bool, default=True
        Whether to validate metadata.
    """

    def __init__(
        self,
        filepath: Union[str, Path],
        id_column: str = 'id',
        delimiter: str = ',',
        validate: bool = True,
    ):
        """
        Initialize metadata reader.

        Parameters
        ----------
        filepath : str or Path
            Path to metadata CSV file.
        id_column : str, default='id'
            Name of column containing sequence IDs.
        delimiter : str, default=','
            CSV delimiter character.
        validate : bool, default=True
            Whether to validate metadata.
        """
        self.filepath = Path(filepath)
        self.id_column = id_column
        self.delimiter = delimiter
        self.validate = validate

        if not self.filepath.exists():
            raise FileNotFoundError(f'File not found: {filepath}')

    @classmethod
    def from_string(
        cls,
        text: str,
        id_column: str = 'id',
        delimiter: str = ',',
        validate: bool = True,
    ) -> 'MetadataReader':
        """
        Create a reader over in-memory CSV text instead of a file.

        Parameters
        ----------
        text : str
            Complete CSV content.
        id_column : str, default='id'
            Name of the sequence-id column.
        delimiter : str, default=','
            Field delimiter.
        validate : bool, default=True
            Whether to validate for duplicate ids.

        Returns
        -------
        MetadataReader
            Reader that parses the given text.
        """
        reader = cls.__new__(cls)
        reader.filepath = None
        reader.id_column = id_column
        reader.delimiter = delimiter
        reader.validate = validate
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
            return open(self.filepath, 'r', encoding='utf-8')

    def read_metadata(self) -> Dict[str, Dict[str, str]]:
        """
        Read metadata from CSV file.

        Returns
        -------
        Dict[str, Dict[str, str]]
            Dictionary mapping sequence IDs to metadata dictionaries.
        """
        metadata = {}

        with self._open_file() as handle:
            reader = csv.DictReader(handle, delimiter=self.delimiter)

            if self.id_column not in reader.fieldnames:
                raise ValueError(
                    f"ID column '{self.id_column}' not found in CSV. "
                    f'Available columns: {reader.fieldnames}'
                )

            for row in reader:
                seq_id = row[self.id_column]

                # Remove ID column from metadata
                row_metadata = {k: v for k, v in row.items() if k != self.id_column}

                if seq_id in metadata and self.validate:
                    raise ValueError(f'Duplicate sequence ID in metadata: {seq_id}')

                metadata[seq_id] = row_metadata

        return metadata

    def apply_to_alignment(self, alignment: Alignment) -> None:
        """
        Apply metadata to sequences in alignment.

        Parameters
        ----------
        alignment : Alignment
            Alignment object to update.
        """
        metadata = self.read_metadata()

        missing_ids = []
        for seq in alignment:
            if seq.id in metadata:
                # Update sequence metadata
                seq.metadata.update(metadata[seq.id])
            elif self.validate:
                missing_ids.append(seq.id)

        if missing_ids and self.validate:
            raise ValueError(
                f'Metadata missing for {len(missing_ids)} sequences: '
                f'{", ".join(missing_ids[:5])}{"..." if len(missing_ids) > 5 else ""}'
            )


class MetadataWriter:
    """
    Writer for CSV-based metadata files.

    Parameters
    ----------
    filepath : str or Path
        Output file path.
    id_column : str, default='id'
        Name of column for sequence IDs.
    delimiter : str, default=','
        CSV delimiter character.
    compress : str, optional
        Compression format ('gzip' or None).
    """

    def __init__(
        self,
        filepath: Union[str, Path],
        id_column: str = 'id',
        delimiter: str = ',',
        compress: Optional[str] = None,
    ):
        """
        Initialize metadata writer.

        Parameters
        ----------
        filepath : str or Path
            Output file path.
        id_column : str, default='id'
            Name of column for sequence IDs.
        delimiter : str, default=','
            CSV delimiter character.
        compress : str, optional
            Compression format ('gzip' or None).
        """
        self.filepath = Path(filepath)
        self.id_column = id_column
        self.delimiter = delimiter
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
            return gzip.open(self.filepath, 'wt', encoding='utf-8')
        else:
            return open(self.filepath, 'w', encoding='utf-8', newline='')

    def write_metadata(
        self,
        metadata: Dict[str, Dict[str, str]],
        trait_order: Optional[List[str]] = None,
    ) -> None:
        """
        Write metadata to CSV file.

        Parameters
        ----------
        metadata : Dict[str, Dict[str, str]]
            Dictionary mapping sequence IDs to metadata dictionaries.
        trait_order : List[str], optional
            Optional list specifying order of trait columns.
        """
        if not metadata:
            raise ValueError('No metadata to write')

        # Determine all trait keys
        all_traits = set()
        for traits in metadata.values():
            all_traits.update(traits.keys())

        # Order traits
        if trait_order:
            # Use specified order, then add any remaining traits
            fieldnames = [self.id_column] + trait_order
            remaining = sorted(all_traits - set(trait_order))
            fieldnames.extend(remaining)
        else:
            fieldnames = [self.id_column] + sorted(all_traits)

        with self._open_file() as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
                delimiter=self.delimiter,
                extrasaction='ignore',
            )

            writer.writeheader()

            for seq_id, traits in sorted(metadata.items()):
                row = {self.id_column: seq_id}
                row.update(traits)
                writer.writerow(row)

    def write_from_alignment(
        self, alignment: Alignment, trait_order: Optional[List[str]] = None
    ) -> None:
        """
        Extract and write metadata from alignment.

        Parameters
        ----------
        alignment : Alignment
            Alignment object.
        trait_order : List[str], optional
            Optional list specifying order of trait columns.
        """
        metadata = {}

        for seq in alignment:
            if seq.metadata:
                metadata[seq.id] = seq.metadata

        if not metadata:
            raise ValueError('No metadata found in alignment')

        self.write_metadata(metadata, trait_order)
