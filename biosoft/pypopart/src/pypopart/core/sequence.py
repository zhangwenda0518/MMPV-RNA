"""DNA sequence representation and manipulation for PyPopART."""

from typing import Any, Dict, Optional


class Sequence:
    """
    Represents a DNA sequence with metadata.

    Supports IUPAC nucleotide codes including ambiguous characters.

    Parameters
    ----------
    id : str
        Sequence identifier.
    data : str
        DNA sequence string.
    metadata : Dict[str, Any], optional
        Optional metadata dictionary.
    description : str, optional
        Optional sequence description.
    """

    # IUPAC nucleotide codes (uppercase)
    VALID_CHARS = set('ACGTRYSWKMBDHVN-?')

    # Reverse complement mapping
    COMPLEMENT = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'R': 'Y',
        'Y': 'R',
        'S': 'S',
        'W': 'W',
        'K': 'M',
        'M': 'K',
        'B': 'V',
        'V': 'B',
        'D': 'H',
        'H': 'D',
        'N': 'N',
        '-': '-',
        '?': '?',
    }

    def __init__(
        self,
        id: str,
        data: str,
        metadata: Optional[Dict[str, Any]] = None,
        description: Optional[str] = None,
    ):
        """
        Initialize a sequence.

        Parameters
        ----------
        id : str
            Sequence identifier.
        data : str
            DNA sequence string.
        metadata : Dict[str, Any], optional
            Optional metadata dictionary.
        description : str, optional
            Optional sequence description.

        Raises
        ------
        ValueError
            If sequence contains invalid characters.
        """
        self.id = id
        self.data = data.strip().upper()
        self.metadata = metadata if metadata is not None else {}
        self.description = description

        # Validate sequence
        self._validate()

    def _validate(self) -> None:
        """
        Validate sequence data contains only valid IUPAC characters.

        Raises
        ------
        ValueError: If invalid characters found
        """
        invalid = set(self.data) - self.VALID_CHARS
        if invalid:
            raise ValueError(f'Invalid characters in sequence: {sorted(invalid)}')

    def __len__(self) -> int:
        """Return sequence length."""
        return len(self.data)

    def __eq__(self, other: object) -> bool:
        """Check equality based on sequence data."""
        if not isinstance(other, Sequence):
            return False
        return self.data == other.data

    def __hash__(self) -> int:
        """Hash based on sequence data for use in sets/dicts."""
        return hash(self.data)

    def __str__(self) -> str:
        """Return FASTA format string representation."""
        header = f'>{self.id}'
        if self.description:
            header += f' {self.description}'
        return f'{header}\n{self.data}'

    def __repr__(self) -> str:
        """Return detailed representation."""
        return (
            f"Sequence(id='{self.id}', length={len(self)}, data='{self.data[:20]}...')"
        )

    def reverse_complement(self) -> 'Sequence':
        """
        Calculate reverse complement of sequence.

        Returns
        -------
        Sequence
            New Sequence object with reverse complement.
        """
        rev_comp_data = ''.join(self.COMPLEMENT[base] for base in reversed(self.data))
        return Sequence(
            id=f'{self.id}_rev_comp',
            data=rev_comp_data,
            metadata=self.metadata.copy(),
            description=f'Reverse complement of {self.id}',
        )

    def gc_content(self) -> float:
        """
        Calculate GC content as fraction (0.0-1.0).

        Ignores gaps (-) and ambiguous characters (N, ?).

        Returns
        -------
        float
            GC content fraction.
        """
        # Count G, C, and S (which represents G or C)
        gc_count = sum(1 for base in self.data if base in 'GCS')
        # Count only ACGT bases (exclude gaps, ambiguous chars)
        total_count = sum(1 for base in self.data if base in 'ACGT')

        if total_count == 0:
            return 0.0
        return gc_count / total_count

    def count_gaps(self) -> int:
        """
        Count number of gap characters.

        Returns
        -------
        int
            Number of gaps.
        """
        return self.data.count('-')

    def count_ambiguous(self) -> int:
        """
        Count number of ambiguous characters (N and ?).

        Returns
        -------
        int
            Number of ambiguous characters.
        """
        return sum(1 for base in self.data if base in 'N?')

    def remove_gaps(self) -> 'Sequence':
        """
        Create new sequence with gaps removed.

        Returns
        -------
        Sequence
            New Sequence object without gaps.
        """
        ungapped_data = self.data.replace('-', '')
        return Sequence(
            id=self.id,
            data=ungapped_data,
            metadata=self.metadata.copy(),
            description=self.description,
        )

    def slice(self, start: int, end: Optional[int] = None) -> 'Sequence':
        """
        Extract a slice of the sequence.

        Parameters
        ----------
        start : int
            Start position (0-based, inclusive).
        end : int, optional
            End position (0-based, exclusive), None for end of sequence.

        Returns
        -------
        Sequence
            New Sequence object with sliced data.
        """
        sliced_data = self.data[start:end]
        slice_desc = f'slice[{start}:{end if end else "end"}]'
        return Sequence(
            id=f'{self.id}_{slice_desc}',
            data=sliced_data,
            metadata=self.metadata.copy(),
            description=f'Slice {slice_desc} of {self.id}',
        )

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert sequence to dictionary representation.

        Returns
        -------
        Dict[str, Any]
            Dictionary with sequence data and statistics.
        """
        return {
            'id': self.id,
            'data': self.data,
            'metadata': self.metadata,
            'description': self.description,
            'length': len(self),
            'gc_content': self.gc_content(),
            'gaps': self.count_gaps(),
            'ambiguous': self.count_ambiguous(),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Sequence':
        """
        Create sequence from dictionary.

        Parameters
        ----------
        data : Dict[str, Any]
            Dictionary with sequence information.

        Returns
        -------
        Sequence
            New Sequence object.
        """
        return cls(
            id=data['id'],
            data=data['data'],
            metadata=data.get('metadata', {}),
            description=data.get('description'),
        )
