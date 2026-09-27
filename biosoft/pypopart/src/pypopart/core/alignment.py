"""Multiple sequence alignment representation and analysis for PyPopART."""

from collections import Counter
from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional, Union

import numpy as np

from .sequence import Sequence


@dataclass
class AlignmentStats:
    """Statistics about a multiple sequence alignment."""

    length: int
    num_sequences: int
    total_gaps: int
    gap_percentage: float
    conserved_sites: int
    variable_sites: int
    parsimony_informative_sites: int
    gc_content: float


class Alignment:
    """
    Represents a multiple sequence alignment.

    Provides methods for alignment analysis, validation, and manipulation.

    Parameters
    ----------
    sequences : List[Sequence], optional
        List of Sequence objects.
    """

    def __init__(self, sequences: Optional[List[Sequence]] = None):
        """
        Initialize alignment with optional sequences.

        Parameters
        ----------
        sequences : List[Sequence], optional
            List of Sequence objects.
        """
        self._sequences: List[Sequence] = []
        self._sequence_index: Dict[str, int] = {}

        if sequences:
            for seq in sequences:
                self.add_sequence(seq)

    def add_sequence(self, sequence: Sequence) -> None:
        """
        Add a sequence to the alignment.

        Parameters
        ----------
        sequence : Sequence
            Sequence object to add.

        Raises
        ------
        ValueError
            If sequence ID already exists or length doesn't match.
        """
        if sequence.id in self._sequence_index:
            raise ValueError(f"Sequence ID '{sequence.id}' already exists in alignment")

        # Check length consistency (if not first sequence)
        if self._sequences and len(sequence) != self.length:
            raise ValueError(
                f"Sequence length {len(sequence)} doesn't match alignment length {self.length}"
            )

        self._sequence_index[sequence.id] = len(self._sequences)
        self._sequences.append(sequence)

    def remove_sequence(self, sequence_id: str) -> None:
        """
        Remove a sequence from the alignment.

        Parameters
        ----------
        sequence_id : str
            ID of sequence to remove.

        Raises
        ------
        KeyError
            If sequence ID not found.
        """
        if sequence_id not in self._sequence_index:
            raise KeyError(f"Sequence ID '{sequence_id}' not found in alignment")

        index = self._sequence_index[sequence_id]
        del self._sequences[index]

        # Rebuild index
        self._sequence_index = {seq.id: i for i, seq in enumerate(self._sequences)}

    def get_sequence(self, sequence_id: str) -> Sequence:
        """
        Get sequence by ID.

        Parameters
        ----------
        sequence_id : str
            ID of sequence to retrieve.

        Returns
        -------
        Sequence
            Sequence object.
            Raises :
            KeyError :
            If sequence ID not found.
        """
        if sequence_id not in self._sequence_index:
            raise KeyError(f"Sequence ID '{sequence_id}' not found in alignment")

        index = self._sequence_index[sequence_id]
        return self._sequences[index]

    def __len__(self) -> int:
        """Return number of sequences in alignment."""
        return len(self._sequences)

    def __iter__(self) -> Iterator[Sequence]:
        """Iterate over sequences in alignment."""
        return iter(self._sequences)

    def __getitem__(self, key: Union[int, str, slice]) -> Union[Sequence, 'Alignment']:
        """
        Get sequence(s) by index, ID, or slice.

        Parameters
        ----------
        key : int or str or slice
            Index, sequence ID, or slice.

        Returns
        -------
        Sequence or Alignment
            Sequence object or new Alignment object for slices.
        """
        if isinstance(key, str):
            return self.get_sequence(key)
        elif isinstance(key, int):
            return self._sequences[key]
        elif isinstance(key, slice):
            return Alignment(self._sequences[key])
        else:
            raise TypeError(f'Invalid key type: {type(key)}')

    @property
    def length(self) -> int:
        """
        Get alignment length (sequence length).

        Returns
        -------
        int
            Alignment length in sites.
        """
        if not self._sequences:
            return 0
        return len(self._sequences[0])

    @property
    def sequence_ids(self) -> List[str]:
        """
        Get list of sequence IDs.

        Returns
        -------
        List[str]
            List of sequence IDs.
        """
        return [seq.id for seq in self._sequences]

    def is_valid(self) -> bool:
        """
        Check if alignment is valid (all sequences same length).

        Returns
        -------
        bool
            True if all sequences have the same length.
        """
        if not self._sequences:
            return True

        first_length = len(self._sequences[0])
        return all(len(seq) == first_length for seq in self._sequences)

    def validate(self) -> None:
        """
        Validate alignment and raise exception if invalid.

        Raises
        ------
        ValueError: If alignment is invalid
        """
        if not self.is_valid():
            lengths = [len(seq) for seq in self._sequences]
            raise ValueError(
                f'Alignment invalid: sequences have different lengths {set(lengths)}'
            )

    def get_column(self, position: int) -> List[str]:
        """
        Get all characters at a specific position.

        Parameters
        ----------
        position : int
            0-based position in alignment.

        Returns
        -------
        List[str]
            List of characters at that position.
        """
        if position < 0 or position >= self.length:
            raise IndexError(f'Position {position} out of range [0, {self.length})')

        return [seq.data[position] for seq in self._sequences]

    def slice_alignment(self, start: int, end: Optional[int] = None) -> 'Alignment':
        """
        Extract a slice of the alignment (specific columns).

        Parameters
        ----------
        start : int
            Start position (0-based, inclusive).
        end : int, optional
            End position (0-based, exclusive, None for end).

        Returns
        -------
        Alignment
            New Alignment object with sliced sequences.
        """
        sliced_sequences = []
        for seq in self._sequences:
            sliced_seq = seq.slice(start, end)
            sliced_sequences.append(sliced_seq)

        return Alignment(sliced_sequences)

    def remove_gaps_columns(self, gap_threshold: float = 1.0) -> 'Alignment':
        """
        Remove columns with gaps above threshold.

        Parameters
        ----------
        gap_threshold : float, default=1.0
            Fraction of gaps required to remove column (0.0-1.0).

        Returns
        -------
        Alignment
            New Alignment object with gap columns removed.
        """
        columns_to_keep = []

        for pos in range(self.length):
            column = self.get_column(pos)
            gap_fraction = column.count('-') / len(column)

            if gap_fraction < gap_threshold:
                columns_to_keep.append(pos)

        # Build new sequences with kept columns
        new_sequences = []
        for seq in self._sequences:
            new_data = ''.join(seq.data[pos] for pos in columns_to_keep)
            new_seq = Sequence(
                id=seq.id,
                data=new_data,
                metadata=seq.metadata.copy(),
                description=seq.description,
            )
            new_sequences.append(new_seq)

        return Alignment(new_sequences)

    def calculate_stats(self) -> AlignmentStats:
        """
        Calculate comprehensive alignment statistics.

        Returns
        -------
        AlignmentStats
            AlignmentStats object with alignment metrics.
        """
        if not self._sequences:
            return AlignmentStats(0, 0, 0, 0.0, 0, 0, 0, 0.0)

        length = self.length
        num_sequences = len(self._sequences)

        # Count gaps
        total_gaps = sum(seq.count_gaps() for seq in self._sequences)
        gap_percentage = (total_gaps / (length * num_sequences)) * 100

        # Analyze sites
        conserved_sites = 0
        variable_sites = 0
        parsimony_informative_sites = 0

        for pos in range(length):
            column = self.get_column(pos)
            # Remove gaps for site analysis
            column_no_gaps = [char for char in column if char != '-']

            if not column_no_gaps:
                continue

            char_counts = Counter(column_no_gaps)
            unique_chars = len(char_counts)

            if unique_chars == 1:
                conserved_sites += 1
            else:
                variable_sites += 1

                # Parsimony informative: at least 2 different chars, each appearing ≥2 times
                if sum(1 for count in char_counts.values() if count >= 2) >= 2:
                    parsimony_informative_sites += 1

        # Calculate overall GC content
        all_bases = ''.join(seq.data.replace('-', '') for seq in self._sequences)
        gc_count = all_bases.count('G') + all_bases.count('C') + all_bases.count('S')
        total_bases = len(all_bases.replace('N', '').replace('?', ''))
        gc_content = (gc_count / total_bases) * 100 if total_bases > 0 else 0.0

        return AlignmentStats(
            length=length,
            num_sequences=num_sequences,
            total_gaps=total_gaps,
            gap_percentage=gap_percentage,
            conserved_sites=conserved_sites,
            variable_sites=variable_sites,
            parsimony_informative_sites=parsimony_informative_sites,
            gc_content=gc_content,
        )

    def get_distance_matrix(self, distance_func=None) -> np.ndarray:
        """
        Calculate pairwise distance matrix between sequences.

        Parameters
        ----------
        distance_func : callable, optional
            Function to calculate distance between two sequences.
                          If None, uses Hamming distance.

        Returns
        -------
        np.ndarray
            Square numpy array with pairwise distances.
        """
        if distance_func is None:
            distance_func = self._hamming_distance

        n = len(self._sequences)
        matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(i + 1, n):
                dist = distance_func(self._sequences[i], self._sequences[j])
                matrix[i, j] = matrix[j, i] = dist

        return matrix

    def _hamming_distance(self, seq1: Sequence, seq2: Sequence) -> int:
        """
        Calculate Hamming distance between two sequences.

        Parameters
        ----------
        seq1 : Sequence
            First sequence.
        seq2 : Sequence
            Second sequence.

        Returns
        -------
        int
            Number of differing positions.
        """
        if len(seq1) != len(seq2):
            raise ValueError('Sequences must have same length for Hamming distance')

        return sum(
            c1 != c2 for c1, c2 in zip(seq1.data, seq2.data) if c1 != '-' and c2 != '-'
        )

    def identify_haplotypes(self) -> Dict[str, List[str]]:
        """
        Identify unique haplotypes and group sequences.

        Returns
        -------
        Dict[str, List[str]]
            Dictionary mapping unique sequence data to list of sequence IDs.
        """
        haplotypes = {}

        for seq in self._sequences:
            # Use sequence without gaps as haplotype key
            haplotype_key = seq.remove_gaps().data

            if haplotype_key not in haplotypes:
                haplotypes[haplotype_key] = []

            haplotypes[haplotype_key].append(seq.id)

        return haplotypes

    def to_fasta(self) -> str:
        """
        Convert alignment to FASTA format string.

        Returns
        -------
        str
            FASTA formatted string.
        """
        return '\n'.join(str(seq) for seq in self._sequences)

    def __str__(self) -> str:
        """Return string representation of alignment."""
        stats = self.calculate_stats()
        return (
            f'Alignment: {stats.num_sequences} sequences, '
            f'{stats.length} positions, '
            f'{stats.variable_sites} variable sites'
        )

    def __repr__(self) -> str:
        """Detailed representation of alignment."""
        return f'Alignment({len(self._sequences)} sequences, length={self.length})'
