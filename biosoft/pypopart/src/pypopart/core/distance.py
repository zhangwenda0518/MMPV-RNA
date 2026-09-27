"""
Distance calculation and metrics for PyPopART.

Implements various distance metrics for DNA sequences including
Hamming distance and evolutionary models.
"""

from math import log
from typing import TYPE_CHECKING, Callable, List, Optional, Tuple

import numpy as np

from .alignment import Alignment
from .sequence import Sequence
from .site_patterns import AMBIGUOUS_CHARS

if TYPE_CHECKING:
    import matplotlib.figure

# Try to import optimized Numba versions
try:
    from .distance_optimized import hamming_distance_optimized

    _NUMBA_AVAILABLE = True
except ImportError:
    _NUMBA_AVAILABLE = False

#: Ambiguity codes always skipped in distance calculations. The gap
#: character is handled separately via ignore_gaps.
_AMBIGUOUS_NO_GAP = frozenset(AMBIGUOUS_CHARS - {'-'})

#: Byte translation table folding every non-gap ambiguity code to 'N',
#: so the byte-level kernels (which skip N/?) honour the full IUPAC set.
_FOLD_AMBIGUOUS = np.arange(256, dtype=np.uint8)
for _char in _AMBIGUOUS_NO_GAP:
    _FOLD_AMBIGUOUS[ord(_char)] = ord('N')


def _encode_for_kernel(strings: List[str]) -> np.ndarray:
    """
    Encode equal-length sequence strings for the byte-level kernels.

    Parameters
    ----------
    strings : list of str
        Equal-length sequence strings.

    Returns
    -------
    np.ndarray
        (n, L) uint8 array with ambiguity codes folded to 'N'.
    """
    n, length = len(strings), len(strings[0])
    encoded = np.zeros((n, length), dtype=np.uint8)
    for i, s in enumerate(strings):
        encoded[i] = np.frombuffer(s.encode('ascii'), dtype=np.uint8)
    return _FOLD_AMBIGUOUS[encoded]


def hamming_distance(
    seq1: Sequence,
    seq2: Sequence,
    ignore_gaps: bool = True,
    use_numba: bool = True,
) -> int:
    """
    Calculate Hamming distance between two sequences.

    Parameters
    ----------
    seq1 : Sequence
        First sequence.
    seq2 : Sequence
        Second sequence.
    ignore_gaps : bool, default=True
        Whether to ignore gap characters ('-').
    use_numba : bool, default=True
        Use Numba-optimized version if available.

    Returns
    -------
    int
        Int        Number of differing positions.

    Raises
    ------
    ValueError
        If sequences have different lengths

    Notes
    -----
    When use_numba=True and Numba is available, uses JIT-compiled
    optimized version for better performance on large datasets.

    N and ? characters are treated as ambiguous and do not count as
    mutations when compared to any base (A, T, G, C) or to each other.
    """
    if len(seq1) != len(seq2):
        raise ValueError(f'Sequences must have same length: {len(seq1)} vs {len(seq2)}')

    # Use Numba optimized version if available and requested
    if use_numba and _NUMBA_AVAILABLE:
        return hamming_distance_optimized(seq1, seq2, ignore_gaps)

    # Fall back to pure Python implementation
    differences = 0
    for c1, c2 in zip(seq1.data, seq2.data):
        if ignore_gaps and (c1 == '-' or c2 == '-'):
            continue
        # Skip ambiguous positions (N, ?, and IUPAC codes, as in PopART)
        if c1 in _AMBIGUOUS_NO_GAP or c2 in _AMBIGUOUS_NO_GAP:
            continue
        if c1 != c2:
            differences += 1

    return differences


def p_distance(seq1: Sequence, seq2: Sequence, ignore_gaps: bool = True) -> float:
    """
    Calculate p-distance (proportion of differing sites).

    N and ? characters are treated as ambiguous and do not count as
    mutations when compared to any base (A, T, G, C) or to each other.

    Parameters
    ----------
    seq1 : Sequence
        First sequence.
    seq2 : Sequence
        Second sequence.
    ignore_gaps : bool, default=True
        Whether to skip gap positions.

    Returns
    -------
    float
        The p-distance (proportion of differing sites).
    """
    if len(seq1) != len(seq2):
        raise ValueError(f'Sequences must have same length: {len(seq1)} vs {len(seq2)}')

    differences = 0
    compared_sites = 0

    for c1, c2 in zip(seq1.data, seq2.data):
        if ignore_gaps and (c1 == '-' or c2 == '-'):
            continue
        # Skip positions with N or ? (ambiguous bases)
        if c1 in 'N?' or c2 in 'N?':
            continue
        compared_sites += 1
        if c1 != c2:
            differences += 1

    if compared_sites == 0:
        raise ValueError('No valid sites to compare')

    return differences / compared_sites


def jukes_cantor_distance(
    seq1: Sequence, seq2: Sequence, ignore_gaps: bool = True
) -> float:
    """
    Calculate Jukes-Cantor corrected distance.

    Parameters
    ----------
    seq1 : Sequence
        First sequence.
    seq2 : Sequence
        Second sequence.
    ignore_gaps : bool, default=True
        Whether to skip gap positions.

    Returns
    -------
    float
        The Jukes-Cantor corrected distance.
    """
    p = p_distance(seq1, seq2, ignore_gaps)

    if p >= 0.75:
        raise ValueError(
            f'Sequences too divergent for Jukes-Cantor correction (p={p:.3f} >= 0.75)'
        )

    distance = -0.75 * log(1 - (4.0 / 3.0) * p)
    return distance


def kimura_2p_distance(
    seq1: Sequence, seq2: Sequence, ignore_gaps: bool = True
) -> float:
    """
    Calculate Kimura 2-parameter distance.

    Parameters
    ----------
    seq1 : Sequence
        First sequence.
    seq2 : Sequence
        Second sequence.
    ignore_gaps : bool, default=True
        Whether to skip gap positions.

    Returns
    -------
    float
        The Kimura 2-parameter distance.
    """
    if len(seq1) != len(seq2):
        raise ValueError(f'Sequences must have same length: {len(seq1)} vs {len(seq2)}')

    transitions = 0
    transversions = 0
    compared_sites = 0

    transition_pairs = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}

    for c1, c2 in zip(seq1.data, seq2.data):
        if ignore_gaps and (c1 == '-' or c2 == '-'):
            continue
        if c1 in 'N?' or c2 in 'N?':
            continue

        compared_sites += 1

        if c1 != c2:
            if (c1, c2) in transition_pairs:
                transitions += 1
            else:
                transversions += 1

    if compared_sites == 0:
        raise ValueError('No valid sites to compare')

    P = transitions / compared_sites
    Q = transversions / compared_sites

    term1 = 1 - 2 * P - Q
    term2 = 1 - 2 * Q

    if term1 <= 0 or term2 <= 0:
        raise ValueError(
            f'Sequences too divergent for K2P correction (P={P:.3f}, Q={Q:.3f})'
        )

    distance = -0.5 * log(term1 * (term2**0.5))
    return distance


def tamura_nei_distance(
    seq1: Sequence, seq2: Sequence, ignore_gaps: bool = True
) -> float:
    """
    Calculate Tamura-Nei distance.

    The Tamura-Nei model accounts for:
    - Different base frequencies (GC content)
    - Different rates for transitions within purines (A<->G) and pyrimidines (C<->T)
    - Different rate for transversions

    Parameters
    ----------
    seq1 : Sequence
        First sequence.
    seq2 : Sequence
        Second sequence.
    ignore_gaps : bool, default=True
        Whether to ignore gap positions.

    Returns
    -------
    float
        Tamura-Nei corrected distance.

    Raises
    ------
        ValueError: If sequences have different lengths or are too divergent

    Reference:
        Tamura K, Nei M (1993) Mol Biol Evol 10(3):512-526
    """
    if len(seq1) != len(seq2):
        raise ValueError(f'Sequences must have same length: {len(seq1)} vs {len(seq2)}')

    # Count base frequencies and differences
    purine_transitions = 0  # A<->G
    pyrimidine_transitions = 0  # C<->T
    transversions = 0
    compared_sites = 0

    base_counts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}

    for c1, c2 in zip(seq1.data, seq2.data):
        if ignore_gaps and (c1 == '-' or c2 == '-'):
            continue
        if c1 in 'N?' or c2 in 'N?':
            continue

        compared_sites += 1

        # Count base frequencies
        if c1 in base_counts:
            base_counts[c1] += 1
        if c2 in base_counts:
            base_counts[c2] += 1

        # Count differences
        if c1 != c2:
            if (c1, c2) in {('A', 'G'), ('G', 'A')}:
                purine_transitions += 1
            elif (c1, c2) in {('C', 'T'), ('T', 'C')}:
                pyrimidine_transitions += 1
            else:
                transversions += 1

    if compared_sites == 0:
        raise ValueError('No valid sites to compare')

    # Calculate base frequencies
    total_bases = sum(base_counts.values())
    if total_bases == 0:
        raise ValueError('No valid bases found')

    freq_A = base_counts['A'] / total_bases
    freq_G = base_counts['G'] / total_bases
    freq_C = base_counts['C'] / total_bases
    freq_T = base_counts['T'] / total_bases

    # Calculate GC content components
    freq_purines = freq_A + freq_G
    freq_pyrimidines = freq_C + freq_T

    # Calculate proportions of differences
    P1 = purine_transitions / compared_sites  # A<->G transitions
    P2 = pyrimidine_transitions / compared_sites  # C<->T transitions
    Q = transversions / compared_sites  # Transversions

    # Calculate GC content and related parameters
    freq_G + freq_C

    # Handle edge cases
    if freq_purines == 0 or freq_pyrimidines == 0:
        # Fall back to Kimura 2-parameter
        P = P1 + P2
        term1 = 1 - 2 * P - Q
        term2 = 1 - 2 * Q
        if term1 <= 0 or term2 <= 0:
            raise ValueError(
                f'Sequences too divergent for TN correction (P={P:.3f}, Q={Q:.3f})'
            )
        return -0.5 * log(term1 * (term2**0.5))

    # Calculate correction terms
    h_R = (
        2 * freq_A * freq_G / freq_purines if freq_purines > 0 else 0
    )  # Purine heterozygosity
    h_Y = (
        2 * freq_C * freq_T / freq_pyrimidines if freq_pyrimidines > 0 else 0
    )  # Pyrimidine heterozygosity

    # Calculate w terms for distance calculation
    # Use small epsilon to avoid division by zero
    epsilon = 1e-10

    if h_R > epsilon:
        w1 = 1 - P1 / h_R - Q / (2 * freq_purines * freq_pyrimidines)
    else:
        # When h_R is very small, use a simplified form
        w1 = (
            1 - P1 - Q / (2 * freq_purines * freq_pyrimidines)
            if freq_purines * freq_pyrimidines > 0
            else 1
        )

    if h_Y > epsilon:
        w2 = 1 - P2 / h_Y - Q / (2 * freq_purines * freq_pyrimidines)
    else:
        # When h_Y is very small, use a simplified form
        w2 = (
            1 - P2 - Q / (2 * freq_purines * freq_pyrimidines)
            if freq_purines * freq_pyrimidines > 0
            else 1
        )

    w3 = (
        1 - Q / (2 * freq_pyrimidines * freq_purines)
        if freq_purines * freq_pyrimidines > 0
        else 1
    )

    # Check for divergence
    if w1 <= 0 or w2 <= 0 or w3 <= 0:
        raise ValueError(
            f'Sequences too divergent for Tamura-Nei correction '
            f'(P1={P1:.3f}, P2={P2:.3f}, Q={Q:.3f})'
        )

    # Calculate Tamura-Nei distance
    if h_R > epsilon and h_Y > epsilon:
        # Full Tamura-Nei formula
        distance = (
            -h_R * log(w1)
            - h_Y * log(w2)
            - (freq_purines * freq_pyrimidines - h_R * h_Y / (h_R + h_Y)) * log(w3)
        )
    else:
        # Simplified form when heterozygosity is low
        distance = -0.5 * log(w1 * w2 * w3)

    return distance


class DistanceMatrix:
    """
    Store and manage a pairwise distance matrix.

    Parameters
    ----------
    labels : List[str]
        Sequence labels, in matrix order.
    matrix : np.ndarray, optional
        Square matrix of pairwise distances. A zero matrix of the
        right shape is created when omitted.

    Raises
    ------
    ValueError
        If the matrix shape does not match the number of labels.
    """

    def __init__(self, labels: List[str], matrix: Optional[np.ndarray] = None):
        """
        Initialize distance matrix.

        Parameters
        ----------
        labels : List[str]
            Sequence labels, in matrix order.
        matrix : np.ndarray, optional
            Square matrix of pairwise distances.
        """
        self.labels = labels
        self.n = len(labels)
        self._label_index = {label: i for i, label in enumerate(labels)}

        if matrix is not None:
            if matrix.shape != (self.n, self.n):
                raise ValueError(
                    f"Matrix shape {matrix.shape} doesn't match labels {self.n}"
                )
            self.matrix = matrix
        else:
            self.matrix = np.zeros((self.n, self.n))

    def get_distance(self, label1: str, label2: str) -> float:
        """
        Get distance between two sequences by label.

        Parameters
        ----------
        label1 : str
            First sequence label.
        label2 : str
            Second sequence label.

        Returns
        -------
        float
            The distance between the two labelled sequences.
        """
        i = self._label_index[label1]
        j = self._label_index[label2]
        return self.matrix[i, j]

    def set_distance(self, label1: str, label2: str, distance: float) -> None:
        """
        Set distance between two sequences.

        Parameters
        ----------
        label1 : str
            First sequence label.
        label2 : str
            Second sequence label.
        distance : float
            Distance value.
        """
        i = self._label_index[label1]
        j = self._label_index[label2]
        self.matrix[i, j] = distance
        self.matrix[j, i] = distance

    def get_row(self, label: str) -> np.ndarray:
        """
        Get all distances for a sequence.

        Parameters
        ----------
        label : str
            Sequence label.

        Returns
        -------
        np.ndarray
            All distances for a sequence.
        """
        i = self._label_index[label]
        return self.matrix[i, :]

    def get_min_distance(self, exclude_zero: bool = True) -> float:
        """
        Get minimum distance in matrix.

        Parameters
        ----------
        exclude_zero : bool, default=True
            Whether to ignore zero distances.

        Returns
        -------
        float
            The smallest pairwise distance in the matrix.
        """
        if exclude_zero:
            mask = np.triu(np.ones_like(self.matrix, dtype=bool), k=1)
            return self.matrix[mask].min()
        else:
            return self.matrix.min()

    def get_max_distance(self) -> float:
        """
        Get maximum distance in matrix.

        Returns
        -------
        float
            The largest pairwise distance in the matrix.
        """
        return self.matrix.max()

    def to_dict(self) -> dict:
        """
        Convert matrix to dictionary representation.

        Returns
        -------
        dict
            Dictionary representation of the matrix.
        """
        return {'labels': self.labels, 'matrix': self.matrix.tolist()}

    @classmethod
    def from_dict(cls, data: dict) -> 'DistanceMatrix':
        """
        Create distance matrix from dictionary.

        Parameters
        ----------
        data : dict
            Dictionary with 'labels' and 'matrix' entries.

        Returns
        -------
        DistanceMatrix
            The reconstructed distance matrix.
        """
        return cls(labels=data['labels'], matrix=np.array(data['matrix']))

    def visualize(
        self,
        title: str = 'Distance Matrix',
        cmap: str = 'viridis',
        figsize: Optional[Tuple[int, int]] = None,
        show_values: bool = True,
        save_path: Optional[str] = None,
    ) -> 'matplotlib.figure.Figure':
        """
        Visualize distance matrix as a heatmap.

        Parameters
        ----------
        title : str, default='Distance Matrix'
            Plot title.
        cmap : str, default='viridis'
            Matplotlib colormap name.
        figsize : Tuple[int, int], optional
            Figure size (width, height), auto-calculated if None.
        show_values : bool, default=True
            Whether to show distance values in cells.
        save_path : str, optional
            Optional path to save figure.

        Returns
        -------
        matplotlib.figure.Figure
            Matplotlib figure object.
            Raises :
            ImportError :
            If matplotlib is not installed.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError(
                'Matplotlib is required for visualization. '
                'Install it with: pip install matplotlib'
            ) from None

        # Auto-calculate figure size based on matrix size
        if figsize is None:
            size = max(8, min(20, self.n * 0.5))
            figsize = (size, size)

        fig, ax = plt.subplots(figsize=figsize)

        # Create heatmap
        im = ax.imshow(self.matrix, cmap=cmap, aspect='auto')

        # Set ticks and labels
        ax.set_xticks(np.arange(self.n))
        ax.set_yticks(np.arange(self.n))
        ax.set_xticklabels(self.labels, rotation=45, ha='right')
        ax.set_yticklabels(self.labels)

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Distance', rotation=270, labelpad=15)

        # Show values in cells if requested and matrix is not too large
        if show_values and self.n <= 20:
            for i in range(self.n):
                for j in range(self.n):
                    value = self.matrix[i, j]
                    if np.isfinite(value):
                        text_color = (
                            'white' if value > self.matrix.max() * 0.5 else 'black'
                        )
                        ax.text(
                            j,
                            i,
                            f'{value:.2f}',
                            ha='center',
                            va='center',
                            color=text_color,
                            fontsize=8,
                        )

        ax.set_title(title)
        fig.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')

        return fig

    def to_csv(self, filepath: str, delimiter: str = ',') -> None:
        """
        Export distance matrix to CSV file.

        Parameters
        ----------
        filepath : str
            Path to output CSV file.
        delimiter : str, default=','
            Delimiter character (default: comma).
        """
        with open(filepath, 'w') as f:
            # Write header
            f.write(delimiter.join([''] + self.labels) + '\n')

            # Write data rows
            for i, label in enumerate(self.labels):
                row_data = [label] + [f'{val:.6f}' for val in self.matrix[i]]
                f.write(delimiter.join(row_data) + '\n')

    @classmethod
    def from_csv(cls, filepath: str, delimiter: str = ',') -> 'DistanceMatrix':
        """
        Import distance matrix from CSV file.

        Parameters
        ----------
        filepath : str
            Path to CSV file.
        delimiter : str, default=','
            Delimiter character (default: comma).

        Returns
        -------
        DistanceMatrix
            DistanceMatrix object.
        """
        with open(filepath, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]

        # Parse header (column labels)
        header = lines[0].split(delimiter)
        labels = header[1:]  # Skip first empty cell

        # Parse data
        n = len(labels)
        matrix = np.zeros((n, n))

        for i, line in enumerate(lines[1:]):
            parts = line.split(delimiter)
            parts[0]
            values = [float(v) for v in parts[1:]]
            matrix[i] = values

        return cls(labels, matrix)

    def __str__(self) -> str:
        """Return string representation."""
        return f'DistanceMatrix({self.n} sequences)'

    def __repr__(self) -> str:
        """Detailed representation."""
        return f'DistanceMatrix(n={self.n}, min={self.get_min_distance():.4f}, max={self.get_max_distance():.4f})'


def calculate_distance_matrix(
    alignment: Alignment,
    distance_func: Optional[Callable[[Sequence, Sequence], float]] = None,
    **kwargs,
) -> DistanceMatrix:
    """
    Calculate pairwise distance matrix for alignment.

    Parameters
    ----------
    alignment : Alignment
        Multiple sequence alignment.
    distance_func : Callable[[Sequence, Sequence], float], optional
        Pairwise distance function to apply.
    **kwargs : dict
        Additional arguments passed to the distance function.

    Returns
    -------
    DistanceMatrix
        The pairwise distance matrix for the alignment.
    """
    if distance_func is None:
        distance_func = hamming_distance

    n = len(alignment)
    labels = alignment.sequence_ids
    matrix = np.zeros((n, n))

    sequences = list(alignment)
    for i in range(n):
        for j in range(i + 1, n):
            try:
                dist = distance_func(sequences[i], sequences[j], **kwargs)
                matrix[i, j] = matrix[j, i] = dist
            except ValueError:
                matrix[i, j] = matrix[j, i] = np.inf

    return DistanceMatrix(labels, matrix)


#: Canonical distance-method names and their accepted aliases.
DISTANCE_METHOD_ALIASES = {
    'hamming': 'hamming',
    'p': 'p',
    'jc': 'jc',
    'jukes-cantor': 'jc',
    'jukes_cantor': 'jc',
    'k2p': 'k2p',
    'kimura': 'k2p',
    'kimura_2p': 'k2p',
    'tn': 'tn',
    'tamura-nei': 'tn',
    'tamura_nei': 'tn',
}


def normalize_distance_method(method: str) -> str:
    """
    Resolve a distance-method name or alias to its canonical form.

    Parameters
    ----------
    method : str
        Method name or alias (case-insensitive), e.g. 'tamura_nei' or 'tn'.

    Returns
    -------
    str
        Canonical method name: 'hamming', 'p', 'jc', 'k2p', or 'tn'.

    Raises
    ------
    ValueError
        If the method is not recognised.
    """
    canonical = DISTANCE_METHOD_ALIASES.get(method.lower())
    if canonical is None:
        valid = sorted(set(DISTANCE_METHOD_ALIASES))
        raise ValueError(f'Unknown distance method: {method}. Valid: {valid}')
    return canonical


def _distance_function(method: str) -> Callable:
    """
    Return the scalar distance function for a canonical method name.

    Parameters
    ----------
    method : str
        Canonical method name from normalize_distance_method.

    Returns
    -------
    Callable
        Function taking (seq1, seq2, ignore_gaps=...) and returning a float.
    """
    return {
        'hamming': hamming_distance,
        'p': p_distance,
        'jc': jukes_cantor_distance,
        'k2p': kimura_2p_distance,
        'tn': tamura_nei_distance,
    }[method]


def sequence_distance(
    seq1, seq2, method: str = 'hamming', ignore_gaps: bool = True
) -> float:
    """
    Calculate the distance between two sequences with a named method.

    Parameters
    ----------
    seq1 : Sequence or Haplotype
        First sequence.
    seq2 : Sequence or Haplotype
        Second sequence.
    method : str, default='hamming'
        Distance method name or alias.
    ignore_gaps : bool, default=True
        Whether to ignore gap positions.

    Returns
    -------
    float
        Distance between the two sequences.
    """
    func = _distance_function(normalize_distance_method(method))
    return float(func(seq1, seq2, ignore_gaps=ignore_gaps))


def pairwise_distance_matrix(
    sequences,
    method: str = 'hamming',
    ignore_gaps: bool = True,
    mask: Optional[np.ndarray] = None,
    site_weights: Optional[np.ndarray] = None,
) -> DistanceMatrix:
    """
    Calculate a pairwise distance matrix for sequences or haplotypes.

    This is the single shared distance path for all network algorithms.
    For Hamming distances it uses a whole-matrix numba kernel when numba
    is installed, or a vectorised numpy fallback otherwise; both match the
    scalar function's gap and IUPAC-ambiguity handling. Other methods fall
    back to a per-pair loop over the scalar distance functions.

    Parameters
    ----------
    sequences : Alignment or list of Sequence or list of Haplotype
        Sequences to compare. Items need `id` and `data` attributes.
    method : str, default='hamming'
        Distance method name or alias (see normalize_distance_method).
    ignore_gaps : bool, default=True
        Whether to ignore gap positions.
    mask : np.ndarray of bool, optional
        Per-column mask; columns where the mask is False are excluded
        from the calculation (PopART's character masking).
    site_weights : np.ndarray, optional
        Per-column weights, e.g. from site-pattern condensation; a
        mismatch at column i contributes site_weights[i]. Hamming only.

    Returns
    -------
    DistanceMatrix
        Symmetric matrix of pairwise distances, labelled by sequence id.

    Raises
    ------
    ValueError
        If site_weights is combined with a non-Hamming method, or the
        mask/weights lengths don't match the alignment length.
    """
    method = normalize_distance_method(method)
    items = list(sequences)
    labels = [item.id for item in items]
    n = len(items)
    strings = [item.data for item in items]

    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        if strings and len(mask) != len(strings[0]):
            raise ValueError('mask length must match sequence length')
        kept = np.flatnonzero(mask)
        strings = [''.join(s[i] for i in kept) for s in strings]
        if site_weights is not None:
            site_weights = np.asarray(site_weights)[kept]

    if site_weights is not None and method != 'hamming':
        raise ValueError('site_weights are only supported for hamming distances')

    equal_length = n > 1 and strings and all(len(s) == len(strings[0]) for s in strings)

    if method == 'hamming' and equal_length:
        length = len(strings[0])
        if site_weights is not None and len(site_weights) != length:
            raise ValueError('site_weights length must match sequence length')

        if site_weights is None and _NUMBA_AVAILABLE:
            from .distance_optimized import pairwise_hamming_matrix_numba

            matrix = pairwise_hamming_matrix_numba(
                _encode_for_kernel(strings), ignore_gaps
            )
            return DistanceMatrix(labels, matrix.astype(float))

        # Vectorised numpy path (also handles weighted distances)
        encoded = _encode_for_kernel(strings)
        invalid = (encoded == ord('N')) | (encoded == ord('?'))
        if ignore_gaps:
            invalid |= encoded == ord('-')
        valid = ~invalid

        weights = (
            np.ones(length) if site_weights is None else np.asarray(site_weights, float)
        )
        matrix = np.zeros((n, n))
        for i in range(n):
            comparable = valid[i] & valid
            diffs = (encoded[i] != encoded) & comparable
            matrix[i] = diffs @ weights
        return DistanceMatrix(labels, matrix)

    # Per-pair scalar loop for corrected distances (and tiny inputs)
    if mask is not None:
        # Rebuild lightweight records over masked strings
        items = [Sequence(id=label, data=s) for label, s in zip(labels, strings)]
    func = _distance_function(method)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            try:
                dist = func(items[i], items[j], ignore_gaps=ignore_gaps)
            except ValueError:
                dist = np.inf
            matrix[i, j] = matrix[j, i] = dist

    return DistanceMatrix(labels, matrix)


def calculate_pairwise_distances(
    alignment: Alignment, method: str = 'hamming', ignore_gaps: bool = True
) -> DistanceMatrix:
    """
    Calculate pairwise distances using specified method.

    Parameters
    ----------
    alignment : Alignment
        Multiple sequence alignment.
    method : str, default='hamming'
        Distance method name or alias.
    ignore_gaps : bool, default=True
        Whether to skip gap positions.

    Returns
    -------
    DistanceMatrix
        The pairwise distance matrix.
    """
    return pairwise_distance_matrix(alignment, method=method, ignore_gaps=ignore_gaps)


class DistanceCalculator:
    """
    Convenience class for calculating distance matrices.

    Provides a simple interface for distance calculation with
    different evolutionary models.

    Parameters
    ----------
    method : str, default='hamming'
        Distance method: 'hamming', 'jc', 'k2p', 'tamura_nei'.
    ignore_gaps : bool, default=True
        Whether to ignore gaps in calculations.
    """

    def __init__(self, method: str = 'hamming', ignore_gaps: bool = True):
        """
        Initialize distance calculator.

        Parameters
        ----------
        method : str, default='hamming'
            Distance method: 'hamming', 'jc', 'k2p', 'tamura_nei'.
        ignore_gaps : bool, default=True
            Whether to ignore gaps in calculations.
        """
        self.method = method.lower()
        self.ignore_gaps = ignore_gaps

        # Map method names to functions
        if self.method == 'hamming':
            self.distance_func = hamming_distance
        elif self.method in ('jc', 'jukes_cantor', 'jukes-cantor'):
            self.distance_func = jukes_cantor_distance
        elif self.method in ('k2p', 'kimura', 'kimura_2p'):
            self.distance_func = kimura_2p_distance
        elif self.method in ('tn', 'tamura_nei', 'tamura-nei'):
            self.distance_func = tamura_nei_distance
        elif self.method == 'p':
            self.distance_func = p_distance
        else:
            raise ValueError(f'Unknown distance method: {method}')

    def calculate(self, seq1: Sequence, seq2: Sequence) -> float:
        """
        Calculate distance between two sequences.

        Parameters
        ----------
        seq1 : Sequence
            First sequence.
        seq2 : Sequence
            Second sequence.

        Returns
        -------
        float
            Float            Distance value.
        """
        return self.distance_func(seq1, seq2, ignore_gaps=self.ignore_gaps)

    def calculate_matrix(self, alignment: Alignment) -> np.ndarray:
        """
        Calculate pairwise distance matrix for alignment.

        Parameters
        ----------
        alignment : Alignment
            Sequence alignment.

        Returns
        -------
        np.ndarray
            Np.ndarray            Square distance matrix.
        """
        dist_matrix = calculate_pairwise_distances(
            alignment, method=self.method, ignore_gaps=self.ignore_gaps
        )
        return dist_matrix.matrix

    def __repr__(self) -> str:
        """Return detailed string representation."""
        return f"DistanceCalculator(method='{self.method}', ignore_gaps={self.ignore_gaps})"
