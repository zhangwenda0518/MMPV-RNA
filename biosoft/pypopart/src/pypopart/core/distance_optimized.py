"""
Optimized distance calculations using Numba JIT compilation.

This module provides high-performance distance calculations for large datasets
using Numba's just-in-time compilation to native machine code.
"""

import numba
import numpy as np

#: str.translate table folding non-gap IUPAC ambiguity codes to 'N', so
#: the byte kernels (which skip N/?) honour the full ambiguity set.
_AMBIGUITY_FOLD = str.maketrans(dict.fromkeys('YRMSVWKDHBX', 'N'))


@numba.jit(nopython=True, cache=True)
def hamming_distance_numba(
    seq1_bytes: np.ndarray, seq2_bytes: np.ndarray, ignore_gaps: bool = True
) -> int:
    """
    Calculate Hamming distance between two sequences using Numba JIT.

    Optimized for performance using compiled native code. Significantly faster
    than pure Python implementation for large sequences or many comparisons.

    Parameters
    ----------
    seq1_bytes : np.ndarray
        First sequence as numpy array of bytes.
    seq2_bytes : np.ndarray
        Second sequence as numpy array of bytes.
    ignore_gaps : bool, default=True
        Whether to ignore gap characters ('-').

    Returns
    -------
    int
        Int        Number of differing positions.

    Notes
    -----
    This function is JIT-compiled and cached for maximum performance.
    First call may be slower due to compilation overhead.

    N and ? characters are treated as ambiguous and do not count as
    mutations when compared to any base or to each other.
    """
    if len(seq1_bytes) != len(seq2_bytes):
        return -1  # Error indicator

    differences = 0
    gap_byte = ord('-')
    N_byte = ord('N')
    question_byte = ord('?')

    for i in range(len(seq1_bytes)):
        c1 = seq1_bytes[i]
        c2 = seq2_bytes[i]

        if ignore_gaps and (c1 == gap_byte or c2 == gap_byte):
            continue

        # Skip positions with N or ? (ambiguous bases)
        if c1 == N_byte or c1 == question_byte or c2 == N_byte or c2 == question_byte:
            continue

        if c1 != c2:
            differences += 1

    return differences


@numba.jit(nopython=True, cache=True, parallel=True)
def pairwise_hamming_matrix_numba(
    sequences: np.ndarray, ignore_gaps: bool = True
) -> np.ndarray:
    """
    Calculate pairwise Hamming distance matrix for multiple sequences.

    Uses parallel computation for improved performance on multi-core systems.

    Parameters
    ----------
    sequences : np.ndarray
        2D array where each row is a sequence (as bytes).
    ignore_gaps : bool, default=True
        Whether to ignore gap characters.

    Returns
    -------
    np.ndarray
        Np.ndarray        Symmetric distance matrix of shape (n_sequences, n_sequences).

    Notes
    -----
    This function uses Numba's parallel execution for significant speedup
    on multi-core CPUs. The parallel pragma distributes the outer loop
    iterations across available cores.
    """
    n_seqs = sequences.shape[0]
    matrix = np.zeros((n_seqs, n_seqs), dtype=np.int32)

    for i in numba.prange(n_seqs):
        for j in range(i + 1, n_seqs):
            dist = hamming_distance_numba(sequences[i], sequences[j], ignore_gaps)
            matrix[i, j] = dist
            matrix[j, i] = dist

    return matrix


def hamming_distance_optimized(seq1, seq2, ignore_gaps: bool = True) -> int:
    """
    Calculate Hamming distance with automatic Numba optimization.

    This is a wrapper around the Numba-optimized function that handles
    Sequence objects and string conversion automatically.

    Parameters
    ----------
    seq1 : Sequence or str
        First sequence.
    seq2 : Sequence or str
        Second sequence.
    ignore_gaps : bool, default=True
        Whether to ignore gap characters.

    Returns
    -------
    int
        Int        Hamming distance.

    Notes
    -----
    For single comparisons, the overhead of conversion may outweigh
    the performance benefit. Use this primarily for batch operations
    or repeated calls where JIT compilation is amortized.
    """
    # Convert to strings if Sequence objects
    s1 = seq1.data if hasattr(seq1, 'data') else str(seq1)
    s2 = seq2.data if hasattr(seq2, 'data') else str(seq2)

    # Fold IUPAC ambiguity codes to 'N' so the kernel skips them,
    # then convert to NumPy byte arrays
    s1 = s1.translate(_AMBIGUITY_FOLD)
    s2 = s2.translate(_AMBIGUITY_FOLD)
    seq1_bytes = np.frombuffer(s1.encode('ascii'), dtype=np.uint8)
    seq2_bytes = np.frombuffer(s2.encode('ascii'), dtype=np.uint8)

    return hamming_distance_numba(seq1_bytes, seq2_bytes, ignore_gaps)
