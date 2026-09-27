"""
Site-pattern condensation, ported from PopART's HapNet::condenseSitePats.

PopART condenses alignment columns before computing distances: columns
that are entirely ambiguous are dropped, and columns equivalent up to a
bijective character relabelling are collapsed into one representative
column whose *weight* records how many original columns it stands for.
Distances then count a weighted mismatch per condensed column, which is
what makes PopART's distances differ from naive per-site Hamming counts.
"""

from typing import List, Tuple
from typing import Sequence as SequenceType

#: Characters PopART treats as ambiguous for DNA (Sequence::isAmbiguousChar),
#: extended with '?' (binary/missing) and 'X' (amino acid) which PyPopART's
#: readers can emit for DNA data.
AMBIGUOUS_CHARS = frozenset('-NYRMSVWKDHB?X')


def is_ambiguous(char: str) -> bool:
    """
    Check whether a character is ambiguous (gap, N, or IUPAC code).

    Parameters
    ----------
    char : str
        Single sequence character (upper case).

    Returns
    -------
    bool
        True when the character is a gap or ambiguity code.
    """
    return char in AMBIGUOUS_CHARS


def condense_site_patterns(
    sequences: SequenceType[str],
) -> Tuple[List[str], List[int]]:
    """
    Collapse equivalent alignment columns into weighted site patterns.

    Port of ``HapNet::condenseSitePats``. Two columns are equivalent when
    a bijective character mapping converts one into the other row-by-row,
    with ambiguous characters required to match exactly. Columns that are
    ambiguous in every sequence are dropped entirely.

    Parameters
    ----------
    sequences : sequence of str
        Equal-length (aligned) sequence strings.

    Returns
    -------
    tuple of (list of str, list of int)
        The condensed sequence strings and, per condensed column, the
        number of original columns it represents (the site weight).

    Raises
    ------
    ValueError
        If sequences differ in length.
    """
    if not sequences:
        return [], []

    n_sites = len(sequences[0])
    if any(len(seq) != n_sites for seq in sequences):
        raise ValueError('Sequences must have the same length')

    dropped = n_sites  # sentinel meaning "column is dropped"
    same_pos_as = list(range(n_sites))

    for i in range(n_sites):
        column_i = [seq[i] for seq in sequences]

        # A column that is ambiguous in every sequence carries no signal
        if all(is_ambiguous(c) for c in column_i):
            same_pos_as[i] = dropped

        # Link the first later column equivalent to this one; chains of
        # links propagate each equivalence class's root (or the dropped
        # sentinel) across all members, as in the C++ implementation.
        for j in range(i + 1, n_sites):
            i2j: dict = {}
            j2i: dict = {}
            same = True

            for seq in sequences:
                char_i = seq[i]
                char_j = seq[j]

                if (is_ambiguous(char_i) or is_ambiguous(char_j)) and (
                    char_i != char_j
                ):
                    same = False
                    break
                if char_i not in i2j:
                    i2j[char_i] = char_j
                    if char_j not in j2i:
                        j2i[char_j] = char_i
                    elif j2i[char_j] != char_i:
                        same = False
                        break
                elif i2j[char_i] != char_j:
                    same = False
                    break

            if same:
                same_pos_as[j] = same_pos_as[i]
                break

    # Emit one representative per class, counting weights
    original_to_condensed = [dropped] * n_sites
    kept_columns: List[int] = []
    for i in range(n_sites):
        if same_pos_as[i] == dropped:
            continue
        if same_pos_as[i] < i:
            original_to_condensed[i] = original_to_condensed[same_pos_as[i]]
        else:
            original_to_condensed[i] = len(kept_columns)
            kept_columns.append(i)

    weights = [0] * len(kept_columns)
    for i in range(n_sites):
        if original_to_condensed[i] != dropped:
            weights[original_to_condensed[i]] += 1

    condensed = [''.join(seq[i] for i in kept_columns) for seq in sequences]
    return condensed, weights
