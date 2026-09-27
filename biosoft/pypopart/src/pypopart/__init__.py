"""
PyPopART - Pure Python implementation of PopART haplotype network analysis.

The top-level package exposes the common workflow directly::

    from pypopart import load_alignment, build_network

    alignment = load_alignment('sequences.fasta')
    network = build_network('mjn', alignment)

Subpackages hold the full API: :mod:`pypopart.core` (domain model),
:mod:`pypopart.algorithms` (network construction), :mod:`pypopart.io`
(file formats), :mod:`pypopart.stats`, :mod:`pypopart.layout`, and
:mod:`pypopart.visualization`.
"""

try:
    from ._version import __version__
except ImportError:
    __version__ = '0.0.0.dev0'

from .algorithms import ALGORITHMS, build, list_algorithms
from .core.alignment import Alignment
from .core.distance import DistanceMatrix, pairwise_distance_matrix
from .core.graph import HaplotypeNetwork
from .core.haplotype import Haplotype, identify_haplotypes_from_alignment
from .core.sequence import Sequence
from .io import load_alignment, load_network, save_alignment, save_network


def build_network(algorithm: str, alignment, **params) -> HaplotypeNetwork:
    """
    Construct a haplotype network in one call.

    Parameters
    ----------
    algorithm : str
        Algorithm name: mst, msn, tcs, mjn, pn, or tsw.
    alignment : Alignment
        Multiple sequence alignment.
    **params : dict
        Algorithm-specific parameters (see pypopart.algorithms.build).

    Returns
    -------
    HaplotypeNetwork
        The constructed network.
    """
    return build(algorithm, **params).build_network(alignment)


__all__ = [
    '__version__',
    'ALGORITHMS',
    'Alignment',
    'DistanceMatrix',
    'Haplotype',
    'HaplotypeNetwork',
    'Sequence',
    'build',
    'build_network',
    'identify_haplotypes_from_alignment',
    'list_algorithms',
    'load_alignment',
    'load_network',
    'pairwise_distance_matrix',
    'save_alignment',
    'save_network',
]
