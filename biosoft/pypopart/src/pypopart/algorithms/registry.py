"""
Algorithm registry for PyPopART.

Single source of truth for available network construction algorithms,
used by the CLI and GUI instead of duplicated if/elif dispatch chains.
"""

from typing import Any, Dict, List, Type

from .base import NetworkAlgorithm
from .mjn import MedianJoiningNetwork
from .msn import MinimumSpanningNetwork
from .mst import MinimumSpanningTree
from .parsimony_net import ParsimonyNetwork
from .tcs import TCS
from .tsw import TightSpanWalker

#: Canonical algorithm name -> implementing class.
ALGORITHMS: Dict[str, Type[NetworkAlgorithm]] = {
    'mst': MinimumSpanningTree,
    'msn': MinimumSpanningNetwork,
    'tcs': TCS,
    'mjn': MedianJoiningNetwork,
    'pn': ParsimonyNetwork,
    'tsw': TightSpanWalker,
}

#: One-line description per algorithm, for CLI listings and GUI dropdowns.
DESCRIPTIONS: Dict[str, str] = {
    'mst': 'Minimum Spanning Tree',
    'msn': 'Minimum Spanning Network',
    'tcs': 'Statistical Parsimony (TCS)',
    'mjn': 'Median-Joining Network',
    'pn': 'Parsimony Network (consensus from sampled trees)',
    'tsw': 'Tight Span Walker',
}


def build(
    name: str, distance_method: str = 'hamming', **params: Any
) -> NetworkAlgorithm:
    """
    Construct a network algorithm instance by name.

    Parameters
    ----------
    name : str
        Algorithm name (case-insensitive): mst, msn, tcs, mjn, pn, tsw.
    distance_method : str, default='hamming'
        Distance metric name or alias.
    **params : dict
        Algorithm-specific parameters (e.g. epsilon for msn/mjn,
        confidence for tcs, random_seed for pn). Unknown parameter
        names raise TypeError from the algorithm constructor.

    Returns
    -------
    NetworkAlgorithm
        Configured algorithm instance.

    Raises
    ------
    ValueError
        If the algorithm name is not registered.
    """
    key = name.lower()
    if key not in ALGORITHMS:
        raise ValueError(
            f'Unknown algorithm: {name}. Available: {", ".join(sorted(ALGORITHMS))}'
        )
    return ALGORITHMS[key](distance_method=distance_method, **params)


def list_algorithms() -> List[Dict[str, str]]:
    """
    List registered algorithms with their descriptions.

    Returns
    -------
    list of dict
        One entry per algorithm with 'name' and 'description' keys,
        in canonical listing order.
    """
    return [{'name': name, 'description': DESCRIPTIONS[name]} for name in ALGORITHMS]
