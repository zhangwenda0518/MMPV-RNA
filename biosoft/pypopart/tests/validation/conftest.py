"""
Shared helpers for golden-network validation tests.

Golden files (goldens/*.json) are self-contained: they carry the input
sequences, the algorithm and parameters to run, and the expected network
as a canonical node set and edge list. Expected values are hand-traced
from the C++ PopART algorithms in ./popart-current (see each file's
"derivation" note); they are the ground truth the Python port must match.

Median/intermediate vertices get generated ids, so edges are compared
after canonicalisation: sampled nodes keep their haplotype ids, inferred
nodes are matched by degree and adjacent sampled nodes.
"""

import json
from pathlib import Path

import pytest

from pypopart.algorithms import build
from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence

GOLDEN_DIR = Path(__file__).parent / 'goldens'


def golden_files():
    """
    Enumerate golden JSON files for parametrisation.

    Returns
    -------
    list of Path
        All golden files, sorted by name.
    """
    return sorted(GOLDEN_DIR.glob('*.json'))


def load_golden(path):
    """
    Load one golden file.

    Parameters
    ----------
    path : Path
        Golden JSON file.

    Returns
    -------
    dict
        Parsed golden specification.
    """
    return json.loads(path.read_text())


def run_golden(golden):
    """
    Run the algorithm a golden file specifies on its input sequences.

    Parameters
    ----------
    golden : dict
        Parsed golden specification.

    Returns
    -------
    HaplotypeNetwork
        The constructed network.
    """
    alignment = Alignment(
        [Sequence(sid, data) for sid, data in golden['sequences'].items()]
    )
    algo = build(golden['algorithm'], **golden.get('params', {}))
    return algo.build_network(alignment)


def canonical_edges(network):
    """
    Canonicalise a network's edges for comparison.

    Parameters
    ----------
    network : HaplotypeNetwork
        Network to canonicalise.

    Returns
    -------
    set of tuple
        Edges as (sorted node pair, distance) tuples.
    """
    edges = set()
    for u, v, attrs in network.graph.edges(data=True):
        pair = tuple(sorted((str(u), str(v))))
        edges.add((pair[0], pair[1], float(attrs.get('distance', 0))))
    return edges


@pytest.fixture(params=golden_files(), ids=lambda p: p.stem)
def golden(request):
    """Parametrised fixture yielding each parsed golden file."""
    return load_golden(request.param)
