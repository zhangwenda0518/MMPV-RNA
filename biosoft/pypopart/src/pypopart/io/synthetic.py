"""
Synthetic alignments with known structure, for tests and demos.

A real dataset rarely exercises every part of the pipeline at once. The
alignment built here does: nested clades, shared and repeated
haplotypes, a reticulation that forces a median vector, long branches
that push an edge past the tick-mark threshold, and a few gaps and
ambiguity codes. Everything is seeded, so the same call always returns
the same sequences.
"""

import random
from typing import Dict, List, Tuple

from ..core.alignment import Alignment
from ..core.sequence import Sequence

#: Populations in the structured dataset, in output order.
POPULATIONS = ('Coastal', 'Highland', 'Riverine', 'Desert', 'Island', 'Alpine')

#: Samples drawn per population.
SAMPLES_PER_POPULATION = 20

#: Mutations on the two long branches: past the default tick threshold
#: (10) but under the hard tick ceiling (30), so those edges draw a
#: numeral while every other edge draws a comb.
LONG_BRANCH_STEPS = 13


class _Mutator:
    """
    Apply fresh point mutations to a template sequence.

    Every mutation lands on a site not yet used anywhere in the dataset,
    so two haplotypes differ by exactly the mutations that separate them
    on the intended tree. The one reticulation is built explicitly from
    reused mutation sets, not by chance.

    Parameters
    ----------
    rng : random.Random
        Seeded source of sites and bases.
    length : int
        Alignment length.
    """

    def __init__(self, rng: random.Random, length: int):
        self.rng = rng
        self.length = length
        self.used: set = set()

    def mutation(self) -> Tuple[int, str]:
        """
        Draw a new (site, base) pair on an unused site.

        Returns
        -------
        Tuple[int, str]
            Zero-based site and the base to put there.
        """
        while True:
            site = self.rng.randrange(self.length)
            if site not in self.used:
                self.used.add(site)
                return site, self.rng.choice('ACGT')

    def apply(self, template: List[str], steps: int) -> List[str]:
        """
        Copy a template and add ``steps`` new mutations to it.

        Parameters
        ----------
        template : List[str]
            Sequence to mutate, as a list of bases.
        steps : int
            Number of mutations to add.

        Returns
        -------
        List[str]
            The mutated copy.
        """
        seq = list(template)
        for _ in range(steps):
            site, base = self.mutation()
            # Guarantee a change even when the drawn base matches.
            if seq[site] == base:
                base = 'ACGT'[('ACGT'.index(base) + 1) % 4]
            seq[site] = base
        return seq


def simulate_structured_alignment(
    seed: int = 7, length: int = 300
) -> Tuple[Alignment, Dict[str, str]]:
    """
    Build a six-population alignment with known network structure.

    Each population descends from its own founder a few steps off a
    shared ancestor; within a population some samples repeat the founder
    (node frequency), some sit on a sub-founder one or two steps on
    (chains, not just stars), and some carry private mutations. On top
    of that:

    - ``Island`` shares the ``Coastal`` founder haplotype, so that node
      is drawn as a two-population pie.
    - ``Desert`` holds three haplotypes two steps apart pairwise whose
      common neighbour was never sampled, which MJN reconstructs as a
      median vector.
    - ``Alpine`` sits :data:`LONG_BRANCH_STEPS` mutations from the
      ancestor and one ``Island`` sample is as far again from its
      founder, giving two edges labelled with a numeral.
    - Three ``Riverine`` samples carry a two-base gap and one
      ``Highland`` sample an ``N``.

    Parameters
    ----------
    seed : int, default=7
        RNG seed; the same seed always gives the same dataset.
    length : int, default=300
        Alignment length in bases. Must leave room for roughly seventy
        distinct mutation sites.

    Returns
    -------
    Tuple[Alignment, Dict[str, str]]
        The alignment and a ``{sample_id: population}`` mapping.
    """
    if length < 120:
        raise ValueError('length must be at least 120 to fit the mutations')

    rng = random.Random(seed)
    mutate = _Mutator(rng, length)
    ancestor = [rng.choice('ACGT') for _ in range(length)]

    sequences: List[Sequence] = []
    populations: Dict[str, str] = {}

    def add(pop: str, seq: List[str]) -> None:
        """
        Record one sample of a population.

        Parameters
        ----------
        pop : str
            Population name.
        seq : List[str]
            Sequence data.
        """
        index = sum(1 for p in populations.values() if p == pop) + 1
        sample_id = f'{pop[:3].upper()}_{index:02d}'
        populations[sample_id] = pop
        sequences.append(Sequence(sample_id, ''.join(seq)))

    founders: Dict[str, List[str]] = {}
    for pop in POPULATIONS:
        steps = LONG_BRANCH_STEPS if pop == 'Alpine' else rng.randrange(2, 5)
        founders[pop] = mutate.apply(ancestor, steps)

    for pop in POPULATIONS:
        founder = founders[pop]
        sub_founder = mutate.apply(founder, rng.randrange(1, 3))

        # Six on the founder, four on the sub-founder: repeated
        # haplotypes give the nodes size.
        for _ in range(6):
            add(pop, founder)
        for _ in range(4):
            add(pop, sub_founder)

        # Private variation: five singletons off each, one pair of
        # which is duplicated so a tip node also has frequency two.
        for base in (founder, sub_founder):
            for i in range(5):
                private = mutate.apply(base, rng.randrange(1, 3))
                add(pop, private)
                if i == 0:
                    add(pop, private)

    # Island shares the Coastal founder: replace three Island privates.
    island_ids = [s for s, p in populations.items() if p == 'Island']
    for sample_id in island_ids[-3:]:
        _replace(sequences, sample_id, founders['Coastal'])

    # One Island sample a long way from its founder.
    _replace(
        sequences,
        island_ids[-4],
        mutate.apply(founders['Island'], LONG_BRANCH_STEPS),
    )

    # Desert reticulation: X = f+a+b, Y = f+a+c, Z = f+b+c, with the
    # median f+a+b+c unsampled.
    desert = founders['Desert']
    a, b, c = (mutate.mutation() for _ in range(3))
    triad = []
    for pair in ((a, b), (a, c), (b, c)):
        seq = list(desert)
        for site, base in pair:
            seq[site] = (
                base if seq[site] != base else 'ACGT'[('ACGT'.index(base) + 1) % 4]
            )
        triad.append(seq)
    desert_ids = [s for s, p in populations.items() if p == 'Desert']
    for sample_id, seq in zip(desert_ids[-3:], triad):
        _replace(sequences, sample_id, seq)

    # Gaps and an ambiguity code.
    riverine_ids = [s for s, p in populations.items() if p == 'Riverine']
    gap_site = length // 2
    for sample_id in riverine_ids[-3:]:
        seq = list(_data(sequences, sample_id))
        seq[gap_site : gap_site + 2] = ['-', '-']
        _replace(sequences, sample_id, seq)
    highland_ids = [s for s, p in populations.items() if p == 'Highland']
    seq = list(_data(sequences, highland_ids[-1]))
    seq[gap_site + 10] = 'N'
    _replace(sequences, highland_ids[-1], seq)

    return Alignment(sequences), populations


def _data(sequences: List[Sequence], sample_id: str) -> str:
    """
    Look up a sequence's data by sample id.

    Parameters
    ----------
    sequences : List[Sequence]
        Sequences built so far.
    sample_id : str
        Id to find.

    Returns
    -------
    str
        The sequence string.
    """
    return next(s.data for s in sequences if s.id == sample_id)


def _replace(sequences: List[Sequence], sample_id: str, seq: List[str]) -> None:
    """
    Swap the data of one sequence in place.

    Parameters
    ----------
    sequences : List[Sequence]
        Sequences built so far.
    sample_id : str
        Id of the sequence to replace.
    seq : List[str]
        New sequence data.
    """
    for i, existing in enumerate(sequences):
        if existing.id == sample_id:
            sequences[i] = Sequence(sample_id, ''.join(seq))
            return
    raise KeyError(sample_id)
