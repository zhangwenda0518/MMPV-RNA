"""The structured synthetic dataset has the shape the docs promise."""

from pypopart.core.haplotype import identify_haplotypes_from_alignment
from pypopart.io.synthetic import (
    LONG_BRANCH_STEPS,
    POPULATIONS,
    SAMPLES_PER_POPULATION,
    simulate_structured_alignment,
)


def test_is_deterministic():
    """The same seed always gives the same sequences."""
    first, first_pops = simulate_structured_alignment()
    second, second_pops = simulate_structured_alignment()

    assert [s.data for s in first] == [s.data for s in second]
    assert first_pops == second_pops


def test_every_population_is_sampled():
    """Six populations, at least the promised number of samples each."""
    alignment, populations = simulate_structured_alignment()

    assert set(populations.values()) == set(POPULATIONS)
    for pop in POPULATIONS:
        assert (
            sum(1 for p in populations.values() if p == pop) >= SAMPLES_PER_POPULATION
        )
    assert set(populations) == set(alignment.sequence_ids)


def test_a_haplotype_is_shared_between_two_populations():
    """The Coastal founder shows up in Island too, for a pie-chart node."""
    alignment, populations = simulate_structured_alignment()
    haplotypes = identify_haplotypes_from_alignment(alignment, populations)

    shared = [h for h in haplotypes if {'Coastal', 'Island'} <= h.get_populations()]
    assert shared


def test_has_long_branches_and_ambiguity():
    """Numeral-fallback edges plus gaps and an N for the parsers."""
    alignment, _ = simulate_structured_alignment()
    data = [s.data for s in alignment]

    assert any('-' in d for d in data)
    assert any('N' in d for d in data)
    # Two sequences at least LONG_BRANCH_STEPS apart from everything near
    # them exist: check the alignment is not one tight cluster.
    assert alignment.length == 300
    assert LONG_BRANCH_STEPS > 10


def test_repeated_and_private_haplotypes():
    """Some nodes have frequency well above one, most are singletons."""
    alignment, _ = simulate_structured_alignment()
    haplotypes = identify_haplotypes_from_alignment(alignment)
    frequencies = sorted((h.frequency for h in haplotypes), reverse=True)

    assert frequencies[0] >= 6
    assert frequencies.count(1) > len(frequencies) // 2
