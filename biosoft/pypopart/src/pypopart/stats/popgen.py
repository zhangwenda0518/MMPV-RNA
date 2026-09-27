"""
Population genetics statistics for PyPopART.

This module provides functions for calculating population genetics
measures including Tajima's D, Fu's Fs, FST, and AMOVA.
"""

from collections import defaultdict
from dataclasses import dataclass
import math
from typing import Dict, List, Optional, Tuple

import numpy as np

from ..core.alignment import Alignment
from ..core.graph import HaplotypeNetwork


@dataclass
class TajimaDResult:
    """Result of Tajima's D test."""

    D: float
    pi: float
    theta_w: float
    n_segregating_sites: int
    n_samples: int


@dataclass
class FuFsResult:
    """Result of Fu's Fs test."""

    Fs: float
    theta_pi: float
    n_haplotypes: int
    n_samples: int


@dataclass
class FstResult:
    """Result of pairwise FST calculation."""

    fst: float
    pop1: str
    pop2: str
    hs: float  # Within-population diversity
    ht: float  # Total diversity


@dataclass
class AMOVAResult:
    """Result of AMOVA (Analysis of Molecular Variance)."""

    phi_st: float  # Among-population variance
    phi_sc: float  # Among-groups-within-populations variance (if groups defined)
    phi_ct: float  # Among-groups variance (if groups defined)
    variance_among_pops: float
    variance_within_pops: float
    percent_among_pops: float
    percent_within_pops: float


def calculate_tajimas_d(
    alignment: Alignment, populations: Optional[Dict[str, str]] = None
) -> TajimaDResult:
    """
    Calculate Tajima's D statistic.

    Tajima's D tests for deviation from neutral evolution by comparing
    two estimates of theta (population mutation rate):
    - π (nucleotide diversity)
    - θ_W (Watterson's estimator based on segregating sites)

    D = 0: neutral evolution
    D > 0: balancing selection or population contraction
    D < 0: directional selection or population expansion

    Parameters
    ----------
    alignment : Alignment
        Alignment object.
    populations : Dict[str, str], optional
        Optional dict mapping sequence_id -> population.

    Returns
    -------
    TajimaDResult
        TajimaDResult with D statistic and related values.
    """
    n = len(alignment)  # Number of sequences
    L = alignment.length  # Sequence length

    if n < 2:
        return TajimaDResult(
            D=0.0, pi=0.0, theta_w=0.0, n_segregating_sites=0, n_samples=n
        )

    # Calculate number of segregating sites (S)
    S = 0
    for pos in range(L):
        column = alignment.get_column(pos)
        # Remove gaps
        column_no_gaps = [c for c in column if c != '-']
        if len(set(column_no_gaps)) > 1:
            S += 1

    # Calculate Watterson's estimator (θ_W)
    # θ_W = S / a1
    a1 = sum(1.0 / i for i in range(1, n))
    theta_w = S / a1 if a1 > 0 else 0.0

    # Calculate nucleotide diversity (π)
    # π = average pairwise differences
    total_diff = 0
    comparisons = 0

    for i in range(n):
        for j in range(i + 1, n):
            seq_i = alignment[i].data
            seq_j = alignment[j].data

            diff = sum(
                1 for a, b in zip(seq_i, seq_j) if a != b and a != '-' and b != '-'
            )
            total_diff += diff
            comparisons += 1

    pi = total_diff / comparisons if comparisons > 0 else 0.0

    # Calculate Tajima's D
    # D = (π - θ_W) / sqrt(V(π - θ_W))

    # Calculate variance terms
    a2 = sum(1.0 / (i**2) for i in range(1, n))

    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n**2 + n + 3) / (9 * n * (n - 1))

    c1 = b1 - 1 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1**2)

    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)

    variance = e1 * S + e2 * S * (S - 1)

    if variance > 0:
        D = (pi - theta_w) / math.sqrt(variance)
    else:
        D = 0.0

    return TajimaDResult(
        D=D, pi=pi, theta_w=theta_w, n_segregating_sites=S, n_samples=n
    )


def calculate_fu_fs(network: HaplotypeNetwork, alignment: Alignment) -> FuFsResult:
    """
    Calculate Fu's Fs statistic.

    Fu's Fs is based on the probability of observing a random neutral sample
    with as many or more haplotypes as observed, given theta estimated from
    nucleotide diversity.

    Negative Fs: excess of recent mutations (population expansion)
    Positive Fs: deficiency of alleles (balancing selection or bottleneck)

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.
    alignment : Alignment
        Alignment object.

    Returns
    -------
    FuFsResult
        FuFsResult with Fs statistic and related values.
    """
    # Get number of haplotypes and samples
    k = len(network.haplotypes)  # Number of haplotypes
    n = network.get_total_samples()  # Total samples

    if n < 2 or k < 1:
        return FuFsResult(Fs=0.0, theta_pi=0.0, n_haplotypes=k, n_samples=n)

    # Calculate nucleotide diversity (π)
    # Use average pairwise differences
    total_diff = 0
    comparisons = 0

    haplotype_list = network.haplotypes
    for i, hap_i in enumerate(haplotype_list):
        for _j, hap_j in enumerate(haplotype_list[i + 1 :], start=i + 1):
            seq_i = hap_i.sequence.data
            seq_j = hap_j.sequence.data

            diff = sum(
                1 for a, b in zip(seq_i, seq_j) if a != b and a != '-' and b != '-'
            )

            # Weight by frequencies
            freq_i = hap_i.frequency
            freq_j = hap_j.frequency

            total_diff += diff * freq_i * freq_j
            comparisons += freq_i * freq_j

    theta_pi = total_diff / comparisons if comparisons > 0 else 0.0

    # Calculate Fs
    # Fs is complex, using simplified approximation
    # Fs ≈ ln(S / (theta_pi)) where S is observed haplotypes

    if theta_pi > 0:
        # Approximate calculation
        expected_k = 1 + theta_pi * sum(1.0 / i for i in range(1, n))

        if expected_k > 0:
            Fs = (
                math.log(k / expected_k)
                if k > expected_k
                else -math.log(expected_k / k)
            )
        else:
            Fs = 0.0
    else:
        Fs = 0.0

    return FuFsResult(Fs=Fs, theta_pi=theta_pi, n_haplotypes=k, n_samples=n)


def calculate_pairwise_fst(
    network: HaplotypeNetwork, pop1: str, pop2: str
) -> FstResult:
    """
    Calculate pairwise FST between two populations.

    FST measures population differentiation based on genetic variance:
    FST = (HT - HS) / HT

    where:
    - HT is the expected heterozygosity in the total population
    - HS is the expected heterozygosity within populations

    FST = 0: no differentiation
    FST = 1: complete differentiation

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.
    pop1 : str
        Name of first population.
    pop2 : str
        Name of second population.

    Returns
    -------
    FstResult
        FstResult with FST and related values.
    """
    # Get haplotypes for each population
    pop1_counts: Dict[str, int] = defaultdict(int)
    pop2_counts: Dict[str, int] = defaultdict(int)

    for hap_id in network.nodes:
        haplotype = network.get_haplotype(hap_id)
        if haplotype is None:
            continue

        freq_info = haplotype.get_frequency_info()
        if pop1 in freq_info.by_population:
            pop1_counts[hap_id] = freq_info.by_population[pop1]
        if pop2 in freq_info.by_population:
            pop2_counts[hap_id] = freq_info.by_population[pop2]

    n1 = sum(pop1_counts.values())
    n2 = sum(pop2_counts.values())

    if n1 == 0 or n2 == 0:
        return FstResult(fst=0.0, pop1=pop1, pop2=pop2, hs=0.0, ht=0.0)

    # Calculate expected heterozygosity for each population
    h1 = 1.0 - sum((count / n1) ** 2 for count in pop1_counts.values())
    h2 = 1.0 - sum((count / n2) ** 2 for count in pop2_counts.values())

    # Average within-population heterozygosity
    hs = (n1 * h1 + n2 * h2) / (n1 + n2)

    # Calculate total heterozygosity
    all_haplotypes = set(pop1_counts.keys()) | set(pop2_counts.keys())
    n_total = n1 + n2

    total_counts = {}
    for hap_id in all_haplotypes:
        total_counts[hap_id] = pop1_counts.get(hap_id, 0) + pop2_counts.get(hap_id, 0)

    ht = 1.0 - sum((count / n_total) ** 2 for count in total_counts.values())

    # Calculate FST
    if ht > 0:
        fst = (ht - hs) / ht
    else:
        fst = 0.0

    return FstResult(fst=fst, pop1=pop1, pop2=pop2, hs=hs, ht=ht)


def calculate_fst_matrix(network: HaplotypeNetwork) -> Dict[Tuple[str, str], float]:
    """
    Calculate pairwise FST for all population pairs.

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.

    Returns
    -------
    Dict[Tuple[str, str], float]
        Dictionary mapping (pop1, pop2) tuples to FST values.
    """
    # Get all populations
    populations = set()
    for hap_id in network.nodes:
        haplotype = network.get_haplotype(hap_id)
        if haplotype is not None:
            freq_info = haplotype.get_frequency_info()
            populations.update(freq_info.by_population.keys())

    populations = sorted(populations)

    # Calculate pairwise FST
    fst_matrix = {}
    for i, pop1 in enumerate(populations):
        for pop2 in populations[i + 1 :]:
            result = calculate_pairwise_fst(network, pop1, pop2)
            fst_matrix[(pop1, pop2)] = result.fst
            fst_matrix[(pop2, pop1)] = result.fst  # Symmetric

    # Add diagonal (self-comparison = 0)
    for pop in populations:
        fst_matrix[(pop, pop)] = 0.0

    return fst_matrix


def calculate_amova(
    network: HaplotypeNetwork,
    alignment: Alignment,
    groups: Optional[Dict[str, str]] = None,
) -> AMOVAResult:
    """
    Calculate AMOVA (Analysis of Molecular Variance).

    AMOVA partitions genetic variance into components:
    - Among populations
    - Within populations
    - (Optionally) Among groups of populations

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.
    alignment : Alignment
        Alignment object.
    groups : Dict[str, str], optional
        Optional dict mapping population -> group.

    Returns
    -------
    AMOVAResult
        AMOVAResult with variance components and phi statistics.
    """
    # Get populations and their samples
    pop_haplotypes: Dict[str, List[str]] = defaultdict(list)

    for hap_id in network.nodes:
        haplotype = network.get_haplotype(hap_id)
        if haplotype is None:
            continue

        freq_info = haplotype.get_frequency_info()
        for pop, count in freq_info.by_population.items():
            # Add haplotype multiple times according to count
            pop_haplotypes[pop].extend([hap_id] * count)

    populations = list(pop_haplotypes.keys())
    n_pops = len(populations)

    if n_pops < 2:
        return AMOVAResult(
            phi_st=0.0,
            phi_sc=0.0,
            phi_ct=0.0,
            variance_among_pops=0.0,
            variance_within_pops=0.0,
            percent_among_pops=0.0,
            percent_within_pops=0.0,
        )

    # Calculate pairwise distances between all haplotypes
    def get_distance(hap_id1: str, hap_id2: str) -> int:
        """
        Get number of differences between two haplotypes.

        Parameters
        ----------
        hap_id1 : str
            First haplotype id.
        hap_id2 : str
            Second haplotype id.

        Returns
        -------
        int
            Number of differing sites between the two haplotypes.
        """
        if hap_id1 == hap_id2:
            return 0

        hap1 = network.get_haplotype(hap_id1)
        hap2 = network.get_haplotype(hap_id2)

        if hap1 is None or hap2 is None:
            return 0

        seq1 = hap1.sequence.data
        seq2 = hap2.sequence.data

        return sum(1 for a, b in zip(seq1, seq2) if a != b and a != '-' and b != '-')

    # Calculate total sum of squares (SST)
    sum(len(haps) for haps in pop_haplotypes.values())

    # Mean pairwise distance over all samples
    all_distances = []
    for pop1_haps in pop_haplotypes.values():
        for i, hap1 in enumerate(pop1_haps):
            for hap2 in pop1_haps[i + 1 :]:
                all_distances.append(get_distance(hap1, hap2))

    for i, pop1 in enumerate(populations):
        for pop2 in populations[i + 1 :]:
            for hap1 in pop_haplotypes[pop1]:
                for hap2 in pop_haplotypes[pop2]:
                    all_distances.append(get_distance(hap1, hap2))

    mean_dist_total = np.mean(all_distances) if all_distances else 0.0

    # Calculate within-population variance
    within_pop_variance = 0.0
    for pop_haps in pop_haplotypes.values():
        if len(pop_haps) < 2:
            continue

        pop_distances = []
        for i, hap1 in enumerate(pop_haps):
            for hap2 in pop_haps[i + 1 :]:
                pop_distances.append(get_distance(hap1, hap2))

        if pop_distances:
            within_pop_variance += np.sum(
                [(d - mean_dist_total) ** 2 for d in pop_distances]
            )

    # Calculate among-population variance
    among_pop_variance = 0.0
    for i, pop1 in enumerate(populations):
        for pop2 in populations[i + 1 :]:
            # Mean distance between populations
            between_distances = []
            for hap1 in pop_haplotypes[pop1]:
                for hap2 in pop_haplotypes[pop2]:
                    between_distances.append(get_distance(hap1, hap2))

            if between_distances:
                mean_between = np.mean(between_distances)
                among_pop_variance += (mean_between - mean_dist_total) ** 2 * len(
                    between_distances
                )

    # Normalize variances
    total_variance = within_pop_variance + among_pop_variance

    if total_variance > 0:
        percent_within = (within_pop_variance / total_variance) * 100
        percent_among = (among_pop_variance / total_variance) * 100
        phi_st = among_pop_variance / total_variance
    else:
        percent_within = 0.0
        percent_among = 0.0
        phi_st = 0.0

    return AMOVAResult(
        phi_st=phi_st,
        phi_sc=0.0,  # Not calculated without groups
        phi_ct=0.0,  # Not calculated without groups
        variance_among_pops=among_pop_variance,
        variance_within_pops=within_pop_variance,
        percent_among_pops=percent_among,
        percent_within_pops=percent_within,
    )


def calculate_mismatch_distribution(network: HaplotypeNetwork) -> Dict[int, int]:
    """
    Calculate the mismatch distribution.

    The mismatch distribution shows the frequency of pairwise differences
    between haplotypes. The shape of this distribution can indicate
    demographic history:
    - Unimodal: recent population expansion
    - Multimodal: population at equilibrium

    Parameters
    ----------
    network : HaplotypeNetwork
        HaplotypeNetwork object.

    Returns
    -------
    Dict[int, int]
        Dictionary mapping number of differences -> frequency.
    """
    mismatch_counts: Dict[int, int] = defaultdict(int)

    haplotype_list = network.haplotypes

    for i, hap_i in enumerate(haplotype_list):
        for hap_j in haplotype_list[i:]:  # Include self-comparison (0 differences)
            seq_i = hap_i.sequence.data
            seq_j = hap_j.sequence.data

            differences = sum(
                1 for a, b in zip(seq_i, seq_j) if a != b and a != '-' and b != '-'
            )

            # Weight by product of frequencies
            freq_i = hap_i.frequency
            freq_j = hap_j.frequency

            if i == haplotype_list.index(hap_j):
                # Self-comparison
                weight = freq_i * (freq_i - 1) // 2 if freq_i > 1 else 0
            else:
                weight = freq_i * freq_j

            mismatch_counts[differences] += weight

    return dict(mismatch_counts)
