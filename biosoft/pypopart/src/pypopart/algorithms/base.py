"""Base classes for network construction algorithms in PyPopART."""

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional

from ..core.alignment import Alignment
from ..core.distance import (
    DistanceMatrix,
    calculate_pairwise_distances,
    pairwise_distance_matrix,
)
from ..core.graph import HaplotypeNetwork


class NetworkAlgorithm(ABC):
    """
    Abstract base class for haplotype network construction algorithms.

    All network construction algorithms should inherit from this class
    and implement the construct_network method.

    Parameters
    ----------
    distance_method : str, default='hamming'
        Method for calculating distances (hamming, p, jc, k2p, tn).
    ignore_gaps : bool, default=True
        Whether to ignore gap positions when calculating distances.
    **kwargs : dict
        Not accepted; present only to give a clear error for typos.
    """

    def __init__(
        self, distance_method: str = 'hamming', ignore_gaps: bool = True, **kwargs
    ):
        """
        Initialize network algorithm.

        Parameters
        ----------
        distance_method : str, default='hamming'
            Method for calculating distances (hamming, p, jc, k2p, tn).
        ignore_gaps : bool, default=True
            Whether to ignore gap positions when calculating distances.
        **kwargs : dict
            Not accepted; present only to give a clear error for typos.

        Raises
        ------
        TypeError
            If unknown keyword arguments are passed (e.g. a misspelled
            parameter name, which would otherwise be silently ignored).
        """
        if kwargs:
            raise TypeError(
                f'{type(self).__name__} got unexpected keyword argument(s): '
                f'{", ".join(sorted(kwargs))}'
            )
        self.distance_method = distance_method
        self.params: Dict[str, Any] = {'ignore_gaps': ignore_gaps}
        self._distance_matrix: Optional[DistanceMatrix] = None

    @abstractmethod
    def construct_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
        Construct haplotype network from sequence alignment.

        Parameters
        ----------
        alignment : Alignment
            Multiple sequence alignment.
        distance_matrix : DistanceMatrix, optional
            Optional pre-computed distance matrix.

        Returns
        -------
        HaplotypeNetwork
            Constructed haplotype network.
        """
        pass

    def calculate_distances(self, alignment: Alignment) -> DistanceMatrix:
        """
        Calculate pairwise distances between sequences.

        Parameters
        ----------
        alignment : Alignment
            Multiple sequence alignment.

        Returns
        -------
        DistanceMatrix
            Pairwise distance matrix for the alignment.
        """
        return calculate_pairwise_distances(
            alignment,
            method=self.distance_method,
            ignore_gaps=self.params.get('ignore_gaps', True),
        )

    def calculate_haplotype_distances(self, haplotypes) -> DistanceMatrix:
        """
        Calculate pairwise distances between unique haplotypes.

        The single shared distance path for all algorithms: honours the
        configured distance_method and ignore_gaps, and uses the fast
        whole-matrix kernel for Hamming distances.

        Parameters
        ----------
        haplotypes : list of Haplotype
            Unique haplotypes to compare.

        Returns
        -------
        DistanceMatrix
            Pairwise distance matrix labelled by haplotype id.
        """
        return pairwise_distance_matrix(
            haplotypes,
            method=self.distance_method,
            ignore_gaps=self.params.get('ignore_gaps', True),
        )

    def build_network(
        self, alignment: Alignment, distance_matrix: Optional[DistanceMatrix] = None
    ) -> HaplotypeNetwork:
        """
        Build haplotype network from sequence alignment.

        This is an alias for construct_network() to maintain backward compatibility
        with CLI and GUI code.

        Parameters
        ----------
        alignment : Alignment
            Multiple sequence alignment.
        distance_matrix : DistanceMatrix, optional
            Optional pre-computed distance matrix.

        Returns
        -------
        HaplotypeNetwork
            Constructed haplotype network.
        """
        return self.construct_network(alignment, distance_matrix)

    def get_parameters(self) -> Dict[str, Any]:
        """
        Get algorithm parameters.

        Returns
        -------
        dict
            Dictionary of algorithm parameters.
        """
        return {'distance_method': self.distance_method, **self.params}

    def __str__(self) -> str:
        """Return string representation."""
        return f'{self.__class__.__name__}(distance={self.distance_method})'

    def __repr__(self) -> str:
        """Detailed representation."""
        params_str = ', '.join(f'{k}={v}' for k, v in self.get_parameters().items())
        return f'{self.__class__.__name__}({params_str})'
