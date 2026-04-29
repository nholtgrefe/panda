"""
Protocol defining the interface for diversity measures.
"""

from __future__ import annotations

from typing import Any, Mapping, Protocol, Set

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooNotImplementedError


class DiversityMeasure(Protocol):
    """
    Protocol defining the interface for diversity measures.
    """

    def compute_diversity(
        self,
        network: DirectedPhyNetwork,
        taxa: Set[str],
        **kwargs: Any,
    ) -> float:
        """
        Compute diversity for a fixed taxon set.

        Parameters
        ----------
        network : DirectedPhyNetwork
            Input phylogenetic network.
        taxa : Set[str]
            Selected taxa.
        **kwargs : Any
            Measure-specific keyword arguments.

        Returns
        -------
        float
            Diversity value.
        """
        ...

    def solve_maximization(
        self,
        network: DirectedPhyNetwork,
        budget: int,
        costs: Mapping[str, int] | None = None,
        **kwargs: Any,
    ) -> tuple[float, Set[str]]:
        """
        Solve diversity maximization under budget constraints.

        Parameters
        ----------
        network : DirectedPhyNetwork
            Input phylogenetic network.
        budget : int
            Total integer budget.
        costs : Mapping[str, int] | None, optional
            Integer costs per taxon. If ``None``, unit costs are assumed.
        **kwargs : Any
            Measure-specific keyword arguments.

        Returns
        -------
        tuple[float, Set[str]]
            Objective value and selected taxa.

        Raises
        ------
        PhyloZooNotImplementedError
            If the concrete measure does not support optimization.
        """
        raise PhyloZooNotImplementedError(
            "This measure does not implement custom optimization."
        )
