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

    Notes
    -----
    Concrete measures implement the operations that the public helpers in
    :mod:`phypanda.base` call on their behalf.
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

        Examples
        --------
        >>> # value = pp.compute_diversity(network, {"a", "b"}, measure=pp.all_paths)
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

        Examples
        --------
        >>> # value, taxa = pp.solve_max_diversity(network, budget=5, measure=pp.all_paths)
        """
        raise PhyloZooNotImplementedError(
            "This measure does not implement custom optimization."
        )
