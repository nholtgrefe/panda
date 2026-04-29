"""Average displayed-tree diversity stubs."""

from __future__ import annotations

from typing import Any, Mapping, Set

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooNotImplementedError


class AverageTreeDiversity:
    """
    Average displayed-tree diversity measure stub.

    Examples
    --------
    >>> import phypanda as pp
    >>> # pp.average_tree.compute_diversity(network, {"a"})
    """

    def compute_diversity(
        self,
        network: DirectedPhyNetwork,
        taxa: Set[str],
        **kwargs: Any,
    ) -> float:
        """
        Compute average displayed-tree diversity.

        Parameters
        ----------
        network : DirectedPhyNetwork
            Input phylogenetic network.
        taxa : Set[str]
            Selected taxa.
        **kwargs : Any
            Unused measure-specific options.

        Returns
        -------
        float
            Average displayed-tree diversity value.
        """
        raise PhyloZooNotImplementedError(
            "AverageTreeDiversity.compute_diversity is not implemented yet."
        )

    def solve_maximization(
        self,
        network: DirectedPhyNetwork,
        budget: int,
        costs: Mapping[str, int] | None = None,
        **kwargs: Any,
    ) -> tuple[float, Set[str]]:
        """
        Solve average displayed-tree diversity maximization.

        Parameters
        ----------
        network : DirectedPhyNetwork
            Input phylogenetic network.
        budget : int
            Integer budget.
        costs : Mapping[str, int] | None, optional
            Optional taxon costs.
        **kwargs : Any
            Unused optimization options.

        Returns
        -------
        tuple[float, Set[str]]
            Objective value and selected taxa.
        """
        raise PhyloZooNotImplementedError(
            "AverageTreeDiversity.solve_maximization is not implemented yet."
        )


average_tree = AverageTreeDiversity()
