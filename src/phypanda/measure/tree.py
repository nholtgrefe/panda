"""Tree-related diversity measure stubs."""

from __future__ import annotations

from typing import Any, Mapping, Set

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooNotImplementedError


class TreeDiversity:
    """
    Tree-only diversity measure stub.

    Examples
    --------
    >>> import phypanda as pp
    >>> # pp.compute_diversity(network, {"a"}, measure=pp.tree)
    """

    def compute_diversity(
        self,
        network: DirectedPhyNetwork,
        taxa: Set[str],
        **kwargs: Any,
    ) -> float:
        """
        Compute diversity on trees only.

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
            Tree-only diversity value.
        """
        raise PhyloZooNotImplementedError(
            "TreeDiversity.compute_diversity is not implemented yet."
        )

    def solve_maximization(
        self,
        network: DirectedPhyNetwork,
        budget: int,
        costs: Mapping[str, int] | None = None,
        **kwargs: Any,
    ) -> tuple[float, Set[str]]:
        """
        Solve tree-only diversity maximization.

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
            "TreeDiversity.solve_maximization is not implemented yet."
        )


tree = TreeDiversity()
