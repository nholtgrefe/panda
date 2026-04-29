"""Network diversity (NetworkPD) stubs."""

from __future__ import annotations

from typing import Any, Mapping, Set

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooNotImplementedError


class NetworkDiversity:
    """
    NetworkPD diversity measure stub.

    Examples
    --------
    >>> import phypanda as pp
    >>> # pp.network.compute_diversity(network, {"a"})
    """

    def compute_diversity(
        self,
        network: DirectedPhyNetwork,
        taxa: Set[str],
        **kwargs: Any,
    ) -> float:
        """
        Compute NetworkPD diversity.

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
            NetworkPD diversity value.
        """
        raise PhyloZooNotImplementedError(
            "NetworkDiversity.compute_diversity is not implemented yet."
        )

    def solve_maximization(
        self,
        network: DirectedPhyNetwork,
        budget: int,
        costs: Mapping[str, int] | None = None,
        **kwargs: Any,
    ) -> tuple[float, Set[str]]:
        """
        Solve NetworkPD maximization.

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
            "NetworkDiversity.solve_maximization is not implemented yet."
        )


network = NetworkDiversity()
