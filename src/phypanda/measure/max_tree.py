"""
Max-tree diversity measure implementation.
"""

from __future__ import annotations

from typing import Any, Mapping, Set

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooNotImplementedError

class MaxTreeDiversity:
    """
    Maximum displayed-tree diversity measure.
    """

    def compute_diversity(
        self,
        network: DirectedPhyNetwork,
        taxa: Set[str],
        **kwargs: Any,
    ) -> float:
        """
        Compute maximum diversity over all displayed trees.

        Parameters
        ----------
        network : DirectedPhyNetwork
            Input phylogenetic network.
        taxa : Set[str]
            Selected taxa.
        **kwargs : Any
            Unused measure-specific keyword arguments.

        Returns
        -------
        float
            Maximum all-paths diversity value across displayed trees.
        """
        raise PhyloZooNotImplementedError(
            "MaxTreeDiversity.compute_diversity is not implemented yet."
        )

    def solve_maximization(
        self,
        network: DirectedPhyNetwork,
        budget: int,
        costs: Mapping[str, int] | None = None,
        **kwargs: Any,
    ) -> tuple[float, Set[str]]:
        """
        Solve budgeted MaxTreePD maximization.

        Notes
        -----
        This exact optimizer is not implemented yet.
        """
        raise PhyloZooNotImplementedError(
            "MaxTreeDiversity.solve_maximization is not implemented yet."
        )


max_tree = MaxTreeDiversity()
