"""
Min-tree diversity measure implementation.
"""

from __future__ import annotations

from typing import Any, Mapping, Set

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooNotImplementedError

class MinTreeDiversity:
    """
    Minimum displayed-tree diversity measure.

    Examples
    --------
    >>> import phypanda as pp
    >>> # pp.min_tree.compute_diversity(network, {"a"})
    """

    def compute_diversity(
        self,
        network: DirectedPhyNetwork,
        taxa: Set[str],
        **kwargs: Any,
    ) -> float:
        """
        Compute minimum diversity over all displayed trees.

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
            Minimum all-paths diversity value across displayed trees.

        Examples
        --------
        >>> import phypanda as pp
        >>> # pp.min_tree.compute_diversity(network, {"a", "b"})  # doctest: +SKIP
        """
        raise PhyloZooNotImplementedError(
            "MinTreeDiversity.compute_diversity is not implemented yet."
        )

    def solve_maximization(
        self,
        network: DirectedPhyNetwork,
        budget: int,
        costs: Mapping[str, int] | None = None,
        **kwargs: Any,
    ) -> tuple[float, Set[str]]:
        """
        Solve budgeted MinTreePD maximization.

        Notes
        -----
        This exact optimizer is not implemented yet.

        Examples
        --------
        >>> import phypanda as pp
        >>> # pp.min_tree.solve_maximization(network, budget=5)  # doctest: +SKIP
        """
        raise PhyloZooNotImplementedError(
            "MinTreeDiversity.solve_maximization is not implemented yet."
        )


min_tree = MinTreeDiversity()
