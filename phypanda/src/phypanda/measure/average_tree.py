"""Average displayed-tree diversity stubs."""

from __future__ import annotations

from typing import Any, Mapping, Set

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooNotImplementedError


class AverageTreeDiversity:
    """Average displayed-tree diversity measure stub."""

    def compute_diversity(
        self,
        network: DirectedPhyNetwork,
        taxa: Set[str],
        **kwargs: Any,
    ) -> float:
        """Compute average displayed-tree diversity."""
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
        """Solve average displayed-tree diversity maximization."""
        raise PhyloZooNotImplementedError(
            "AverageTreeDiversity.solve_maximization is not implemented yet."
        )


average_tree = AverageTreeDiversity()
