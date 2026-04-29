"""Tree-related diversity measure stubs."""

from __future__ import annotations

from typing import Any, Mapping, Set

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooNotImplementedError


class TreeDiversity:
    """Tree-only diversity measure stub."""

    def compute_diversity(
        self,
        network: DirectedPhyNetwork,
        taxa: Set[str],
        **kwargs: Any,
    ) -> float:
        """Compute diversity on trees only."""
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
        """Solve tree-only diversity maximization."""
        raise PhyloZooNotImplementedError(
            "TreeDiversity.solve_maximization is not implemented yet."
        )


tree = TreeDiversity()
