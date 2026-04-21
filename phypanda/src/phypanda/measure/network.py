"""Network diversity (NetworkPD) stubs."""

from __future__ import annotations

from typing import Any, Mapping, Set

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooNotImplementedError


class NetworkDiversity:
    """NetworkPD diversity measure stub."""

    def compute_diversity(
        self,
        network: DirectedPhyNetwork,
        taxa: Set[str],
        **kwargs: Any,
    ) -> float:
        """Compute NetworkPD diversity."""
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
        """Solve NetworkPD maximization."""
        raise PhyloZooNotImplementedError(
            "NetworkDiversity.solve_maximization is not implemented yet."
        )


network = NetworkDiversity()
