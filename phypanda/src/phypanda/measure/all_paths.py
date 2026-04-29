"""All-paths diversity measure implementation."""

from __future__ import annotations

from typing import Any, Mapping, Set

import networkx as nx
from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooNotImplementedError, PhyloZooValueError

from .all_paths_solvers import solve_esw_fpt


def _all_paths_diversity_for_taxa(
    network: DirectedPhyNetwork,
    taxa: Set[str],
) -> float:
    """Compute all-paths diversity for a fixed taxon set."""
    if not taxa:
        return 0.0

    leaf_nodes: Set[Any] = set()
    for taxon in taxa:
        node_id = network.get_node_id(taxon)
        if node_id is None:
            raise PhyloZooValueError(f"Taxon '{taxon}' not found in network.")
        leaf_nodes.add(node_id)

    dag = network._graph._graph
    nodes_to_keep: Set[Any] = set(leaf_nodes)
    for leaf in leaf_nodes:
        nodes_to_keep.update(nx.ancestors(dag, leaf))

    total_diversity = 0.0
    for u, v, _, data in network._graph.edges(keys=True, data=True):
        if u in nodes_to_keep and v in nodes_to_keep:
            total_diversity += data.get("branch_length", 1.0)
    return total_diversity


class AllPathsDiversity:
    """All-paths diversity measure."""

    _ALGORITHM_REGISTRY = {
        "esw_fpt": solve_esw_fpt,
    }

    def compute_diversity(
        self,
        network: DirectedPhyNetwork,
        taxa: Set[str],
        **kwargs: Any,
    ) -> float:
        """Compute all-paths diversity for a taxon set."""
        return _all_paths_diversity_for_taxa(network, taxa)

    def solve_maximization(
        self,
        network: DirectedPhyNetwork,
        budget: int,
        costs: Mapping[str, int] | None = None,
        algorithm: str = "esw_fpt",
        **kwargs: Any,
    ) -> tuple[float, Set[str]]:
        """
        Solve all-paths diversity maximization under budget constraints.

        Parameters
        ----------
        network : DirectedPhyNetwork
            Input network.
        budget : int
            Integer budget.
        costs : Mapping[str, int] | None, optional
            Taxon costs.
        algorithm : str, optional
            Solver backend name. Default is ``"esw_fpt"``.
        **kwargs : Any
            Solver-specific options (for example ``tree_extension``).
        """
        solver = self._ALGORITHM_REGISTRY.get(algorithm)
        if solver is None:
            available = ", ".join(sorted(self._ALGORITHM_REGISTRY))
            raise PhyloZooNotImplementedError(
                f"Unknown all_paths algorithm '{algorithm}'. Available: {available}"
            )
        return solver(
            network=network,
            budget=budget,
            costs=costs,
            **kwargs,
        )


all_paths = AllPathsDiversity()
