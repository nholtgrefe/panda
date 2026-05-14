"""All-paths diversity measure implementation."""

from __future__ import annotations

from typing import Any, Mapping, Set

import networkx as nx
from scanwidth import TreeExtension
from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooNotImplementedError, PhyloZooValueError

from .all_paths_solvers import solve_esw_fpt, solve_nsw_fpt_budget


def _all_paths_diversity_for_taxa(
    network: DirectedPhyNetwork,
    taxa: Set[str],
) -> float:
    """
    Compute all-paths diversity for a fixed taxon set.

    Parameters
    ----------
    network : DirectedPhyNetwork
        Input phylogenetic network.
    taxa : Set[str]
        Taxa to retain.

    Returns
    -------
    float
        Total branch length of the induced ancestral subnetwork.
    """
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
    """
    All-paths diversity measure.

    Examples
    --------
    >>> import phypanda as pp
    >>> # value = pp.compute_diversity(network, {"a", "b"}, measure=pp.all_paths)
    """

    _ALGORITHM_REGISTRY = {
        "esw_fpt": solve_esw_fpt,
        "nsw_fpt_budget": solve_nsw_fpt_budget,
    }

    def compute_diversity(
        self,
        network: DirectedPhyNetwork,
        taxa: Set[str],
        **kwargs: Any,
    ) -> float:
        """
        Compute all-paths diversity for a taxon set.

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
            All-paths diversity value.

        Examples
        --------
        >>> import phypanda as pp
        >>> # value = pp.compute_diversity(network, {"a", "b"}, measure=pp.all_paths)
        """
        return _all_paths_diversity_for_taxa(network, taxa)

    def solve_maximization(
        self,
        network: DirectedPhyNetwork,
        budget: int,
        costs: Mapping[str, int] | None = None,
        algorithm: str = "nsw_fpt_budget",
        tree_extension: TreeExtension | None = None,
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
            Solver backend name. Default is ``"nsw_fpt_budget"``.
        tree_extension : TreeExtension | None, optional
            Optional precomputed tree extension for compatible solvers.
            If ``None``, the solver computes one using ``scanwidth`` defaults.
        **kwargs : Any
            Solver-specific options forwarded to the selected backend (for
            example ``numba=False`` with ``algorithm="nsw_fpt_budget"`` to
            disable Numba in the child-merge DP).

        Returns
        -------
        tuple[float, Set[str]]
            Objective value and selected taxa.

        Examples
        --------
        >>> import phypanda as pp
        >>> # value, taxa = pp.solve_max_diversity(network, budget=5, measure=pp.all_paths)
        >>> # value2, taxa2 = pp.solve_max_diversity(
        ... #     network, budget=5, algorithm="esw_fpt", tree_extension=None
        ... # )
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
            tree_extension=tree_extension,
            **kwargs,
        )


all_paths = AllPathsDiversity()
