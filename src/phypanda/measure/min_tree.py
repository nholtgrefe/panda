"""
Min-tree diversity measure implementation.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Mapping, Set

import networkx as nx
from scanwidth import DAG, TreeExtension, node_scanwidth

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.derivations import subnetwork as induced_subnetwork
from phylozoo.utils.exceptions import (
    PhyloZooError,
    PhyloZooNotImplementedError,
    PhyloZooRuntimeError,
    PhyloZooValueError,
)

from ..utils import _ChildMergeDP, powerset


def prune_to_taxa(
    network: DirectedPhyNetwork,
    taxa: Set[str],
) -> DirectedPhyNetwork:
    """
    Return the PhyloZoo subnetwork induced by ``taxa``.

    If ``taxa`` is exactly the leaf set of ``network``, returns ``network``
    unchanged (same object).
    """
    if not taxa:
        raise PhyloZooValueError("taxa must be non-empty for MinTreePD.")

    if taxa == network.taxa:
        return network

    try:
        return induced_subnetwork(
            network,
            list(taxa),
            identify_parallel_edges=True,
        )
    except PhyloZooError:
        raise PhyloZooRuntimeError(
            "PhyloZoo failed to construct induced subnetwork for the requested taxa."
        )


def adapt_tree_extension_for_pruned_network(
    working_network: DirectedPhyNetwork,
    tree_extension: TreeExtension,
) -> TreeExtension:
    """
    Derive a tree extension for ``working_network`` from one on a supergraph.

    For each Γ-vertex absent from ``working_network``, connect its Γ-parent to
    each Γ-child and remove that vertex.  The result is a valid tree extension
    for the pruned DAG; node scanwidth is at most that of ``tree_extension``.

    If ``working_network`` already uses the full vertex set of
    ``tree_extension.tree``, the shortcut loop does nothing and the returned
    object wraps a copy of Γ with the working DAG.
    """
    vertices = set(working_network._graph.nodes())
    full_nodes = set(tree_extension.tree.nodes())
    if not vertices.issubset(full_nodes):
        raise PhyloZooValueError(
            "Every vertex of working_network must occur in tree_extension.tree."
        )
    gamma = tree_extension.tree.copy()
    while True:
        outside = [n for n in gamma.nodes() if n not in vertices]
        if not outside:
            break
        leaf_outside = sorted(
            (v for v in outside if gamma.out_degree(v) == 0),
            key=str,
        )
        pick_from = leaf_outside if leaf_outside else sorted(outside, key=str)
        v = pick_from[0]
        preds = list(gamma.predecessors(v))
        if not preds:
            raise PhyloZooValueError(
                "Cannot adapt tree extension: removable tree vertex has no Γ-parent."
            )
        if len(preds) != 1:
            raise PhyloZooRuntimeError("Tree extension vertex with indegree ≠ 1.")
        p = preds[0]
        succs = list(gamma.successors(v))
        gamma.remove_node(v)
        for w in succs:
            gamma.add_edge(p, w)
    if set(gamma.nodes()) != vertices:
        raise PhyloZooRuntimeError(
            "Internal error: adapted Γ-tree vertex set does not match pruned network."
        )
    dag = DAG(working_network._graph._graph)
    return TreeExtension(dag, gamma)


class _DPInstance:
    """Dynamic-programming instance for MinTreePD fixed-set computation."""

    def __init__(
        self,
        network: DirectedPhyNetwork,
        tree_extension: TreeExtension,
    ) -> None:
        self.network = network
        self.tree_extension = tree_extension
        self.root = network.root_node
        self.graph = network._graph
        self.leaves = set(network.leaves)

        total_weight = 0.0
        for u, v, _ in self.graph.edges(keys=True):
            total_weight += self.network.get_branch_length(u, v) or 1.0
        self.plus_infinity = 2 * total_weight + 1.0

        self.GW = self._initialize_GW()
        self.table: dict[Any, dict[frozenset[Any], dict[frozenset[Any], float]]] = defaultdict(
            lambda: defaultdict(dict)
        )
        self._child_merge_dp = _ChildMergeDP(
            infeasible_value=self.plus_infinity,
            objective="min",
        )

    def _initialize_GW(self) -> dict[Any, frozenset[Any]]:
        """Build node-scanwidth bags ``GW(v)`` for all nodes ``v``."""
        gw: dict[Any, frozenset[Any]] = {}
        for v in self.graph.nodes():
            gw[v] = frozenset(self.tree_extension.node_scanwidth_bag(v))
        return gw

    def _children_in_tree(self, v: Any) -> list[Any]:
        """Return children of ``v`` in the tree extension."""
        return list(self.tree_extension.tree.successors(v))

    def _dp_get(self, v: Any, y: frozenset[Any], z: frozenset[Any]) -> float:
        """Get DP value with ``+inf`` fallback for missing states."""
        return self.table.get(v, {}).get(y, {}).get(z, self.plus_infinity)

    def _combine_children(
        self,
        v: Any,
        y: frozenset[Any],
        z: frozenset[Any],
    ) -> float:
        """Compute DP'[v, y, z] using child partition merge."""
        children = self._children_in_tree(v)
        child_tables: list[dict[frozenset[Any], dict[int, float]]] = []
        for child in children:
            z_child = frozenset(set(z) & set(self.GW[child]))
            child_table: dict[frozenset[Any], dict[int, float]] = {}
            for y_child in powerset(set(self.GW[child])):
                y_child_frozen = frozenset(y_child)
                child_table[y_child_frozen] = {0: self._dp_get(child, y_child_frozen, z_child)}
            child_tables.append(child_table)
        value, _choice = self._child_merge_dp.combine(
            child_tables=child_tables,
            items=y,
            budget=0,
        )
        return value

    def _process_leaf(self, x: Any) -> None:
        """Fill leaf table entries from MinTreePD base case."""
        parents = list(self.network.parents(x))
        if len(parents) != 1:
            raise PhyloZooValueError(
                f"Leaf {x} must have exactly one parent for the MinTree base case."
            )
        parent = parents[0]
        edge_weight = self.network.get_branch_length(parent, x) or 1.0

        for z in powerset(set(self.GW[x])):
            z_frozen = frozenset(z)
            for y in powerset(set(z_frozen)):
                y_frozen = frozenset(y)
                if z_frozen == frozenset({parent}):
                    self.table[x][y_frozen][z_frozen] = edge_weight
                else:
                    self.table[x][y_frozen][z_frozen] = self.plus_infinity

    def _build_child_tables_for_z(
        self,
        children: list[Any],
        z: frozenset[Any],
    ) -> list[dict[frozenset[Any], dict[int, float]]]:
        """Build child merge tables for a fixed z (budget is always 0)."""
        child_tables = []
        for child in children:
            z_child = frozenset(set(z) & set(self.GW[child]))
            child_table: dict[frozenset[Any], dict[int, float]] = {}
            for y_child in powerset(set(self.GW[child])):
                y_child_frozen = frozenset(y_child)
                child_table[y_child_frozen] = {0: self._dp_get(child, y_child_frozen, z_child)}
            child_tables.append(child_table)
        return child_tables

    def _process_internal(self, v: Any) -> None:
        """Fill internal-node entries from MinTreePD recurrence.

        Builds child tables once per distinct z value (and once for z2 = z | {v}
        when parents are active), then fills all y ⊆ z states via O(1) array
        lookups instead of per-(y, z) kernel calls.
        """
        parents_v = set(self.network.parents(v))
        z_states = [frozenset(z) for z in powerset(set(self.GW[v]))]
        if v == self.root:
            # Required by theorem text: allow Z={rho} even if GW(root)=empty.
            z_states.append(frozenset({self.root}))
        children = self._children_in_tree(v)

        for z_frozen in z_states:
            # Exclude option: batch merge over the full z universe (once per z).
            child_tables_z = self._build_child_tables_for_z(children, z_frozen)
            excl_vals, _, _, excl_sbm = self._child_merge_dp.combine_full_table(
                child_tables_z, z_frozen, 0
            )
            excl_mbs: dict[frozenset[Any], int] = {s: m for m, s in enumerate(excl_sbm)}

            # Include option: z2 = z | {v} is the same for all parents.
            active_parents = set(z_frozen) & parents_v
            if active_parents:
                z2 = frozenset(set(z_frozen) | {v})
                child_tables_z2 = self._build_child_tables_for_z(children, z2)
                incl_vals, _, _, incl_sbm = self._child_merge_dp.combine_full_table(
                    child_tables_z2, z2, 0
                )
                incl_mbs: dict[frozenset[Any], int] = {s: m for m, s in enumerate(incl_sbm)}
            else:
                incl_vals = incl_mbs = None
                z2 = None

            for y in powerset(set(z_frozen)):
                y_frozen = frozenset(y)

                excl_mask = excl_mbs.get(y_frozen)
                option1 = (
                    float(excl_vals[excl_mask, 0])
                    if excl_mask is not None
                    else self.plus_infinity
                )
                best = option1

                for parent in active_parents:
                    edge_weight = self.network.get_branch_length(parent, v) or 1.0
                    y2_frozen = frozenset((set(y_frozen) - {parent}) | {v})
                    incl_mask = incl_mbs.get(y2_frozen) if incl_mbs is not None else None
                    option2_base = (
                        float(incl_vals[incl_mask, 0])
                        if incl_mask is not None
                        else self.plus_infinity
                    )
                    option2 = edge_weight + option2_base
                    if option2 < best:
                        best = option2

                self.table[v][y_frozen][z_frozen] = best

    def _fill_tables(self) -> None:
        """Compute DP tables in postorder of the tree extension."""
        for v in nx.dfs_postorder_nodes(self.tree_extension.tree):
            if v in self.leaves:
                self._process_leaf(v)
            else:
                self._process_internal(v)

    def solve(self) -> float:
        """Return MinTreePD value from final root state."""
        self._fill_tables()
        return self._dp_get(self.root, frozenset(), frozenset({self.root}))


class MinTreeDiversity:
    """
    Minimum displayed-tree diversity measure.

    Examples
    --------
    >>> import phypanda as pp
    >>> # pp.compute_diversity(network, {"a"}, measure=pp.min_tree)
    """

    def compute_diversity(
        self,
        network: DirectedPhyNetwork,
        taxa: Set[str],
        tree_extension: TreeExtension | None = None,
        **kwargs: Any,
    ) -> float:
        """
        Compute minimum displayed-tree diversity for a fixed taxon set.

        Parameters
        ----------
        network : DirectedPhyNetwork
            Input phylogenetic network.
        taxa : Set[str]
            Selected taxa.
        tree_extension : TreeExtension | None, optional
            Optional tree extension for the **full** ``network``.  If the working
            network is ``network`` itself (no pruning), it is used as-is; otherwise
            it is adapted with :func:`adapt_tree_extension_for_pruned_network`.
        **kwargs : Any
            Keyword arguments forwarded to ``scanwidth.node_scanwidth`` when a
            tree extension is computed.

        Returns
        -------
        float
            Minimum all-paths diversity value across displayed trees.

        Examples
        --------
        >>> import phypanda as pp
        >>> # pp.compute_diversity(network, {"a", "b"}, measure=pp.min_tree)  # doctest: +SKIP
        """
        taxa_set = set(taxa)
        working_network = prune_to_taxa(network, taxa_set)

        if tree_extension is None:
            dag = DAG(working_network._graph._graph)
            value, extension = node_scanwidth(dag, **kwargs)
            if value is None:
                raise PhyloZooRuntimeError("Failed to compute node-scanwidth extension")
            te = extension.to_canonical_tree_extension()
        else:
            if not isinstance(tree_extension, TreeExtension):
                raise PhyloZooValueError("tree_extension must be a TreeExtension or None.")
            full_vertices = set(network._graph.nodes())
            if set(tree_extension.tree.nodes()) != full_vertices:
                raise PhyloZooValueError(
                    "tree_extension must be defined on the full network vertex set."
                )
            if working_network is network:
                te = tree_extension
            else:
                te = adapt_tree_extension_for_pruned_network(
                    working_network,
                    tree_extension,
                )

        dp = _DPInstance(
            network=working_network,
            tree_extension=te,
        )
        return dp.solve()

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
        >>> # pp.solve_max_diversity(network, budget=5, measure=pp.min_tree)  # doctest: +SKIP
        """
        raise PhyloZooNotImplementedError(
            "MinTreeDiversity.solve_maximization is not implemented yet."
        )


min_tree = MinTreeDiversity()
