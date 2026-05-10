"""
Node-scanwidth FPT solver backend for budgeted all-paths diversity.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Mapping, Set

import networkx as nx
from scanwidth import DAG, TreeExtension, node_scanwidth

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.classifications import has_parallel_edges
from phylozoo.utils.exceptions import PhyloZooRuntimeError, PhyloZooValueError

from ...utils import _ChildMergeDP, normalize_costs, powerset


class _DPInstance:
    """
    Dynamic-programming instance for the node-scanwidth budgeted recurrence.

    Notes
    -----
    The table state is ``DP[v, Y, b]`` where:

    - ``v`` is a node in the tree extension,
    - ``Y`` is a subset of the node-scanwidth bag ``GW(v)``,
    - ``b`` is a budget-like parameter (``DP1`` or ``DP2`` depending on
      ``self.use_dp2``).

    Internal transitions require combining child states under
    partition/split constraints; that nested operation is delegated to
    :class:`phypanda.utils._ChildMergeDP`.
    """

    def __init__(
        self,
        network: DirectedPhyNetwork,
        tree_extension: TreeExtension,
        budget: int,
        costs: Mapping[str, int],
        *,
        use_numba: bool = True,
    ) -> None:
        self.network = network
        self.tree_extension = tree_extension
        self.root = network.root_node
        self.graph = network._graph
        self.budget = budget
        self.costs = dict(costs)
        self.leaves = set(network.leaves)

        self.edge_pairs: list[tuple[Any, Any]] = []
        total_weight = 0.0
        for u, v, _ in self.graph.edges(keys=True):
            self.edge_pairs.append((u, v))
            total_weight += self.network.get_branch_length(u, v) or 1.0
        self.minus_infinity = -2 * total_weight - 1.0

        self.total_cost = sum(self.costs[taxon] for taxon in network.taxa)
        self.b_bar = self.total_cost - self.budget
        self.b_max = min(self.budget, self.b_bar)
        if self.b_max < 0:
            raise PhyloZooValueError(
                "Budget is larger than total available taxon cost; recurrence undefined."
            )
        self.use_dp2 = self.budget > self.b_bar

        self.leaf_cost: dict[Any, int] = {}
        for leaf in self.leaves:
            label = self.network.get_label(leaf)
            if label is None:
                raise PhyloZooValueError(f"Leaf node {leaf} has no taxon label.")
            self.leaf_cost[leaf] = self.costs[label]

        # Prefer bulk bag extraction from scanwidth for lower Python overhead.
        self.GW = {
            v: frozenset(bag)
            for v, bag in self.tree_extension.node_scanwidth_bags().items()
        }
        self.table: dict[Any, dict[frozenset[Any], dict[int, float]]] = defaultdict(
            lambda: defaultdict(dict)
        )
        self.pointers: dict[Any, dict[frozenset[Any], dict[int, Any]]] = defaultdict(
            lambda: defaultdict(dict)
        )
        self._child_merge_dp = _ChildMergeDP(
            self.minus_infinity,
            use_numba=use_numba,
        )

    def _incoming_weight(self, v: Any) -> float:
        """Return sum of weights of all incoming network edges to ``v``."""
        return sum(
            self.network.get_branch_length(u, v) or 1.0
            for u, _, _ in self.network.incident_parent_edges(v, keys=True)
        )

    def _children_in_tree(self, v: Any) -> list[Any]:
        """Return children of ``v`` in the tree extension."""
        return list(self.tree_extension.tree.successors(v))

    def _dp_get(
        self,
        table: dict[Any, dict[frozenset[Any], dict[int, float]]],
        v: Any,
        y: frozenset[Any],
        b: int,
    ) -> float:
        """Get DP value with ``-inf`` fallback for missing states."""
        return table.get(v, {}).get(y, {}).get(b, self.minus_infinity)

    def _combine_children(
        self,
        table: dict[Any, dict[frozenset[Any], dict[int, float]]],
        v: Any,
        y: frozenset[Any],
        b: int,
    ) -> tuple[float, tuple[tuple[frozenset[Any], ...], tuple[int, ...]] | None]:
        """
        Combine child tables for one parent state.

        Parameters
        ----------
        table : dict[Any, dict[frozenset[Any], dict[int, float]]]
            DP table to read child values from.
        v : Any
            Parent node in the tree extension.
        y : frozenset[Any]
            Parent bag subset that must be partitioned across children.
        b : int
            Parent budget that must be split across children.

        Returns
        -------
        tuple[float, tuple[tuple[frozenset[Any], ...], tuple[int, ...]] | None]
            Best merge value and argmax partition/split choices for backtracking.
        """
        children = self._children_in_tree(v)
        child_tables = [table.get(child, {}) for child in children]
        return self._child_merge_dp.combine(child_tables=child_tables, items=y, budget=b)

    def _process_leaf(self, x: Any) -> None:
        """Fill leaf states for the selected table using base formulas."""
        parents = list(self.network.parents(x))
        if len(parents) != 1:
            raise PhyloZooValueError(
                f"Leaf {x} must have exactly one parent for the leaf base case."
            )
        parent = parents[0]
        y_parent = frozenset({parent})
        edge_weight = self.network.get_branch_length(parent, x) or 1.0
        cost_x = self.leaf_cost[x]

        for y in powerset(set(self.GW[x])):
            y_frozen = frozenset(y)
            for b in range(self.b_max + 1):
                if self.use_dp2:
                    # DP2 leaf base
                    if b == 0:
                        val = edge_weight
                        ptr = ("leaf_include",)
                    else:
                        if y_frozen == y_parent or b > cost_x:
                            val = self.minus_infinity
                            ptr = ("leaf_invalid",)
                        else:
                            val = 0.0
                            ptr = ("leaf_exclude",)
                else:
                    # DP1 leaf base
                    if b >= cost_x:
                        val = edge_weight
                        ptr = ("leaf_include",)
                    else:
                        if y_frozen == y_parent:
                            val = self.minus_infinity
                            ptr = ("leaf_invalid",)
                        else:
                            val = 0.0
                            ptr = ("leaf_exclude",)
                self.table[x][y_frozen][b] = val
                self.pointers[x][y_frozen][b] = ptr

    def _process_internal(self, v: Any) -> None:
        """
        Fill internal-node states from include/exclude recurrence transitions.

        For each state ``(Y, b)``, we compare:

        - excluding ``v``: merge children with ``Y``,
        - including ``v``: merge children with transformed set ``Y_include`` and
          add incoming-edge weight contribution.
        """
        parents_v = set(self.network.parents(v))
        w_in_v = self._incoming_weight(v)

        for y in powerset(set(self.GW[v])):
            y_frozen = frozenset(y)
            y_include = frozenset((set(y_frozen) - parents_v) | {v})

            for b in range(self.b_max + 1):
                # Option 1: v is excluded; children must realize Y directly.
                option1, choice1 = self._combine_children(self.table, v, y_frozen, b)
                # Option 2: v is included; adjust boundary set and add local gain.
                option2_base, choice2 = self._combine_children(self.table, v, y_include, b)
                option2 = (
                    self.minus_infinity
                    if option2_base <= self.minus_infinity
                    else w_in_v + option2_base
                )
                if option2 > option1:
                    self.table[v][y_frozen][b] = option2
                    self.pointers[v][y_frozen][b] = ("include", choice2)
                else:
                    self.table[v][y_frozen][b] = option1
                    self.pointers[v][y_frozen][b] = ("exclude", choice1)

    def _fill_tables(self) -> None:
        """Compute DP tables in postorder of the tree extension."""
        for v in nx.dfs_postorder_nodes(self.tree_extension.tree):
            if v in self.leaves:
                self._process_leaf(v)
            else:
                self._process_internal(v)

    def _backtrack(self, v: Any, y: frozenset[Any], b: int) -> set[Any]:
        """Backtrack selected vertex set ``Z`` from stored pointers."""
        ptr = self.pointers.get(v, {}).get(y, {}).get(b)
        if ptr is None:
            return set()

        tag = ptr[0]
        if tag == "leaf_include":
            return {v}
        if tag in {"leaf_exclude", "leaf_invalid"}:
            return set()
        if tag not in {"include", "exclude"}:
            return set()

        choice = ptr[1]
        if choice is None:
            return {v} if tag == "include" else set()
        y_parts, b_split = choice
        children = self._children_in_tree(v)
        z: set[Any] = {v} if tag == "include" else set()
        for child, y_child, b_child in zip(children, y_parts, b_split):
            z |= self._backtrack(child, y_child, b_child)
        return z

    def solve(self) -> tuple[float, Set[str]]:
        """Solve recurrence and return objective value with selected taxa set."""
        self._fill_tables()
        empty = frozenset()
        final_b = self.b_bar if self.use_dp2 else self.budget
        value = self._dp_get(self.table, self.root, empty, final_b)
        z_vertices = self._backtrack(self.root, empty, final_b)

        selected_taxa: Set[str] = set()
        for v in z_vertices:
            if v in self.leaves:
                label = self.network.get_label(v)
                if label is not None:
                    selected_taxa.add(label)
        return value, selected_taxa


def solve_nsw_fpt_budget(
    network: DirectedPhyNetwork,
    budget: int,
    costs: Mapping[str, int] | None = None,
    tree_extension: TreeExtension | None = None,
    *,
    numba: bool = True,
    **kwargs: Any,
) -> tuple[float, Set[str]]:
    """
    Solve budgeted all-paths maximization using node-scanwidth DP.

    Parameters
    ----------
    network : DirectedPhyNetwork
        Input phylogenetic network.
    budget : int
        Integer budget ``B``.
    costs : Mapping[str, int] | None, optional
        Integer costs per taxon. If ``None`` or missing for some taxa, unit
        costs are used for those taxa.
    tree_extension : TreeExtension | None, optional
        Optional precomputed tree extension. If ``None``, one is computed via
        :func:`scanwidth.node_scanwidth`.
    numba : bool, default=True
        Use Numba-accelerated child-table merges inside the DP.  Set to
        ``False`` to force the pure Python merge (useful for debugging).
    **kwargs : Any
        Keyword arguments forwarded to ``scanwidth.node_scanwidth`` when a
        tree extension is computed.

    Returns
    -------
    tuple[float, Set[str]]
        Objective value and selected taxa.

    Notes
    -----
    This routine computes (or accepts) a node-scanwidth tree extension and then
    runs a bottom-up DP on that tree extension. The child-combination subproblem
    at each internal node is handled by :class:`phypanda.utils._ChildMergeDP`.
    """
    if budget < 0:
        raise PhyloZooValueError(f"budget must be non-negative, got {budget}")
    if has_parallel_edges(network):
        raise PhyloZooValueError(
            "MAPPD algorithm cannot be applied to networks with parallel edges"
        )

    taxa = set(network.taxa)
    normalized_costs = normalize_costs(taxa, costs)

    total_cost = sum(normalized_costs[taxon] for taxon in taxa)
    if budget >= total_cost:
        from ..all_paths import all_paths

        all_taxa = set(taxa)
        return all_paths.compute_diversity(network, all_taxa), all_taxa

    working_network = network

    if tree_extension is None:
        dag = DAG(working_network._graph._graph)
        value, extension = node_scanwidth(dag, **kwargs)
        if value is None:
            raise PhyloZooRuntimeError("Failed to compute node-scanwidth extension")
        tree_extension = extension.to_canonical_tree_extension()
    elif not isinstance(tree_extension, TreeExtension):
        raise PhyloZooValueError("tree_extension must be a TreeExtension or None.")

    dp = _DPInstance(
        network=working_network,
        tree_extension=tree_extension,
        budget=budget,
        costs=normalized_costs,
        use_numba=numba,
    )
    return dp.solve()

