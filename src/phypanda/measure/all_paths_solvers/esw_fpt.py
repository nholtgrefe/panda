"""
Edge-scanwidth FPT solver backend for all-paths diversity.
"""

from __future__ import annotations

from collections import defaultdict
import math
from typing import Any, Mapping, Set

import networkx as nx
from scanwidth import DAG, TreeExtension, edge_scanwidth

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.classifications import has_parallel_edges, is_binary
from phylozoo.core.network.dnetwork.transformations import binary_resolution
from phylozoo.utils.exceptions import (
    PhyloZooNotImplementedError,
    PhyloZooRuntimeError,
    PhyloZooValueError,
)

from ...utils import powerset


class _DPInstance:
    """Dynamic-programming instance for unit-cost all-paths MAPPD."""

    def __init__(self, network: DirectedPhyNetwork, tree_extension: TreeExtension) -> None:
        """Initialize dynamic-programming state."""
        self.network = network
        self.tree_extension = tree_extension
        self.root = network.root_node
        self.graph = network._graph

        self.edge_pairs: list[tuple[Any, Any]] = []
        total_weight = 0.0
        for u, v, _ in self.graph.edges(keys=True):
            self.edge_pairs.append((u, v))
            total_weight += self.network.get_branch_length(u, v) or 1.0
        self.minus_infinity = -2 * total_weight - 1

        self.GW = self._initialize_GW()
        self.table: dict[
            Any,
            dict[int, dict[frozenset[tuple[Any, Any]], tuple[float, Any]]],
        ] = defaultdict(lambda: defaultdict(dict))
        self.m = self._max_edge_offspring_count()

    @staticmethod
    def _solution_set_key(solution_set: frozenset[tuple[Any, Any]]) -> tuple[tuple[str, str], ...]:
        """Build a deterministic comparable key for a solution edge-set."""
        return tuple(sorted((str(e[0]), str(e[1])) for e in solution_set))

    def _get_branch_length(self, u: Any, v: Any) -> float:
        """Return branch length on edge ``(u, v)`` with default 1.0."""
        return self.network.get_branch_length(u, v) or 1.0

    def _initialize_GW(self) -> dict[Any, set[tuple[Any, Any]]]:
        """Build scanwidth bags ``GW_v`` for all vertices ``v``."""
        res: dict[Any, set[tuple[Any, Any]]] = {v: set() for v in self.graph.nodes()}

        for v in self.graph.nodes():
            sink = set(nx.descendants(self.tree_extension.tree, v))
            sink.add(v)

            for u, w in self.edge_pairs:
                if u not in sink and w in sink:
                    res[v].add((u, w))

        return res

    def _max_edge_offspring_count(self) -> int:
        """Compute maximum leaf edge-offspring count used for pruning."""
        indeg = {v: self.network.indegree(v) for v in self.graph.nodes()}
        edge_count: dict[Any, int] = {}

        def count_edges(u: Any) -> int:
            if u in edge_count:
                return edge_count[u]
            total = indeg[u]
            for v in self.network.parents(u):
                total += count_edges(v)
            edge_count[u] = total
            return total

        return max((count_edges(leaf) for leaf in self.network.leaves), default=0)

    def _fill_dp_table(self, k: int) -> None:
        """Fill DP table bottom-up for budget ``k`` in unit-cost mode."""
        self.table = defaultdict(lambda: defaultdict(dict))
        leaves = self.network.leaves
        for v in nx.dfs_postorder_nodes(self.tree_extension.tree):
            if v in leaves:
                self._process_leaf_node(v, k)
            else:
                self._process_non_leaf_node(v, k)

    def _process_leaf_node(self, v: Any, k: int) -> None:
        """Fill DP entries for a leaf node."""
        GW_v = self.GW.get(v, set())
        for phi in powerset(GW_v):
            phi_frozen = frozenset(phi)
            for l in range(k + 1):
                if l == 0:
                    dp = 0.0 if len(phi) == 0 else self.minus_infinity
                    pointer = None
                else:
                    dp = sum(self._get_branch_length(u, w) for (u, w) in phi)
                    pointer = (0,)
                self.table[v][l][phi_frozen] = (dp, pointer)

    def _process_non_leaf_node(self, v: Any, k: int) -> None:
        """Fill DP entries for a non-leaf node."""
        GW_v = self.GW.get(v, set())
        delta_in_v = {
            (u, v) for u, _, _ in self.network.incident_parent_edges(v, keys=True)
        }
        delta_out_v = {
            (v, w) for _, w, _ in self.network.incident_child_edges(v, keys=True)
        }

        children = list(self.tree_extension.tree.successors(v))
        GW_children = [self.GW.get(child, set()) for child in children]

        for phi in powerset(GW_v):
            phi_frozen = frozenset(phi)
            phi_len = len(phi)
            phi_delta_in_v = phi & delta_in_v
            bound = math.ceil(phi_len / self.m) if self.m > 0 else 0
            val3 = sum(self._get_branch_length(u, w) for (u, w) in phi_delta_in_v)
            phi_psi_subsets = [phi | psi for psi in powerset(delta_out_v) if len(psi) > 0]
            phi_psi_GW_subsets: list[list[frozenset[tuple[Any, Any]]]] = []
            for i in range(len(children)):
                GW = GW_children[i]
                child_subsets: list[frozenset[tuple[Any, Any]]] = []
                for phi_psi in phi_psi_subsets:
                    child_subsets.append(frozenset(phi_psi & GW))
                phi_psi_GW_subsets.append(child_subsets)

            for l in range(k + 1):
                if l == 0:
                    dp = 0.0 if phi_len == 0 else self.minus_infinity
                    pointer = None
                elif l < bound:
                    dp = self.minus_infinity
                    pointer = None
                else:
                    dp = self.minus_infinity - 1
                    pointer = None
                    if len(children) == 1:
                        u = children[0]
                        if len(phi_delta_in_v) == 0 and v != self.root:
                            S = frozenset(phi & GW_children[0])
                            dp = self.table.get(u, {}).get(l, {}).get(
                                S, (self.minus_infinity, None)
                            )[0]
                            pointer = (1, (u, l, S))
                        else:
                            sorted_subsets = sorted(
                                phi_psi_GW_subsets[0],
                                key=lambda s: (self._solution_set_key(s), id(s)),
                            )
                            for S in sorted_subsets:
                                val1 = self.table.get(u, {}).get(l, {}).get(
                                    S, (self.minus_infinity, None)
                                )[0]
                                val = val1 + val3
                                should_update = False
                                if val > dp:
                                    should_update = True
                                elif val == dp and pointer is not None:
                                    current_S = pointer[1][2]
                                    if self._solution_set_key(S) < self._solution_set_key(current_S):
                                        should_update = True
                                elif val == dp and pointer is None:
                                    should_update = True
                                if should_update:
                                    dp = val
                                    pointer = (1, (u, l, S))

                    elif len(children) == 2:
                        u, w = children
                        if len(phi_delta_in_v) == 0 and v != self.root:
                            S1 = frozenset(phi & GW_children[0])
                            S2 = frozenset(phi & GW_children[1])
                            for l_prime in range(l + 1):
                                val1 = self.table.get(u, {}).get(l_prime, {}).get(
                                    S1, (self.minus_infinity, None)
                                )[0]
                                val2 = self.table.get(w, {}).get(
                                    l - l_prime, {}
                                ).get(S2, (self.minus_infinity, None))[0]
                                val = val1 + val2
                                should_update = False
                                if val > dp:
                                    should_update = True
                                elif val == dp and pointer is not None:
                                    current_S1 = pointer[1][2]
                                    current_S2 = pointer[2][2]
                                    if (
                                        self._solution_set_key(S1) < self._solution_set_key(current_S1)
                                        or self._solution_set_key(S2) < self._solution_set_key(current_S2)
                                    ):
                                        should_update = True
                                elif val == dp and pointer is None:
                                    should_update = True
                                if should_update:
                                    dp = val
                                    pointer = (2, (u, l_prime, S1), (w, l - l_prime, S2))
                        else:
                            for l_prime in range(l + 1):
                                indices = list(range(len(phi_psi_GW_subsets[0])))
                                indices.sort(
                                    key=lambda i: self._solution_set_key(phi_psi_GW_subsets[0][i])
                                )
                                for i in indices:
                                    S1 = phi_psi_GW_subsets[0][i]
                                    val1 = self.table.get(u, {}).get(l_prime, {}).get(
                                        S1, (self.minus_infinity, None)
                                    )[0]
                                    S2 = phi_psi_GW_subsets[1][i]
                                    val2 = self.table.get(w, {}).get(
                                        l - l_prime, {}
                                    ).get(S2, (self.minus_infinity, None))[0]
                                    val = val1 + val2 + val3
                                    should_update = False
                                    if val > dp:
                                        should_update = True
                                    elif val == dp and pointer is not None:
                                        current_S1 = pointer[1][2]
                                        current_S2 = pointer[2][2]
                                        solution_key = (
                                            self._solution_set_key(S1),
                                            self._solution_set_key(S2),
                                        )
                                        current_solution_key = (
                                            self._solution_set_key(current_S1),
                                            self._solution_set_key(current_S2),
                                        )
                                        if solution_key < current_solution_key:
                                            should_update = True
                                    elif val == dp and pointer is None:
                                        should_update = True
                                    if should_update:
                                        dp = val
                                        pointer = (2, (u, l_prime, S1), (w, l - l_prime, S2))
                self.table[v][l][phi_frozen] = (dp, pointer)

    def _backtrack_solution(
        self,
        v: Any,
        l: int,
        phi_frozen: frozenset[tuple[Any, Any]],
    ) -> set[Any]:
        """Reconstruct selected leaves from DP pointers."""
        pointer = self.table.get(v, {}).get(l, {}).get(phi_frozen, (None, None))[1]
        if pointer is None:
            return set()
        if pointer[0] == 0:
            return {v}
        if pointer[0] == 1:
            child, l_child, phi_child = pointer[1]
            return self._backtrack_solution(child, l_child, phi_child)
        if pointer[0] == 2:
            (u, l1, phi_u) = pointer[1]
            (w, l2, phi_w) = pointer[2]
            return self._backtrack_solution(
                u,
                l1,
                phi_u,
            ) | self._backtrack_solution(w, l2, phi_w)
        return set()

    def solve(self, k: int) -> tuple[float, set[Any]]:
        """Solve the DP and return optimal value and selected leaf node IDs."""
        self._fill_dp_table(k)
        root = self.network.root_node
        pd = self.table[root][k][frozenset()][0]
        solution = self._backtrack_solution(root, k, frozenset())
        return pd, solution


def solve_esw_fpt(
    network: DirectedPhyNetwork,
    budget: int,
    costs: Mapping[str, int] | None = None,
    tree_extension: str | TreeExtension = "xp",
    **kwargs: Any,
) -> tuple[float, Set[str]]:
    """
    Solve unit-cost MAPPD exactly using edge-scanwidth FPT dynamic programming.
    """
    if budget < 0 or budget > len(network.taxa):
        raise PhyloZooValueError(
            f"budget must be between 0 and {len(network.taxa)}, got {budget}"
        )

    if costs is not None and any(cost != 1 for cost in costs.values()):
        raise PhyloZooNotImplementedError(
            "esw_fpt currently supports only unit costs."
        )

    if has_parallel_edges(network):
        raise PhyloZooValueError(
            "MAPPD algorithm cannot be applied to networks with parallel edges"
        )

    working_network = network
    if not is_binary(network):
        working_network = binary_resolution(network)
        if isinstance(tree_extension, TreeExtension):
            tree_extension = "xp"

    if not isinstance(tree_extension, TreeExtension):
        dag = DAG(working_network._graph._graph)
        res = edge_scanwidth(dag, algorithm=tree_extension, **kwargs)
        if res[0] is None:
            raise PhyloZooRuntimeError("Failed to compute tree extension")
        _, extension = res
        tree_extension = extension.to_canonical_tree_extension()

    dp_instance = _DPInstance(working_network, tree_extension)
    pd_value, solution_node_ids = dp_instance.solve(budget)
    solution_taxa: Set[str] = set()
    for node_id in solution_node_ids:
        label = working_network.get_label(node_id)
        if label is None:
            raise PhyloZooRuntimeError(f"Leaf node {node_id} has no label")
        solution_taxa.add(label)

    return pd_value, solution_taxa
