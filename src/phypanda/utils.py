"""
General utility helpers for the phypanda package.
"""

from __future__ import annotations

from itertools import combinations
from typing import Any, Iterator, Literal, Mapping

import numpy as np
from numba import njit
from phylozoo.utils.exceptions import PhyloZooValueError


@njit(cache=True)
def _combine_core_max_numba(
    child_values: np.ndarray,
    infeasible: float,
    budget: int,
    full_mask: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Numba core for maximizing child merge DP."""
    deg = child_values.shape[0]
    max_mask = full_mask + 1
    values = np.full((max_mask, budget + 1), infeasible, dtype=np.float64)
    values[0, 0] = 0.0

    choices = np.full((deg, max_mask, budget + 1, 4), -1, dtype=np.int64)

    for j in range(deg):
        next_values = np.full((max_mask, budget + 1), infeasible, dtype=np.float64)
        for assigned_mask in range(max_mask):
            for used_budget in range(budget + 1):
                base_value = values[assigned_mask, used_budget]
                if base_value == infeasible:
                    continue
                remaining_mask = full_mask ^ assigned_mask
                remaining_budget = budget - used_budget
                submask = remaining_mask
                while True:
                    for allocated_budget in range(remaining_budget + 1):
                        child_value = child_values[j, submask, allocated_budget]
                        if child_value == infeasible:
                            continue
                        new_mask = assigned_mask | submask
                        new_budget = used_budget + allocated_budget
                        candidate = base_value + child_value
                        if candidate > next_values[new_mask, new_budget]:
                            next_values[new_mask, new_budget] = candidate
                            choices[j, new_mask, new_budget, 0] = assigned_mask
                            choices[j, new_mask, new_budget, 1] = used_budget
                            choices[j, new_mask, new_budget, 2] = submask
                            choices[j, new_mask, new_budget, 3] = allocated_budget
                    if submask == 0:
                        break
                    submask = (submask - 1) & remaining_mask
        values = next_values
    return values, choices


@njit(cache=True)
def _combine_core_min_numba(
    child_values: np.ndarray,
    infeasible: float,
    budget: int,
    full_mask: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Numba core for minimizing child merge DP."""
    deg = child_values.shape[0]
    max_mask = full_mask + 1
    values = np.full((max_mask, budget + 1), infeasible, dtype=np.float64)
    values[0, 0] = 0.0

    choices = np.full((deg, max_mask, budget + 1, 4), -1, dtype=np.int64)

    for j in range(deg):
        next_values = np.full((max_mask, budget + 1), infeasible, dtype=np.float64)
        for assigned_mask in range(max_mask):
            for used_budget in range(budget + 1):
                base_value = values[assigned_mask, used_budget]
                if base_value == infeasible:
                    continue
                remaining_mask = full_mask ^ assigned_mask
                remaining_budget = budget - used_budget
                submask = remaining_mask
                while True:
                    for allocated_budget in range(remaining_budget + 1):
                        child_value = child_values[j, submask, allocated_budget]
                        if child_value == infeasible:
                            continue
                        new_mask = assigned_mask | submask
                        new_budget = used_budget + allocated_budget
                        candidate = base_value + child_value
                        if candidate < next_values[new_mask, new_budget]:
                            next_values[new_mask, new_budget] = candidate
                            choices[j, new_mask, new_budget, 0] = assigned_mask
                            choices[j, new_mask, new_budget, 1] = used_budget
                            choices[j, new_mask, new_budget, 2] = submask
                            choices[j, new_mask, new_budget, 3] = allocated_budget
                    if submask == 0:
                        break
                    submask = (submask - 1) & remaining_mask
        values = next_values
    return values, choices


def powerset(s: set[Any]) -> Iterator[set[Any]]:
    """
    Yield all subsets of a set.

    Parameters
    ----------
    s : set[Any]
        Input set.

    Yields
    ------
    Iterator[set[Any]]
        Subsets of ``s``.

    Examples
    --------
    >>> list(powerset({"a", "b"}))
    [set(), {'a'}, {'b'}, {'a', 'b'}]
    """
    s_list = tuple(sorted(s, key=str))
    for r in range(len(s_list) + 1):
        for subset in combinations(s_list, r):
            yield set(subset)


def normalize_costs(
    taxa: set[str],
    costs: Mapping[str, int] | None,
) -> dict[str, int]:
    """
    Normalize and validate integer taxon costs.

    Parameters
    ----------
    taxa : set[str]
        Taxon universe to normalize costs for.
    costs : Mapping[str, int] | None
        Optional taxon cost mapping. If ``None``, unit costs are used. Missing
        taxa in the mapping are assigned unit cost.

    Returns
    -------
    dict[str, int]
        Cost mapping for all taxa in ``taxa``.

    Raises
    ------
    PhyloZooValueError
        If a provided cost is not an integer or is negative.
    """
    if costs is None:
        return {taxon: 1 for taxon in taxa}

    normalized: dict[str, int] = {}
    for taxon in taxa:
        cost = costs.get(taxon, 1)
        if not isinstance(cost, int):
            raise PhyloZooValueError(
                f"Cost for taxon '{taxon}' must be an integer, got {type(cost).__name__}"
            )
        if cost < 0:
            raise PhyloZooValueError(
                f"Cost for taxon '{taxon}' must be non-negative, got {cost}"
            )
        normalized[taxon] = cost
    return normalized


class _ChildMergeDP:
    """
    Reusable DP engine for partition-and-budget child-table merges.

    Parameters
    ----------
    infeasible_value : float
        Sentinel value that denotes infeasible states.
    objective : {"max", "min"}, default="max"
        Optimization direction for the merge recurrence.
    use_numba : bool, default=True
        If ``True``, use the Numba-accelerated merge kernels; otherwise use the
        pure Python implementation (same semantics, slower).
    """

    def __init__(
        self,
        infeasible_value: float,
        objective: Literal["max", "min"] = "max",
        *,
        use_numba: bool = True,
    ) -> None:
        self.infeasible_value = infeasible_value
        self.objective = objective
        self._use_numba = use_numba
        self._subset_mask_cache: dict[frozenset[Any], tuple[int, list[frozenset[Any]]]] = {}
        if objective not in {"max", "min"}:
            raise ValueError(f"objective must be 'max' or 'min', got {objective!r}")

    def _subset_mask_data(self, items: frozenset[Any]) -> tuple[int, list[frozenset[Any]]]:
        """Return subset-mask encoding for a fixed item set."""
        cached = self._subset_mask_cache.get(items)
        if cached is not None:
            return cached
        elements = tuple(sorted(items, key=str))
        full_mask = (1 << len(elements)) - 1
        subset_by_mask: list[frozenset[Any]] = [frozenset() for _ in range(full_mask + 1)]
        for mask in range(full_mask + 1):
            subset_by_mask[mask] = frozenset(
                elements[idx] for idx in range(len(elements)) if (mask >> idx) & 1
            )
        out = (full_mask, subset_by_mask)
        self._subset_mask_cache[items] = out
        return out

    def combine(
        self,
        child_tables: list[dict[frozenset[Any], dict[int, float]]],
        items: frozenset[Any],
        budget: int,
    ) -> tuple[float, tuple[tuple[frozenset[Any], ...], tuple[int, ...]] | None]:
        """
        Combine child DP tables using one unified k-children recurrence.

        Returns
        -------
        tuple[float, tuple[tuple[frozenset[Any], ...], tuple[int, ...]] | None]
            Best objective value and argopt partition/split; returns
            ``(infeasible_value, None)`` if no feasible assignment exists.
        """
        if self.objective == "max":
            if self._use_numba:
                return self._combine_max_numba(
                    child_tables=child_tables,
                    items=items,
                    budget=budget,
                )
            return self._combine_max(child_tables=child_tables, items=items, budget=budget)
        if self._use_numba:
            return self._combine_min_numba(
                child_tables=child_tables,
                items=items,
                budget=budget,
            )
        return self._combine_min(child_tables=child_tables, items=items, budget=budget)

    def combine_full_table(
        self,
        child_tables: list[dict[frozenset[Any], dict[int, float]]],
        items: frozenset[Any],
        budget: int,
    ) -> tuple[np.ndarray, np.ndarray, int, list[frozenset[Any]]]:
        """Compute the full merge table for all (mask, budget) pairs at once.

        Returns
        -------
        tuple[np.ndarray, np.ndarray, int, list[frozenset[Any]]]
            ``(values, choices, full_mask, subset_by_mask)`` where:
            - ``values`` has shape ``[full_mask+1, budget+1]``
            - ``choices`` has shape ``[deg, full_mask+1, budget+1, 4]``
        """
        infeasible = self.infeasible_value
        full_mask, subset_by_mask = self._subset_mask_data(items)
        max_mask = full_mask + 1
        deg = len(child_tables)

        if deg == 0:
            values = np.full((max_mask, budget + 1), infeasible, dtype=np.float64)
            values[0, 0] = 0.0
            choices = np.full((0, max_mask, budget + 1, 4), -1, dtype=np.int64)
            return values, choices, full_mask, subset_by_mask

        if deg == 1:
            mask_by_subset = {s: m for m, s in enumerate(subset_by_mask)}
            values = np.full((max_mask, budget + 1), infeasible, dtype=np.float64)
            choices = np.full((1, max_mask, budget + 1, 4), -1, dtype=np.int64)
            for subset, budget_map in child_tables[0].items():
                mask = mask_by_subset.get(subset)
                if mask is None:
                    continue
                for b_val, val in budget_map.items():
                    if 0 <= b_val <= budget and val != infeasible:
                        values[mask, b_val] = val
                        choices[0, mask, b_val, 0] = 0
                        choices[0, mask, b_val, 1] = 0
                        choices[0, mask, b_val, 2] = mask
                        choices[0, mask, b_val, 3] = b_val
            return values, choices, full_mask, subset_by_mask

        child_values = self._dense_child_values(child_tables, subset_by_mask, budget)
        if self._use_numba:
            if self.objective == "max":
                values, choices = _combine_core_max_numba(
                    child_values, infeasible, budget, full_mask
                )
            else:
                values, choices = _combine_core_min_numba(
                    child_values, infeasible, budget, full_mask
                )
        else:
            values, choices = self._combine_full_python(
                child_tables, full_mask, subset_by_mask, budget
            )
        return values, choices, full_mask, subset_by_mask

    def _dense_child_values(
        self,
        child_tables: list[dict[frozenset[Any], dict[int, float]]],
        subset_by_mask: list[frozenset[Any]],
        budget: int,
    ) -> np.ndarray:
        """Build dense child value tensor ``[child, subset_mask, budget]``."""
        infeasible = self.infeasible_value
        max_mask = len(subset_by_mask)
        child_values = np.full(
            (len(child_tables), max_mask, budget + 1),
            infeasible,
            dtype=np.float64,
        )
        mask_by_subset = {subset: mask for mask, subset in enumerate(subset_by_mask)}
        for child_idx, child_table in enumerate(child_tables):
            for subset, budget_map in child_table.items():
                mask = mask_by_subset.get(subset)
                if mask is None:
                    continue
                for allocated_budget, value in budget_map.items():
                    if 0 <= allocated_budget <= budget:
                        child_values[child_idx, mask, allocated_budget] = value
        return child_values

    def _combine_max_numba(
        self,
        child_tables: list[dict[frozenset[Any], dict[int, float]]],
        items: frozenset[Any],
        budget: int,
    ) -> tuple[float, tuple[tuple[frozenset[Any], ...], tuple[int, ...]] | None]:
        """Numba-accelerated maximization backend."""
        deg = len(child_tables)
        if deg == 0:
            if len(items) == 0 and budget == 0:
                return 0.0, (tuple(), tuple())
            return self.infeasible_value, None
        if deg == 1:
            val = child_tables[0].get(items, {}).get(budget, self.infeasible_value)
            if val == self.infeasible_value:
                return self.infeasible_value, None
            return val, ((items,), (budget,))

        full_mask, subset_by_mask = self._subset_mask_data(items)
        child_values = self._dense_child_values(child_tables, subset_by_mask, budget)
        values, choices = _combine_core_max_numba(
            child_values=child_values,
            infeasible=self.infeasible_value,
            budget=budget,
            full_mask=full_mask,
        )
        value = float(values[full_mask, budget])
        if value == self.infeasible_value:
            return self.infeasible_value, None

        part_masks_rev: list[int] = []
        budget_parts_rev: list[int] = []
        mask = full_mask
        used_budget = budget
        for j in range(deg - 1, -1, -1):
            prev_mask = int(choices[j, mask, used_budget, 0])
            prev_budget = int(choices[j, mask, used_budget, 1])
            picked_submask = int(choices[j, mask, used_budget, 2])
            picked_budget = int(choices[j, mask, used_budget, 3])
            if prev_mask < 0:
                return self.infeasible_value, None
            part_masks_rev.append(picked_submask)
            budget_parts_rev.append(picked_budget)
            mask, used_budget = prev_mask, prev_budget

        part_masks_rev.reverse()
        budget_parts_rev.reverse()
        item_parts = tuple(subset_by_mask[mask_value] for mask_value in part_masks_rev)
        budget_parts = tuple(budget_parts_rev)
        return value, (item_parts, budget_parts)

    def _combine_min_numba(
        self,
        child_tables: list[dict[frozenset[Any], dict[int, float]]],
        items: frozenset[Any],
        budget: int,
    ) -> tuple[float, tuple[tuple[frozenset[Any], ...], tuple[int, ...]] | None]:
        """Numba-accelerated minimization backend."""
        deg = len(child_tables)
        if deg == 0:
            if len(items) == 0 and budget == 0:
                return 0.0, (tuple(), tuple())
            return self.infeasible_value, None
        if deg == 1:
            val = child_tables[0].get(items, {}).get(budget, self.infeasible_value)
            if val == self.infeasible_value:
                return self.infeasible_value, None
            return val, ((items,), (budget,))

        full_mask, subset_by_mask = self._subset_mask_data(items)
        child_values = self._dense_child_values(child_tables, subset_by_mask, budget)
        values, choices = _combine_core_min_numba(
            child_values=child_values,
            infeasible=self.infeasible_value,
            budget=budget,
            full_mask=full_mask,
        )
        value = float(values[full_mask, budget])
        if value == self.infeasible_value:
            return self.infeasible_value, None

        part_masks_rev: list[int] = []
        budget_parts_rev: list[int] = []
        mask = full_mask
        used_budget = budget
        for j in range(deg - 1, -1, -1):
            prev_mask = int(choices[j, mask, used_budget, 0])
            prev_budget = int(choices[j, mask, used_budget, 1])
            picked_submask = int(choices[j, mask, used_budget, 2])
            picked_budget = int(choices[j, mask, used_budget, 3])
            if prev_mask < 0:
                return self.infeasible_value, None
            part_masks_rev.append(picked_submask)
            budget_parts_rev.append(picked_budget)
            mask, used_budget = prev_mask, prev_budget

        part_masks_rev.reverse()
        budget_parts_rev.reverse()
        item_parts = tuple(subset_by_mask[mask_value] for mask_value in part_masks_rev)
        budget_parts = tuple(budget_parts_rev)
        return value, (item_parts, budget_parts)

    def _combine_max(
        self,
        child_tables: list[dict[frozenset[Any], dict[int, float]]],
        items: frozenset[Any],
        budget: int,
    ) -> tuple[float, tuple[tuple[frozenset[Any], ...], tuple[int, ...]] | None]:
        """Combine child tables in maximization mode."""
        deg = len(child_tables)
        if deg == 0:
            if len(items) == 0 and budget == 0:
                return 0.0, (tuple(), tuple())
            return self.infeasible_value, None
        if deg == 1:
            val = child_tables[0].get(items, {}).get(budget, self.infeasible_value)
            if val == self.infeasible_value:
                return self.infeasible_value, None
            return val, ((items,), (budget,))

        full_mask, subset_by_mask = self._subset_mask_data(items)
        infeasible = self.infeasible_value

        max_mask = full_mask + 1
        values: list[list[float]] = [[infeasible] * (budget + 1) for _ in range(max_mask)]
        values[0][0] = 0.0
        active_states: list[tuple[int, int]] = [(0, 0)]
        layer_choices: list[list[list[tuple[int, int, int, int] | None]]] = []

        for child_table in child_tables:
            next_values: list[list[float]] = [[infeasible] * (budget + 1) for _ in range(max_mask)]
            next_choice: list[list[tuple[int, int, int, int] | None]] = [
                [None] * (budget + 1) for _ in range(max_mask)
            ]
            next_active_states: list[tuple[int, int]] = []
            next_seen: list[list[bool]] = [[False] * (budget + 1) for _ in range(max_mask)]
            for assigned_mask, used_budget in active_states:
                base_value = values[assigned_mask][used_budget]
                remaining_mask = full_mask ^ assigned_mask
                remaining_budget = budget - used_budget
                submask = remaining_mask
                while True:
                    child_budget_map = child_table.get(subset_by_mask[submask])
                    if child_budget_map is None:
                        if submask == 0:
                            break
                        submask = (submask - 1) & remaining_mask
                        continue
                    new_mask = assigned_mask | submask
                    for allocated_budget, child_value in child_budget_map.items():
                        if allocated_budget > remaining_budget or child_value == infeasible:
                            continue
                        new_budget = used_budget + allocated_budget
                        candidate_value = base_value + child_value
                        current = next_values[new_mask][new_budget]
                        if candidate_value > current:
                            next_values[new_mask][new_budget] = candidate_value
                            next_choice[new_mask][new_budget] = (
                                assigned_mask,
                                used_budget,
                                submask,
                                allocated_budget,
                            )
                            if not next_seen[new_mask][new_budget]:
                                next_seen[new_mask][new_budget] = True
                                next_active_states.append((new_mask, new_budget))
                    if submask == 0:
                        break
                    submask = (submask - 1) & remaining_mask
            values = next_values
            active_states = next_active_states
            layer_choices.append(next_choice)
            if not active_states:
                return infeasible, None

        value = values[full_mask][budget]
        if value == infeasible:
            return infeasible, None

        part_masks_rev: list[int] = []
        budget_parts_rev: list[int] = []
        key = (full_mask, budget)
        for j in range(deg - 1, -1, -1):
            prev = layer_choices[j][key[0]][key[1]]
            if prev is None:
                return infeasible, None
            prev_mask, prev_budget, picked_submask, picked_budget = prev
            part_masks_rev.append(picked_submask)
            budget_parts_rev.append(picked_budget)
            key = (prev_mask, prev_budget)

        part_masks_rev.reverse()
        budget_parts_rev.reverse()
        item_parts = tuple(subset_by_mask[mask] for mask in part_masks_rev)
        budget_parts = tuple(budget_parts_rev)
        return value, (item_parts, budget_parts)

    def _combine_full_python(
        self,
        child_tables: list[dict[frozenset[Any], dict[int, float]]],
        full_mask: int,
        subset_by_mask: list[frozenset[Any]],
        budget: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Pure-Python full-table merge returning dense numpy arrays."""
        infeasible = self.infeasible_value
        deg = len(child_tables)
        max_mask = full_mask + 1
        compare = (lambda a, b: a > b) if self.objective == "max" else (lambda a, b: a < b)

        values: list[list[float]] = [[infeasible] * (budget + 1) for _ in range(max_mask)]
        values[0][0] = 0.0
        active_states: list[tuple[int, int]] = [(0, 0)]
        layer_choices: list[list[list[tuple[int, int, int, int] | None]]] = []

        for child_table in child_tables:
            next_values: list[list[float]] = [[infeasible] * (budget + 1) for _ in range(max_mask)]
            next_choice: list[list[tuple[int, int, int, int] | None]] = [
                [None] * (budget + 1) for _ in range(max_mask)
            ]
            next_active_states: list[tuple[int, int]] = []
            next_seen: list[list[bool]] = [[False] * (budget + 1) for _ in range(max_mask)]
            for assigned_mask, used_budget in active_states:
                base_value = values[assigned_mask][used_budget]
                remaining_mask = full_mask ^ assigned_mask
                remaining_budget = budget - used_budget
                submask = remaining_mask
                while True:
                    child_budget_map = child_table.get(subset_by_mask[submask])
                    if child_budget_map is None:
                        if submask == 0:
                            break
                        submask = (submask - 1) & remaining_mask
                        continue
                    new_mask = assigned_mask | submask
                    for allocated_budget, child_value in child_budget_map.items():
                        if allocated_budget > remaining_budget or child_value == infeasible:
                            continue
                        new_budget = used_budget + allocated_budget
                        candidate_value = base_value + child_value
                        current = next_values[new_mask][new_budget]
                        if compare(candidate_value, current):
                            next_values[new_mask][new_budget] = candidate_value
                            next_choice[new_mask][new_budget] = (
                                assigned_mask,
                                used_budget,
                                submask,
                                allocated_budget,
                            )
                            if not next_seen[new_mask][new_budget]:
                                next_seen[new_mask][new_budget] = True
                                next_active_states.append((new_mask, new_budget))
                    if submask == 0:
                        break
                    submask = (submask - 1) & remaining_mask
            values = next_values
            active_states = next_active_states
            layer_choices.append(next_choice)

        values_arr = np.array(values, dtype=np.float64)
        choices_arr = np.full((deg, max_mask, budget + 1, 4), -1, dtype=np.int64)
        for j, layer in enumerate(layer_choices):
            for m in range(max_mask):
                for bb in range(budget + 1):
                    entry = layer[m][bb]
                    if entry is not None:
                        choices_arr[j, m, bb, 0] = entry[0]
                        choices_arr[j, m, bb, 1] = entry[1]
                        choices_arr[j, m, bb, 2] = entry[2]
                        choices_arr[j, m, bb, 3] = entry[3]
        return values_arr, choices_arr

    def _combine_min(
        self,
        child_tables: list[dict[frozenset[Any], dict[int, float]]],
        items: frozenset[Any],
        budget: int,
    ) -> tuple[float, tuple[tuple[frozenset[Any], ...], tuple[int, ...]] | None]:
        """Combine child tables in minimization mode."""
        deg = len(child_tables)
        if deg == 0:
            if len(items) == 0 and budget == 0:
                return 0.0, (tuple(), tuple())
            return self.infeasible_value, None
        if deg == 1:
            val = child_tables[0].get(items, {}).get(budget, self.infeasible_value)
            if val == self.infeasible_value:
                return self.infeasible_value, None
            return val, ((items,), (budget,))

        full_mask, subset_by_mask = self._subset_mask_data(items)
        infeasible = self.infeasible_value

        max_mask = full_mask + 1
        values: list[list[float]] = [[infeasible] * (budget + 1) for _ in range(max_mask)]
        values[0][0] = 0.0
        active_states: list[tuple[int, int]] = [(0, 0)]
        layer_choices: list[list[list[tuple[int, int, int, int] | None]]] = []

        for child_table in child_tables:
            next_values: list[list[float]] = [[infeasible] * (budget + 1) for _ in range(max_mask)]
            next_choice: list[list[tuple[int, int, int, int] | None]] = [
                [None] * (budget + 1) for _ in range(max_mask)
            ]
            next_active_states: list[tuple[int, int]] = []
            next_seen: list[list[bool]] = [[False] * (budget + 1) for _ in range(max_mask)]
            for assigned_mask, used_budget in active_states:
                base_value = values[assigned_mask][used_budget]
                remaining_mask = full_mask ^ assigned_mask
                remaining_budget = budget - used_budget
                submask = remaining_mask
                while True:
                    child_budget_map = child_table.get(subset_by_mask[submask])
                    if child_budget_map is None:
                        if submask == 0:
                            break
                        submask = (submask - 1) & remaining_mask
                        continue
                    new_mask = assigned_mask | submask
                    for allocated_budget, child_value in child_budget_map.items():
                        if allocated_budget > remaining_budget or child_value == infeasible:
                            continue
                        new_budget = used_budget + allocated_budget
                        candidate_value = base_value + child_value
                        current = next_values[new_mask][new_budget]
                        if candidate_value < current:
                            next_values[new_mask][new_budget] = candidate_value
                            next_choice[new_mask][new_budget] = (
                                assigned_mask,
                                used_budget,
                                submask,
                                allocated_budget,
                            )
                            if not next_seen[new_mask][new_budget]:
                                next_seen[new_mask][new_budget] = True
                                next_active_states.append((new_mask, new_budget))
                    if submask == 0:
                        break
                    submask = (submask - 1) & remaining_mask
            values = next_values
            active_states = next_active_states
            layer_choices.append(next_choice)
            if not active_states:
                return infeasible, None

        value = values[full_mask][budget]
        if value == infeasible:
            return infeasible, None

        part_masks_rev: list[int] = []
        budget_parts_rev: list[int] = []
        key = (full_mask, budget)
        for j in range(deg - 1, -1, -1):
            prev = layer_choices[j][key[0]][key[1]]
            if prev is None:
                return infeasible, None
            prev_mask, prev_budget, picked_submask, picked_budget = prev
            part_masks_rev.append(picked_submask)
            budget_parts_rev.append(picked_budget)
            key = (prev_mask, prev_budget)

        part_masks_rev.reverse()
        budget_parts_rev.reverse()
        item_parts = tuple(subset_by_mask[mask] for mask in part_masks_rev)
        budget_parts = tuple(budget_parts_rev)
        return value, (item_parts, budget_parts)
