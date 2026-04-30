"""
General utility helpers for the phypanda package.
"""

from __future__ import annotations

from itertools import combinations
from typing import Any, Iterator, Literal

import numpy as np

try:
    from numba import njit
except Exception:  # pragma: no cover - fallback when numba unavailable
    njit = None


if njit is not None:

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


class _ChildMergeDP:
    """
    Reusable DP engine for partition-and-budget child-table merges.

    Parameters
    ----------
    infeasible_value : float
        Sentinel value that denotes infeasible states.
    objective : {"max", "min"}, default="max"
        Optimization direction for the merge recurrence.
    """

    def __init__(
        self,
        infeasible_value: float,
        objective: Literal["max", "min"] = "max",
    ) -> None:
        self.infeasible_value = infeasible_value
        self.objective = objective
        self._subset_mask_cache: dict[frozenset[Any], tuple[int, list[frozenset[Any]]]] = {}
        if objective not in {"max", "min"}:
            raise ValueError(f"objective must be 'max' or 'min', got {objective!r}")
        self._use_numba = njit is not None

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
        max_mask = full_mask + 1
        infeasible = self.infeasible_value

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
            next_seen: set[tuple[int, int]] = set()
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
                    for allocated_budget in range(remaining_budget + 1):
                        child_value = child_budget_map.get(allocated_budget, infeasible)
                        if child_value == infeasible:
                            continue
                        new_mask = assigned_mask | submask
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
                            key = (new_mask, new_budget)
                            if key not in next_seen:
                                next_seen.add(key)
                                next_active_states.append(key)
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
        max_mask = full_mask + 1
        infeasible = self.infeasible_value

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
            next_seen: set[tuple[int, int]] = set()
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
                    for allocated_budget in range(remaining_budget + 1):
                        child_value = child_budget_map.get(allocated_budget, infeasible)
                        if child_value == infeasible:
                            continue
                        new_mask = assigned_mask | submask
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
                            key = (new_mask, new_budget)
                            if key not in next_seen:
                                next_seen.add(key)
                                next_active_states.append(key)
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
