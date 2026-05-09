"""Tests for child-merge DP and shared utility helpers."""

from __future__ import annotations

import pytest

from phypanda.utils import _ChildMergeDP, powerset


def _bruteforce_ordered_partitions(items: frozenset[str], parts: int) -> list[tuple[frozenset[str], ...]]:
    """Return all ordered set partitions with empty parts allowed."""
    if parts == 0:
        return [tuple()] if len(items) == 0 else []
    elements = sorted(items)
    buckets: list[set[str]] = [set() for _ in range(parts)]
    out: list[tuple[frozenset[str], ...]] = []

    def recurse(index: int) -> None:
        if index == len(elements):
            out.append(tuple(frozenset(bucket) for bucket in buckets))
            return
        element = elements[index]
        for bucket in buckets:
            bucket.add(element)
            recurse(index + 1)
            bucket.remove(element)

    recurse(0)
    return out


def _bruteforce_splits(total: int, parts: int) -> list[tuple[int, ...]]:
    """Return all ordered non-negative integer splits of ``total``."""
    if parts == 0:
        return [tuple()] if total == 0 else []
    out: list[tuple[int, ...]] = []

    def recurse(index: int, remaining: int, current: list[int]) -> None:
        if index == parts - 1:
            current.append(remaining)
            out.append(tuple(current))
            current.pop()
            return
        for value in range(remaining + 1):
            current.append(value)
            recurse(index + 1, remaining - value, current)
            current.pop()

    recurse(0, total, [])
    return out


def test_child_merge_dp_matches_bruteforce() -> None:
    """Child-merge DP matches exhaustive partition-and-split search."""
    children = ("c1", "c2", "c3")
    items = frozenset({"a", "b"})
    budget = 2
    minus_infinity = -10_000.0
    engine = _ChildMergeDP(infeasible_value=minus_infinity, objective="max")

    score_map: dict[tuple[str, frozenset[str], int], float] = {}
    for child in children:
        for subset_set in powerset(set(items)):
            subset = frozenset(subset_set)
            for budget_value in range(budget + 1):
                score_map[(child, subset, budget_value)] = float(
                    (len(subset) * 10)
                    + (2 * budget_value)
                    + (1 if child == "c1" else 0)
                )

    child_tables: list[dict[frozenset[str], dict[int, float]]] = []
    for child in children:
        table: dict[frozenset[str], dict[int, float]] = {}
        for subset_set in powerset(set(items)):
            subset = frozenset(subset_set)
            table[subset] = {}
            for budget_value in range(budget + 1):
                table[subset][budget_value] = score_map[(child, subset, budget_value)]
        child_tables.append(table)

    value, choice = engine.combine(child_tables=child_tables, items=items, budget=budget)

    best = minus_infinity
    for partition in _bruteforce_ordered_partitions(items, len(children)):
        for split in _bruteforce_splits(budget, len(children)):
            candidate = sum(
                score_map[(child, subset, budget_value)]
                for child, subset, budget_value in zip(children, partition, split)
            )
            if candidate > best:
                best = candidate
    assert value == best
    assert choice is not None
    part_choice, split_choice = choice
    assert len(part_choice) == len(children)
    assert len(split_choice) == len(children)
    assert frozenset().union(*part_choice) == items
    assert sum(split_choice) == budget


def test_child_merge_dp_handles_infeasible_states() -> None:
    """Child-merge DP returns ``minus_infinity`` when no assignment is feasible."""
    minus_infinity = -1234.0
    engine = _ChildMergeDP(infeasible_value=minus_infinity, objective="max")
    child_tables: list[dict[frozenset[str], dict[int, float]]] = [
        {frozenset(): {0: minus_infinity}, frozenset({"x"}): {0: minus_infinity, 1: minus_infinity}},
        {frozenset(): {0: minus_infinity}, frozenset({"x"}): {0: minus_infinity, 1: minus_infinity}},
    ]
    value, choice = engine.combine(
        child_tables=child_tables,
        items=frozenset({"x"}),
        budget=1,
    )
    assert value == minus_infinity
    assert choice is None


def test_child_merge_dp_zero_children_empty_state() -> None:
    """Zero children with empty/zero target is feasible with value 0."""
    engine = _ChildMergeDP(infeasible_value=-9999.0, objective="max")
    value, choice = engine.combine(child_tables=[], items=frozenset(), budget=0)
    assert value == 0.0
    assert choice == (tuple(), tuple())


def test_child_merge_dp_zero_children_non_empty_is_infeasible() -> None:
    """Zero children cannot realize non-empty item or non-zero budget targets."""
    engine = _ChildMergeDP(infeasible_value=-9999.0, objective="max")
    value1, choice1 = engine.combine(child_tables=[], items=frozenset({"x"}), budget=0)
    value2, choice2 = engine.combine(child_tables=[], items=frozenset(), budget=1)
    assert value1 == -9999.0 and choice1 is None
    assert value2 == -9999.0 and choice2 is None


def test_child_merge_dp_single_child_direct_lookup() -> None:
    """Single-child case returns direct table lookup."""
    infeasible = -5000.0
    engine = _ChildMergeDP(infeasible_value=infeasible, objective="max")
    table = {frozenset({"a"}): {1: 7.0}}
    value, choice = engine.combine(child_tables=[table], items=frozenset({"a"}), budget=1)
    assert value == 7.0
    assert choice == ((frozenset({"a"}),), (1,))


def test_child_merge_dp_minimize_matches_bruteforce() -> None:
    """Minimization mode matches exhaustive search."""
    children = ("c1", "c2", "c3")
    items = frozenset({"a", "b"})
    budget = 2
    plus_infinity = 10_000.0
    engine = _ChildMergeDP(infeasible_value=plus_infinity, objective="min")

    score_map: dict[tuple[str, frozenset[str], int], float] = {}
    for child in children:
        for subset_set in powerset(set(items)):
            subset = frozenset(subset_set)
            for budget_value in range(budget + 1):
                score_map[(child, subset, budget_value)] = float(
                    (5 * len(subset))
                    + budget_value
                    + (2 if child == "c2" else 0)
                )

    child_tables: list[dict[frozenset[str], dict[int, float]]] = []
    for child in children:
        table: dict[frozenset[str], dict[int, float]] = {}
        for subset_set in powerset(set(items)):
            subset = frozenset(subset_set)
            table[subset] = {}
            for budget_value in range(budget + 1):
                table[subset][budget_value] = score_map[(child, subset, budget_value)]
        child_tables.append(table)

    value, choice = engine.combine(child_tables=child_tables, items=items, budget=budget)

    best = plus_infinity
    for partition in _bruteforce_ordered_partitions(items, len(children)):
        for split in _bruteforce_splits(budget, len(children)):
            candidate = sum(
                score_map[(child, subset, budget_value)]
                for child, subset, budget_value in zip(children, partition, split)
            )
            if candidate < best:
                best = candidate
    assert value == best
    assert choice is not None
    part_choice, split_choice = choice
    assert len(part_choice) == len(children)
    assert len(split_choice) == len(children)
    assert frozenset().union(*part_choice) == items
    assert sum(split_choice) == budget


def test_child_merge_dp_invalid_objective_raises() -> None:
    """Invalid objective mode raises a ValueError."""
    with pytest.raises(ValueError):
        _ChildMergeDP(infeasible_value=-1.0, objective="median")  # type: ignore[arg-type]


def test_child_merge_dp_numba_false_matches_numba_true() -> None:
    """Pure Python merge matches Numba merge on the same structured instance."""
    children = ("c1", "c2", "c3")
    items = frozenset({"a", "b"})
    budget = 2
    minus_infinity = -10_000.0

    score_map: dict[tuple[str, frozenset[str], int], float] = {}
    for child in children:
        for subset_set in powerset(set(items)):
            subset = frozenset(subset_set)
            for budget_value in range(budget + 1):
                score_map[(child, subset, budget_value)] = float(
                    (len(subset) * 10)
                    + (2 * budget_value)
                    + (1 if child == "c1" else 0)
                )

    child_tables: list[dict[frozenset[str], dict[int, float]]] = []
    for child in children:
        table: dict[frozenset[str], dict[int, float]] = {}
        for subset_set in powerset(set(items)):
            subset = frozenset(subset_set)
            table[subset] = {}
            for budget_value in range(budget + 1):
                table[subset][budget_value] = score_map[(child, subset, budget_value)]
        child_tables.append(table)

    fast = _ChildMergeDP(
        infeasible_value=minus_infinity,
        objective="max",
        use_numba=True,
    )
    slow = _ChildMergeDP(
        infeasible_value=minus_infinity,
        objective="max",
        use_numba=False,
    )
    v1, ch1 = fast.combine(child_tables=child_tables, items=items, budget=budget)
    v2, ch2 = slow.combine(child_tables=child_tables, items=items, budget=budget)
    assert v1 == v2
    assert ch1 is not None and ch2 is not None

    best = minus_infinity
    for partition in _bruteforce_ordered_partitions(items, len(children)):
        for split in _bruteforce_splits(budget, len(children)):
            candidate = sum(
                score_map[(child, subset, budget_value)]
                for child, subset, budget_value in zip(children, partition, split)
            )
            if candidate > best:
                best = candidate
    assert v1 == best
    for choice in (ch1, ch2):
        part_choice, split_choice = choice
        assert sum(split_choice) == budget
        assert frozenset().union(*part_choice) == items
        achieved = sum(
            score_map[(child, subset, budget_value)]
            for child, subset, budget_value in zip(children, part_choice, split_choice)
        )
        assert achieved == best
