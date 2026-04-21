"""Regression tests for updated phypanda API."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pytest
from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooNotImplementedError

import phypanda as pp


def _network_1() -> DirectedPhyNetwork:
    """Return a small reticulate network matching legacy benchmark."""
    edges = [
        {"u": "r", "v": "x", "branch_length": 1.0},
        {"u": "r", "v": "y", "branch_length": 1.0},
        {"u": "x", "v": "a", "branch_length": 6.0},
        {"u": "x", "v": "h", "branch_length": 2.0},
        {"u": "y", "v": "h", "branch_length": 2.0},
        {"u": "y", "v": "c", "branch_length": 2.0},
        {"u": "h", "v": "b", "branch_length": 4.0},
    ]
    nodes = [("a", {"label": "a"}), ("b", {"label": "b"}), ("c", {"label": "c"})]
    return DirectedPhyNetwork(edges=edges, nodes=nodes)


def _network_2() -> DirectedPhyNetwork:
    """Return a small tree network."""
    edges = [
        {"u": "r", "v": "x", "branch_length": 3.0},
        {"u": "r", "v": "c", "branch_length": 4.0},
        {"u": "x", "v": "a", "branch_length": 1.0},
        {"u": "x", "v": "b", "branch_length": 2.0},
    ]
    nodes = [("a", {"label": "a"}), ("b", {"label": "b"}), ("c", {"label": "c"})]
    return DirectedPhyNetwork(edges=edges, nodes=nodes)


def test_legacy_baseline_network1_budget1() -> None:
    """Match legacy MAPPD baseline for network 1, budget 1."""
    value, solution = pp.solve_max_diversity(_network_1(), budget=1, measure=pp.all_paths)
    assert value == 10.0
    assert solution == {"b"}


def test_legacy_baseline_network1_budget2() -> None:
    """Match legacy MAPPD baseline for network 1, budget 2."""
    value, solution = pp.solve_max_diversity(_network_1(), budget=2, measure=pp.all_paths)
    assert value == 16.0
    assert solution == {"a", "b"}


def test_legacy_baseline_network2_budget1() -> None:
    """Match legacy MAPPD baseline for network 2, budget 1."""
    value, solution = pp.solve_max_diversity(_network_2(), budget=1, measure=pp.all_paths)
    assert value == 5.0
    assert solution == {"b"}


def test_legacy_baseline_network2_budget2() -> None:
    """Match legacy MAPPD baseline for network 2, budget 2."""
    value, solution = pp.solve_max_diversity(_network_2(), budget=2, measure=pp.all_paths)
    assert value == 9.0
    assert solution == {"b", "c"}


def test_non_unit_costs_not_implemented_for_all_paths() -> None:
    """Reject non-unit costs in exact all-paths solver."""
    with pytest.raises(PhyloZooNotImplementedError):
        pp.solve_max_diversity(
            _network_1(),
            budget=2,
            costs={"a": 2, "b": 1, "c": 1},
            measure=pp.all_paths,
        )
