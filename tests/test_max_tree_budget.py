"""Tests for MaxTreePD budgeted DP solver."""

from __future__ import annotations

import pytest
from scanwidth import DAG, TreeExtension, node_scanwidth
from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooValueError

import phypanda as pp
from tests.networks_exp1 import build_network


def _small_tree_network() -> DirectedPhyNetwork:
    """Create a small weighted rooted tree with four taxa."""
    edges = [
        {"u": "r", "v": "u", "branch_length": 1.0},
        {"u": "r", "v": "v", "branch_length": 1.0},
        {"u": "u", "v": "a", "branch_length": 4.0},
        {"u": "u", "v": "b", "branch_length": 2.0},
        {"u": "v", "v": "c", "branch_length": 3.0},
        {"u": "v", "v": "d", "branch_length": 1.5},
    ]
    nodes = [
        ("a", {"label": "a"}),
        ("b", {"label": "b"}),
        ("c", {"label": "c"}),
        ("d", {"label": "d"}),
    ]
    return DirectedPhyNetwork(edges=edges, nodes=nodes)


def test_max_tree_matches_all_paths_on_tree_instance() -> None:
    """On trees, MaxTreePD objective matches AllPaths objective."""
    network = _small_tree_network()
    costs = {"a": 2, "b": 1, "c": 2, "d": 1}
    budget = 3

    max_tree_value, max_tree_taxa = pp.solve_max_diversity(
        network,
        budget=budget,
        costs=costs,
        measure=pp.max_tree,
    )
    all_paths_value, all_paths_taxa = pp.solve_max_diversity(
        network,
        budget=budget,
        costs=costs,
        measure=pp.all_paths,
        algorithm="nsw_fpt_budget",
    )

    assert max_tree_value == all_paths_value
    assert max_tree_taxa == all_paths_taxa


def test_max_tree_accepts_explicit_tree_extension() -> None:
    """MaxTree solver accepts explicit node-scanwidth tree extension."""
    network = build_network("00001")
    budget = 2
    dag = DAG(network._graph._graph)
    _, extension = node_scanwidth(dag)
    tree_extension: TreeExtension = extension.to_canonical_tree_extension()

    value_default, taxa_default = pp.solve_max_diversity(
        network,
        budget=budget,
        measure=pp.max_tree,
    )
    value_explicit, taxa_explicit = pp.solve_max_diversity(
        network,
        budget=budget,
        measure=pp.max_tree,
        tree_extension=tree_extension,
    )

    assert value_explicit == value_default
    assert taxa_explicit == taxa_default


def test_max_tree_rejects_zero_costs() -> None:
    """Zero taxon costs are rejected (costs must be strictly positive)."""
    network = _small_tree_network()
    with pytest.raises(PhyloZooValueError, match="positive integer|positive"):
        pp.solve_max_diversity(
            network,
            budget=2,
            costs={"a": 0, "b": 1, "c": 1, "d": 1},
            measure=pp.max_tree,
        )


def test_max_tree_defaults_missing_costs_to_one() -> None:
    """Missing taxon costs default to unit cost."""
    network = _small_tree_network()
    budget = 3
    partial_costs = {"a": 2}
    explicit_costs = {"a": 2, "b": 1, "c": 1, "d": 1}

    value_partial, taxa_partial = pp.solve_max_diversity(
        network,
        budget=budget,
        costs=partial_costs,
        measure=pp.max_tree,
    )
    value_explicit, taxa_explicit = pp.solve_max_diversity(
        network,
        budget=budget,
        costs=explicit_costs,
        measure=pp.max_tree,
    )
    assert value_partial == value_explicit
    assert taxa_partial == taxa_explicit
