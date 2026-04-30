"""Tests for the node-scanwidth budgeted all-paths solver."""

from __future__ import annotations

import pytest
from scanwidth import DAG, node_scanwidth
from scanwidth import TreeExtension
from phylozoo.core.network.dnetwork import DirectedPhyNetwork

import phypanda as pp
from tests.apd_baselines import APD_BASELINES, APD_BUDGETS
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


def test_nsw_fpt_budget_matches_esw_on_unit_costs_exp1_networks() -> None:
    """For unit costs, NSW and ESW should match values on selected exp1 networks."""
    for network_id in sorted(APD_BASELINES):
        for budget in APD_BUDGETS:
            esw_value, _ = pp.solve_max_diversity(
                build_network(network_id),
                budget=budget,
                measure=pp.all_paths,
                algorithm="esw_fpt",
            )
            nsw_value, _ = pp.solve_max_diversity(
                build_network(network_id),
                budget=budget,
                measure=pp.all_paths,
                algorithm="nsw_fpt_budget",
            )
            assert nsw_value == pytest.approx(esw_value)


def test_nsw_fpt_budget_with_explicit_tree_extension() -> None:
    """Accept explicit tree extension and match default-extension result."""
    network = _small_tree_network()
    costs = {"a": 2, "b": 1, "c": 2, "d": 1}
    budget = 3

    dag = DAG(network._graph._graph)
    _, extension = node_scanwidth(dag)
    tree_extension: TreeExtension = extension.to_canonical_tree_extension()

    default_value, default_taxa = pp.solve_max_diversity(
        network,
        budget=budget,
        costs=costs,
        measure=pp.all_paths,
        algorithm="nsw_fpt_budget",
    )
    explicit_value, explicit_taxa = pp.solve_max_diversity(
        network,
        budget=budget,
        costs=costs,
        measure=pp.all_paths,
        algorithm="nsw_fpt_budget",
        tree_extension=tree_extension,
    )

    assert explicit_value == default_value
    assert explicit_taxa == default_taxa

