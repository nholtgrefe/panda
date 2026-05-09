"""Tests for the node-scanwidth budgeted all-paths solver."""

from __future__ import annotations

import pytest
from scanwidth import DAG, TreeExtension, node_scanwidth
from phylozoo.utils.exceptions import PhyloZooValueError

import phypanda as pp
from tests.baselines import APD_BASELINES, APD_BUDGETS, SMALL_TREE_BUDGET_COSTS
from tests.test_networks import build_network, build_small_tree_network


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
    network = build_small_tree_network()
    costs = SMALL_TREE_BUDGET_COSTS
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


def test_nsw_fpt_budget_rejects_negative_costs() -> None:
    """Negative taxon costs are rejected."""
    network = build_small_tree_network()
    with pytest.raises(PhyloZooValueError, match="non-negative"):
        pp.solve_max_diversity(
            network,
            budget=2,
            costs={"a": -1, "b": 1, "c": 1, "d": 1},
            measure=pp.all_paths,
            algorithm="nsw_fpt_budget",
        )


def test_nsw_fpt_budget_defaults_missing_costs_to_one() -> None:
    """Missing taxon costs default to unit cost."""
    network = build_small_tree_network()
    budget = 3
    partial_costs = {"a": 2}
    explicit_costs = {"a": 2, "b": 1, "c": 1, "d": 1}

    value_partial, taxa_partial = pp.solve_max_diversity(
        network,
        budget=budget,
        costs=partial_costs,
        measure=pp.all_paths,
        algorithm="nsw_fpt_budget",
    )
    value_explicit, taxa_explicit = pp.solve_max_diversity(
        network,
        budget=budget,
        costs=explicit_costs,
        measure=pp.all_paths,
        algorithm="nsw_fpt_budget",
    )
    assert value_partial == value_explicit
    assert taxa_partial == taxa_explicit
