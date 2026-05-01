"""Regression checks for PD baselines on fixed experiment networks."""

from __future__ import annotations

import pytest

import phypanda as pp

from tests.baselines import (
    APD_BASELINES,
    APD_BUDGETS,
    EXP1_FIXED_SUBSET_DIVERSITY_BASELINES,
    EXP1_MIN_TREE_FULL_LEAF_DIVERSITY,
    MAX_TREE_BUDGET_BASELINES,
)
from tests.test_networks import EXP1_NETWORK_IDS, build_network


@pytest.mark.parametrize("network_id", EXP1_NETWORK_IDS)
def test_selected_networks_have_20_leaves(network_id: str) -> None:
    """Validate that each exp1 benchmark network has 20 leaves.

    Parameters
    ----------
    network_id : str
        Identifier of the selected benchmark network.
    """
    network = build_network(network_id)
    assert len(network.leaves) == 20


@pytest.mark.parametrize("network_id", sorted(APD_BASELINES))
@pytest.mark.parametrize("budget", APD_BUDGETS)
def test_all_paths_baselines_match(network_id: str, budget: int) -> None:
    """Check computed APD values against stored baselines.

    Parameters
    ----------
    network_id : str
        Identifier of the selected benchmark network.
    budget : int
        Taxon budget used for the maximization run.
    """
    expected = APD_BASELINES[network_id][budget]
    network = build_network(network_id)
    value, solution = pp.solve_max_diversity(network, budget=budget, measure=pp.all_paths)

    assert value == pytest.approx(expected["value"])
    assert len(solution) == budget


@pytest.mark.parametrize("network_id", sorted(MAX_TREE_BUDGET_BASELINES))
@pytest.mark.parametrize("budget", APD_BUDGETS)
def test_max_tree_budget_baselines_match(network_id: str, budget: int) -> None:
    """MaxTreePD unit-cost maximization matches stored baselines on exp1 networks."""
    expected_value = MAX_TREE_BUDGET_BASELINES[network_id][budget]
    network = build_network(network_id)
    value, solution = pp.solve_max_diversity(network, budget=budget, measure=pp.max_tree)

    assert value == pytest.approx(expected_value)
    assert len(solution) == budget


@pytest.mark.parametrize("network_id", sorted(EXP1_FIXED_SUBSET_DIVERSITY_BASELINES))
@pytest.mark.parametrize("budget", APD_BUDGETS)
def test_fixed_subset_all_paths_and_min_tree_baselines(
    network_id: str,
    budget: int,
) -> None:
    """Fixed-set all-paths and MinTreePD values match baselines (explicit taxa per row)."""
    row = EXP1_FIXED_SUBSET_DIVERSITY_BASELINES[network_id][budget]
    taxa = set(row["taxa"])
    network = build_network(network_id)

    ap = pp.all_paths.compute_diversity(network, taxa)
    mt = pp.min_tree.compute_diversity(network, taxa)

    assert ap == pytest.approx(row["all_paths"])
    assert mt == pytest.approx(row["min_tree"])


@pytest.mark.parametrize("network_id", sorted(EXP1_MIN_TREE_FULL_LEAF_DIVERSITY))
def test_min_tree_full_leaf_baselines(network_id: str) -> None:
    """MinTreePD over all leaves matches stored baseline."""
    network = build_network(network_id)
    value = pp.min_tree.compute_diversity(network, set(network.taxa))
    assert value == pytest.approx(EXP1_MIN_TREE_FULL_LEAF_DIVERSITY[network_id])
