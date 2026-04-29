"""Regression checks for APD baselines on fixed experiment networks."""

from __future__ import annotations

import pytest

import phypanda as pp

from tests.apd_baselines import APD_BASELINES, APD_BUDGETS
from tests.networks_exp1 import NETWORK_DEFINITIONS, build_network


@pytest.mark.parametrize("network_id", sorted(NETWORK_DEFINITIONS))
def test_selected_networks_have_20_leaves(network_id: str) -> None:
    """Validate that each selected benchmark network has 20 leaves.

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

