"""Regression tests: edges with length 0 must count as length 0, not as length 1.

Up to version 2.1.1, the solvers read branch lengths with
``get_branch_length(u, v) or 1.0``, which also turned an explicit length of 0
into 1 (only a missing length should default to 1).
"""

from __future__ import annotations

import pytest
from phylozoo import DirectedPhyNetwork

import phypanda as pp
from phypanda.utils import branch_length

# Root r with children p1 and p2; p1 -> a, p2 -> c, and both p1 and p2 are
# parents of the reticulation H (two edges of length 0), whose child is b.
#   all-paths:  {a} = 2, {c} = 2, {b} = 1 + 0 + 1 + 0 + 1 = 3
#   max-tree:   {b} = 1 + 0 + 1 = 2, so the best single taxon has value 2
#   min-tree of all taxa: every displayed tree drops one 0-length edge, total 5
ZERO_LENGTH_NEWICK = "((a:1,(b:1)#H1:0):1,(#H1:0,c:1):1);"


@pytest.fixture
def network() -> DirectedPhyNetwork:
    """Network with two reticulation edges of length 0."""
    return DirectedPhyNetwork.from_string(ZERO_LENGTH_NEWICK)


def test_branch_length_keeps_zero(network: DirectedPhyNetwork) -> None:
    """A length of 0 stays 0, and every edge of the test network has a length."""
    lengths = [branch_length(network, u, v) for u, v, _ in network._graph.edges(keys=True)]
    assert sorted(lengths) == [0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0]


def test_all_paths_counts_zero_length_edges_as_zero(network: DirectedPhyNetwork) -> None:
    """All-paths maximization (both solvers) agrees with ``compute_diversity``."""
    for algorithm in ("nsw_fpt_budget", "esw_fpt"):
        value, taxa = pp.solve_max_diversity(
            network, budget=1, measure=pp.all_paths, algorithm=algorithm
        )
        assert value == pytest.approx(3.0)
        assert taxa == {"b"}
        assert pp.all_paths.compute_diversity(network, taxa) == pytest.approx(value)


def test_max_tree_counts_zero_length_edges_as_zero(network: DirectedPhyNetwork) -> None:
    """MaxTreePD of the best single taxon is 2."""
    value, _ = pp.solve_max_diversity(network, budget=1, measure=pp.max_tree)
    assert value == pytest.approx(2.0)


def test_min_tree_counts_zero_length_edges_as_zero(network: DirectedPhyNetwork) -> None:
    """MinTreePD of all taxa is the total length 5."""
    value = pp.min_tree.compute_diversity(network, set(network.taxa))
    assert value == pytest.approx(5.0)
