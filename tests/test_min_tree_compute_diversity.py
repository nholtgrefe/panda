"""Tests for MinTreePD fixed-set diversity computation."""

from __future__ import annotations

from scanwidth import DAG, TreeExtension, node_scanwidth

import phypanda as pp
from tests.baselines import SMALL_TREE_FIXED_DIVERSITY
from tests.test_networks import build_small_tree_network


def test_min_tree_compute_diversity_matches_baseline() -> None:
    """MinTreePD matches stored baseline and AllPaths on the small tree."""
    network = build_small_tree_network()
    taxa = {"a", "c", "d"}
    expected = SMALL_TREE_FIXED_DIVERSITY[frozenset(taxa)]
    min_tree_value = pp.min_tree.compute_diversity(network, taxa)
    all_paths_value = pp.all_paths.compute_diversity(network, taxa)
    assert min_tree_value == expected
    assert min_tree_value == all_paths_value


def test_min_tree_compute_diversity_with_explicit_tree_extension() -> None:
    """Accept explicit tree extension and match default-extension result."""
    network = build_small_tree_network()
    taxa = {"a", "b", "c", "d"}
    expected = SMALL_TREE_FIXED_DIVERSITY[frozenset(taxa)]

    dag = DAG(network._graph._graph)
    _, extension = node_scanwidth(dag)
    tree_extension: TreeExtension = extension.to_canonical_tree_extension()

    default_value = pp.min_tree.compute_diversity(network, taxa)
    explicit_value = pp.min_tree.compute_diversity(
        network,
        taxa,
        tree_extension=tree_extension,
    )
    assert default_value == expected
    assert explicit_value == default_value
