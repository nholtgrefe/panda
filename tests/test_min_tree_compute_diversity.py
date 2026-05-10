"""Tests for MinTreePD fixed-set diversity computation."""

from __future__ import annotations

import math

from scanwidth import DAG, TreeExtension, node_scanwidth

import phypanda as pp
from phypanda.measure.min_tree import adapt_tree_extension_for_pruned_network, prune_to_taxa
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


def test_prune_to_taxa_all_leaves_returns_same_network() -> None:
    """Inducing on the full leaf set must return the original network object."""
    network = build_small_tree_network()
    out = prune_to_taxa(network, set(network.taxa))
    assert out is network


def test_adapted_full_network_extension_is_valid_for_pruned_network() -> None:
    """Γ shortcutting yields a scanwidth-valid tree on the pruned vertex set."""
    network = build_small_tree_network()
    taxa = {"a", "c", "d"}
    dag = DAG(network._graph._graph)
    _, extension = node_scanwidth(dag)
    tree_extension: TreeExtension = extension.to_canonical_tree_extension()

    working = prune_to_taxa(network, taxa)
    adapted = adapt_tree_extension_for_pruned_network(working, tree_extension)
    assert set(adapted.tree.nodes()) == set(working._graph.nodes())


def test_min_tree_with_full_network_extension_on_subset_taxa() -> None:
    """Explicit full-network extension is adapted when taxa form a proper subset."""
    network = build_small_tree_network()
    taxa = {"a", "c", "d"}
    dag = DAG(network._graph._graph)
    _, extension = node_scanwidth(dag)
    tree_extension: TreeExtension = extension.to_canonical_tree_extension()

    baseline = pp.min_tree.compute_diversity(network, taxa)
    with_full_te = pp.min_tree.compute_diversity(
        network,
        taxa,
        tree_extension=tree_extension,
    )
    assert math.isfinite(baseline)
    assert math.isfinite(with_full_te)
    assert with_full_te == baseline
