"""Tests for :func:`~phypanda.measure.min_tree.adapt_tree_extension_for_pruned_network`."""

from __future__ import annotations

import math

import networkx as nx
from scanwidth import DAG, TreeExtension, node_scanwidth

import phypanda as pp
from phypanda.measure.min_tree import adapt_tree_extension_for_pruned_network, prune_to_taxa
from tests.baselines import SMALL_TREE_FIXED_DIVERSITY
from tests.test_networks import build_small_tree_network


def _canonical_extension(network) -> TreeExtension:
    dag = DAG(network._graph._graph)
    _, extension = node_scanwidth(dag)
    return extension.to_canonical_tree_extension()


def test_adapt_full_network_is_no_op_on_gamma_structure() -> None:
    """When the working network is the full graph, Γ keeps the same edges."""
    network = build_small_tree_network()
    te = _canonical_extension(network)
    adapted = adapt_tree_extension_for_pruned_network(network, te)
    assert set(adapted.tree.nodes()) == set(te.tree.nodes())
    assert set(adapted.tree.edges()) == set(te.tree.edges())
    assert adapted.node_scanwidth() == te.node_scanwidth()


def test_adapt_full_network_is_isomorphic_copy_not_same_graph_object() -> None:
    """Adaptation returns a new TreeExtension with a copied Γ (not ``te.tree``)."""
    network = build_small_tree_network()
    te = _canonical_extension(network)
    adapted = adapt_tree_extension_for_pruned_network(network, te)
    assert adapted.tree is not te.tree
    assert nx.is_isomorphic(
        te.tree.to_undirected(),
        adapted.tree.to_undirected(),
    )


def test_adapt_pruned_network_vertex_set_and_scanwidth() -> None:
    """Pruned adaptation yields a valid extension on the working vertex set."""
    network = build_small_tree_network()
    te = _canonical_extension(network)
    taxa = {"a", "c", "d"}
    working = prune_to_taxa(network, taxa)
    adapted = adapt_tree_extension_for_pruned_network(working, te)
    assert set(adapted.tree.nodes()) == set(working._graph.nodes())
    nsw = adapted.node_scanwidth()
    assert isinstance(nsw, int)
    assert nsw >= 0
    assert nsw <= te.node_scanwidth()


def test_adapt_pruned_min_tree_matches_baseline() -> None:
    """MinTreePD with adapted Γ still matches the baseline (small tree)."""
    network = build_small_tree_network()
    te = _canonical_extension(network)
    taxa = {"a", "c", "d"}
    expected = SMALL_TREE_FIXED_DIVERSITY[frozenset(taxa)]
    value = pp.min_tree.compute_diversity(
        network,
        taxa,
        tree_extension=te,
    )
    assert value == expected
    assert math.isfinite(value)


def test_adapt_full_taxa_min_tree_matches_baseline() -> None:
    """All taxa: adapt is a no-op structurally; MinTreePD matches baseline."""
    network = build_small_tree_network()
    taxa = set(network.taxa)
    te = _canonical_extension(network)
    expected = SMALL_TREE_FIXED_DIVERSITY[frozenset(taxa)]
    value = pp.min_tree.compute_diversity(
        network,
        taxa,
        tree_extension=te,
    )
    assert value == expected
