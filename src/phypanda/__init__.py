"""
Panda module for phylogenetic diversity calculations.

This module provides a framework for computing various phylogenetic diversity
measures with a common interface.
"""

from .base import (
    diversity,
    greedy_max_diversity,
    marginal_diversities,
    solve_max_diversity,
)
from .measure import (
    average_tree,
    all_paths,
    AllPathsDiversity,
    network,
    NetworkDiversity,
    min_tree,
    MinTreeDiversity,
    max_tree,
    MaxTreeDiversity,
    tree,
    TreeDiversity,
    AverageTreeDiversity,
)
from .protocol import DiversityMeasure

__all__ = [
    "diversity",
    "marginal_diversities",
    "greedy_max_diversity",
    "solve_max_diversity",
    "all_paths",
    "AllPathsDiversity",
    "min_tree",
    "MinTreeDiversity",
    "max_tree",
    "MaxTreeDiversity",
    "tree",
    "TreeDiversity",
    "network",
    "NetworkDiversity",
    "average_tree",
    "AverageTreeDiversity",
    "DiversityMeasure",
]

__version__ = "2.0.0"