"""
General utility helpers for the phypanda package.
"""

from __future__ import annotations

from itertools import combinations
from typing import Any, Iterator


def powerset(s: set[Any]) -> Iterator[set[Any]]:
    """
    Yield all subsets of a set.

    Parameters
    ----------
    s : set[Any]
        Input set.

    Yields
    ------
    Iterator[set[Any]]
        Subsets of ``s``.
    """
    s_list = tuple(sorted(s, key=str))
    for r in range(len(s_list) + 1):
        for subset in combinations(s_list, r):
            yield set(subset)
