"""
Base functions for diversity calculations.
"""

from __future__ import annotations

from typing import Any, Dict, Mapping, Set

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooValueError

from .protocol import DiversityMeasure


def _normalize_costs(
    network: DirectedPhyNetwork,
    costs: Mapping[str, int] | None,
) -> dict[str, int]:
    """
    Normalize and validate integer taxon costs.

    Parameters
    ----------
    network : DirectedPhyNetwork
        Phylogenetic network providing the taxon universe.
    costs : Mapping[str, int] | None
        Optional taxon cost mapping. If ``None``, unit costs are used.

    Returns
    -------
    dict[str, int]
        Cost mapping for all taxa in ``network``.

    Raises
    ------
    PhyloZooValueError
        If costs are missing for any taxon, non-integer, or negative.
    """
    taxa = set(network.taxa)
    if costs is None:
        return {taxon: 1 for taxon in taxa}

    missing = taxa - set(costs.keys())
    if missing:
        raise PhyloZooValueError(f"Missing costs for taxa: {missing}")

    normalized: dict[str, int] = {}
    for taxon in taxa:
        cost = costs[taxon]
        if not isinstance(cost, int):
            raise PhyloZooValueError(
                f"Cost for taxon '{taxon}' must be an integer, got {type(cost).__name__}"
            )
        if cost < 0:
            raise PhyloZooValueError(
                f"Cost for taxon '{taxon}' must be non-negative, got {cost}"
            )
        normalized[taxon] = cost
    return normalized


def _validate_budget(budget: int) -> None:
    """
    Validate budget as a non-negative integer.

    Parameters
    ----------
    budget : int
        Total available budget.

    Raises
    ------
    PhyloZooValueError
        If ``budget`` is not an integer or is negative.
    """
    if not isinstance(budget, int):
        raise PhyloZooValueError(
            f"budget must be an integer, got {type(budget).__name__}"
        )
    if budget < 0:
        raise PhyloZooValueError(f"budget must be non-negative, got {budget}")


def diversity(
    network: DirectedPhyNetwork,
    taxa: set[str],
    measure: DiversityMeasure,
    **kwargs: Any,
) -> float:
    """
    Compute diversity for a fixed taxon set.

    Parameters
    ----------
    network : DirectedPhyNetwork
        Input phylogenetic network.
    taxa : set[str]
        Taxa for which diversity is evaluated.
    measure : DiversityMeasure
        Diversity measure implementation.
    **kwargs : Any
        Measure-specific keyword arguments.

    Returns
    -------
    float
        Diversity value for ``taxa``.

    Raises
    ------
    PhyloZooValueError
        If one or more taxa are not present in the network.

    Examples
    --------
    >>> import phypanda as pp
    >>> # value = pp.diversity(network, {"a", "b"}, measure=pp.all_paths)
    """
    network_taxa = set(network.taxa)
    if not taxa.issubset(network_taxa):
        missing = taxa - network_taxa
        raise PhyloZooValueError(f"Taxa not found in network: {missing}")
    return measure.compute_diversity(network, taxa, **kwargs)


def marginal_diversities(
    network: DirectedPhyNetwork,
    saved_taxa: set[str],
    measure: DiversityMeasure,
    **kwargs: Any,
) -> Dict[str, float]:
    """
    Compute marginal diversity contributions for all taxa.

    Parameters
    ----------
    network : DirectedPhyNetwork
        Input phylogenetic network.
    saved_taxa : set[str]
        Current selected taxa.
    measure : DiversityMeasure
        Diversity measure implementation.
    **kwargs : Any
        Measure-specific keyword arguments.

    Returns
    -------
    Dict[str, float]
        Mapping from taxon to marginal gain/loss relative to ``saved_taxa``.

    Examples
    --------
    >>> import phypanda as pp
    >>> # marg = pp.marginal_diversities(network, {"a"}, measure=pp.all_paths)
    """
    total_div = diversity(network, saved_taxa, measure, **kwargs)
    marginal: Dict[str, float] = {}
    all_taxa = set(network.taxa)

    for taxon in all_taxa:
        if taxon in saved_taxa:
            div_minus = diversity(network, saved_taxa - {taxon}, measure, **kwargs)
            marginal[taxon] = div_minus - total_div
        else:
            div_plus = diversity(network, saved_taxa | {taxon}, measure, **kwargs)
            marginal[taxon] = div_plus - total_div

    return marginal


def greedy_max_diversity(
    network: DirectedPhyNetwork,
    budget: int,
    measure: DiversityMeasure,
    costs: Mapping[str, int] | None = None,
    **kwargs: Any,
) -> tuple[float, Set[str]]:
    """
    Greedily maximize diversity under an integer budget.

    The selected taxon at each step is the affordable taxon with highest
    normalized gain ``marginal_diversity / cost``.

    Parameters
    ----------
    network : DirectedPhyNetwork
        Input phylogenetic network.
    budget : int
        Total integer budget.
    measure : DiversityMeasure
        Diversity measure implementation.
    costs : Mapping[str, int] | None, optional
        Integer cost per taxon. If ``None``, unit costs are used.
    **kwargs : Any
        Measure-specific keyword arguments.

    Returns
    -------
    tuple[float, Set[str]]
        Greedy objective value and selected taxa.

    Examples
    --------
    >>> import phypanda as pp
    >>> # value, taxa = pp.greedy_max_diversity(network, 5, measure=pp.all_paths)
    >>> # sorted(taxa)  # doctest: +SKIP
    """
    _validate_budget(budget)
    normalized_costs = _normalize_costs(network, costs)

    saved_taxa: set[str] = set()
    remaining_budget = budget
    all_taxa = set(network.taxa)

    while True:
        marginal = marginal_diversities(network, saved_taxa, measure, **kwargs)
        best_taxon: str | None = None
        best_score = -float("inf")

        for taxon in all_taxa:
            if taxon in saved_taxa:
                continue
            cost = normalized_costs[taxon]
            if cost > remaining_budget:
                continue
            gain = marginal[taxon]
            score = float("inf") if cost == 0 and gain > 0 else (gain / cost if cost else gain)
            if score > best_score:
                best_score = score
                best_taxon = taxon

        if best_taxon is None:
            break

        saved_taxa.add(best_taxon)
        remaining_budget -= normalized_costs[best_taxon]
        if remaining_budget == 0:
            break

    div_value = diversity(network, saved_taxa, measure, **kwargs)
    return div_value, saved_taxa


def solve_max_diversity(
    network: DirectedPhyNetwork,
    budget: int,
    measure: DiversityMeasure,
    costs: Mapping[str, int] | None = None,
    **kwargs: Any,
) -> tuple[float, Set[str]]:
    """
    Solve maximum diversity under an integer budget.

    Parameters
    ----------
    network : DirectedPhyNetwork
        Input phylogenetic network.
    budget : int
        Total integer budget.
    measure : DiversityMeasure
        Diversity measure implementation.
    costs : Mapping[str, int] | None, optional
        Integer cost per taxon. If ``None``, unit costs are used.
    **kwargs : Any
        Measure-specific keyword arguments.

    Returns
    -------
    tuple[float, Set[str]]
        Optimal objective value and selected taxa.

    Examples
    --------
    >>> import phypanda as pp
    >>> # value, taxa = pp.solve_max_diversity(network, 5, measure=pp.all_paths)
    >>> # len(taxa) <= 5  # doctest: +SKIP
    """
    _validate_budget(budget)
    normalized_costs = _normalize_costs(network, costs)
    return measure.solve_maximization(
        network,
        budget=budget,
        costs=normalized_costs,
        **kwargs,
    )
