"""
Contract a share of the shortest eligible branches on a PhyloZoo directed network.

Eligibility (for a directed edge ``(u, v)``): skip if ``v`` is a leaf; skip if
``indegree(u) >= 2`` and ``outdegree(v) >= 2``. A contraction is applied only if
merging ``v`` into ``u`` keeps every internal node in a valid tree/hybrid degree
class for :class:`~phylozoo.core.network.dnetwork.DirectedPhyNetwork` and does
not introduce parallel arcs (multiple edges with the same head and tail).

This script was run with phylozoo version 0.1.2.
"""

from __future__ import annotations

from typing import Any
import csv
from pathlib import Path

import networkx as nx
from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.conversions import dnetwork_from_graph
from phylozoo.core.primitives.d_multigraph.conversions import (
    multidigraph_to_directedmultigraph,
)
from phylozoo.core.primitives.d_multigraph.features import has_parallel_edges


def _directed_phylo_internal_degrees_ok(graph: nx.MultiDiGraph) -> bool:
    """True if each non-root internal node is tree-type or hybrid-type."""
    for n in graph.nodes():
        ins, outs = graph.in_degree(n), graph.out_degree(n)
        if outs == 0 or ins == 0:
            continue
        if ins == 1 and outs >= 2:
            continue
        if ins >= 2 and outs == 1:
            continue
        return False
    return True


def _has_parallel_arcs(graph: nx.MultiDiGraph) -> bool:
    """True if some ordered pair has more than one directed edge (PhyloZoo helper)."""
    return has_parallel_edges(multidigraph_to_directedmultigraph(graph))


def _eligible_edges(graph: nx.MultiDiGraph) -> list[tuple[Any, Any, float]]:
    """(u, v, length) for edges that pass the screening rule (one entry per arc)."""
    out: list[tuple[Any, Any, float]] = []
    for u, v, _k, data in graph.edges(keys=True, data=True):
        if graph.out_degree(v) == 0:
            continue
        if graph.in_degree(u) >= 2 and graph.out_degree(v) >= 2:
            continue
        bl = data.get("branch_length")
        length = 1.0 if bl is None else float(bl)
        out.append((u, v, length))
    return out


def _contract_edge_safe(graph: nx.MultiDiGraph, u: Any, v: Any) -> bool:
    """Merge ``v`` into ``u`` if some ``(u,v)`` arc exists and the result is valid."""
    if not graph.has_edge(u, v) or graph.out_degree(v) == 0:
        return False
    trial = graph.copy()
    nx.contracted_nodes(trial, u, v, self_loops=False, copy=False)
    if trial.has_edge(u, u):
        trial.remove_edge(u, u)
    if _has_parallel_arcs(trial):
        return False
    if not _directed_phylo_internal_degrees_ok(trial):
        return False
    nx.contracted_nodes(graph, u, v, self_loops=False, copy=False)
    if graph.has_edge(u, u):
        graph.remove_edge(u, u)
    return True


def contract_percent_shortest_edges(
    network: DirectedPhyNetwork,
    percent: float,
) -> DirectedPhyNetwork:
    """
    Contract a percentage of the shortest eligible branches.

    Eligible edges are sorted by ``branch_length`` (default 1 if missing). Up to
    ``floor(n * percent / 100)`` are attempted in that order; skipped if invalid.

    Parameters
    ----------
    network : DirectedPhyNetwork
        Input network.
    percent : float
        Percentage in ``[0, 100]`` of eligible edges to try to contract.

    Returns
    -------
    DirectedPhyNetwork
        New network (validated on construction).

    Raises
    ------
    ValueError
        If ``percent`` is not finite or not in ``[0, 100]``.
    """
    if not (0.0 <= percent <= 100.0) or percent != percent:
        raise ValueError("percent must be a finite number in [0, 100].")

    graph = network._graph._graph.copy()
    eligible = _eligible_edges(graph)
    eligible.sort(key=lambda t: t[2])
    target = int(len(eligible) * percent / 100.0)

    done = 0
    for u, v, _ in eligible:
        if done >= target:
            break
        if _contract_edge_safe(graph, u, v):
            done += 1

    return dnetwork_from_graph(graph)


def _repo_root() -> Path:
    p = Path(__file__).resolve().parent
    for _ in range(8):
        if (p / "src" / "phypanda").is_dir():
            return p
        if p.parent == p:
            break
        p = p.parent
    raise RuntimeError("Cannot find repository root (expected src/phypanda).")


def _run_exp1_nonbinary_csv() -> None:
    """Write ``nonbinary_nets.csv``: exp1 networks after 10% shortest-edge contraction."""

    root = _repo_root()
    src = root / "experiments" / "MAPPD" / "exp1_simulated_networks.csv"
    dst = Path(__file__).resolve().parent / "nonbinary_nets.csv"

    with src.open(encoding="utf-8") as f_in, dst.open(
        "w", encoding="utf-8", newline=""
    ) as f_out:
        writer = csv.writer(f_out, delimiter=" ")
        writer.writerow(["id", "ntips", "level", "newick"])
        for raw in f_in:
            line = raw.strip()
            if not line:
                continue
            parts = line.split(maxsplit=6)
            if parts[0] == "id" or len(parts) < 7:
                continue
            net_id, _ntips_in, level_s, _nu, _mu, _lam, newick = parts
            try:
                net = DirectedPhyNetwork.from_string(newick.strip())
                out = contract_percent_shortest_edges(net, 10.0)
            except Exception as exc:
                print(f"{net_id}: skipped ({exc})")
                continue
            writer.writerow(
                [net_id, str(len(out.taxa)), level_s, out.to_string()]
            )


if __name__ == "__main__":
    _run_exp1_nonbinary_csv()
