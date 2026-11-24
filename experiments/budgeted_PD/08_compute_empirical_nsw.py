"""
XP node-scanwidth of the empirical networks, with batch export of optimal extensions.

Reads ``empirical_networks.csv``: the 27 real phylogenetic networks used in the
experiments of the ``scanwidth`` paper (Holtgrefe, van Iersel and Jones, JCSS 2026;
originally from http://phylnet.univ-mlv.fr/), taken from the edge-list files in
``experiments/networks/real_networks`` of the ``scanwidth`` repository. Each row
contains ``name`` and ``arcs`` (a comma-separated list of ``(tail,head)`` arcs).
The networks are stored as arcs rather than eNewick strings, because several of
them are not phylogenetic networks in the sense of PhyloZoo (for example, they
contain vertices with indegree and outdegree 1).

The optimal node-scanwidth is computed on each network as given, with the same
settings as ``03_compute_simulated_nsw.py``. Output ``empirical_nsw_extensions.csv``
uses space-separated columns ``name``, ``nsw`` (node-scanwidth value),
``comp_time_xp`` (computation time in seconds), and ``extension`` (the vertex order
as comma-separated vertex names).
"""

from __future__ import annotations

import csv
import re
import time
from pathlib import Path

import networkx as nx
from scanwidth import DAG, node_scanwidth
from scanwidth.node_scanwidth.reduction.config import ReducerConfig

_ARC_RE = re.compile(r"\(([^,()]+),([^,()]+)\)")


def _dag_from_arcs(arcs: str) -> nx.DiGraph:
    """Build a directed graph from ``(tail,head)`` tokens."""
    graph = nx.DiGraph()
    graph.add_edges_from(_ARC_RE.findall(arcs))
    return graph


def _run_empirical_extensions_csv() -> None:
    """Read ``empirical_networks.csv`` and write ``empirical_nsw_extensions.csv``."""
    base = Path(__file__).resolve().parent
    src = base / "empirical_networks.csv"
    dst = base / "empirical_nsw_extensions.csv"

    with src.open(encoding="utf-8") as f_in, dst.open(
        "w", encoding="utf-8", newline=""
    ) as f_out:
        writer = csv.writer(f_out, delimiter=" ")
        writer.writerow(["name", "nsw", "comp_time_xp", "extension"])
        for raw in f_in:
            line = raw.strip()
            if not line:
                continue
            parts = line.split(maxsplit=1)
            if parts[0] == "name" or len(parts) < 2:
                continue
            name, arcs = parts
            graph = _dag_from_arcs(arcs)

            t0 = time.perf_counter()
            nsw, extension = node_scanwidth(
                DAG(graph),
                algorithm="xp",
                reduce=True,
                reducer_config=ReducerConfig(parallel_sblocks=False),
            )
            t_xp = time.perf_counter() - t0

            writer.writerow(
                [name, nsw, f"{t_xp:.9f}", ",".join(str(v) for v in extension.ordering)]
            )


if __name__ == "__main__":
    _run_empirical_extensions_csv()
