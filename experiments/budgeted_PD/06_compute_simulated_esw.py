"""
Binary resolutions and optimal edge-scanwidth extensions for the PaNDA comparison.

This script must be run with ``phypanda`` version ``1.0.3`` (the PaNDA software of
the first paper), installed in a separate environment: ``pip install phypanda==1.0.3``.

A subsample of ``simulated_networks.csv`` is used: for each group of networks with the same
number of leaves and level (100 networks per group), ``SAMPLE_PER_GROUP`` networks
are drawn uniformly at random (seed ``RNG_SEED``), giving 1600 networks in total.

For each sampled network:
  - The eNewick string is parsed with ``DirectedNetwork.load_from_enewick``.
  - A binary resolution is computed with :func:`binary_resolution` below. This is
    the caterpillar construction of ``DirectedNetwork.binary_resolution`` from
    phypanda 1.0.3 (new edges get length 0), with two changes: the resolution of
    nodes with indegree larger than 2 is corrected (phypanda 1.0.3 raises an error
    on such nodes, which occur in 513 of the 6400 networks), and the new nodes get
    deterministic names ``res0``, ``res1``, ... instead of random ones.
  - An optimal edge-scanwidth extension of the binary resolution is computed with
    ``DAG.optimal_scanwidth`` (the default ``"optimal_XP"`` method that
    ``solve_MAPPD`` of phypanda 1.0.3 uses).

Output: ``simulated_esw_extensions.csv`` (space-delimited) with columns ``id``, ``esw``
(edge-scanwidth of the binary resolution), ``time_net`` (parsing the eNewick string),
``time_res`` (binary resolution), ``comp_time_xp`` (edge-scanwidth computation),
``resolution`` (the binary resolution as a list of ``(tail,head,length)`` arcs) and
``extension`` (a comma-separated ordering of nodes of the resolution forming the
extension). All times are in seconds.
"""

from __future__ import annotations

import csv
import random
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import networkx as nx
from phypanda import DAG, DirectedNetwork

SAMPLE_PER_GROUP: int = 25
RNG_SEED: int = 42


def sample_ids(csv_path: Path) -> list[str]:
    """Draw ``SAMPLE_PER_GROUP`` network ids per ``(ntips, level)`` group, sorted."""
    groups: dict[tuple[int, int], list[str]] = defaultdict(list)
    with csv_path.open(encoding="utf-8") as f:
        next(f, None)
        for line in f:
            parts = line.strip().split(None, 3)
            if len(parts) < 4:
                continue
            groups[(int(parts[1]), int(parts[2]))].append(parts[0])

    rng = random.Random(RNG_SEED)
    sampled: list[str] = []
    for key in sorted(groups):
        sampled.extend(rng.sample(sorted(groups[key]), SAMPLE_PER_GROUP))
    return sorted(sampled)


def binary_resolution(network: DirectedNetwork) -> DirectedNetwork:
    """
    Return a binary resolution of ``network``, with new edges of length 0.

    Every node with outdegree ``m > 2`` is replaced by a caterpillar: its first
    child stays, and the other children hang off a path of ``m - 2`` new nodes.
    Every node with indegree ``m > 2`` is resolved symmetrically: its first parent
    stays, and the other parents enter a path of ``m - 2`` new nodes that ends in
    the node. Children and parents are ordered by their name, and new nodes are
    named ``res0``, ``res1``, ..., so the result is deterministic.

    Parameters
    ----------
    network : DirectedNetwork
        Input network (phypanda 1.0.3 object); it is not modified.

    Returns
    -------
    DirectedNetwork
        Binary network with the same leaves and the same all-paths diversity.
    """
    edges = {(u, v): network.weight(u, v) for (u, v) in network.edges}
    counter = 0

    def new_node() -> str:
        nonlocal counter
        name = f"res{counter}"
        counter += 1
        assert name not in network.nodes
        return name

    for u in sorted(network.nodes, key=str):
        children = sorted(network.successors(u), key=str)
        if len(children) <= 2:
            continue
        last = u
        for i, c in enumerate(children[1:], start=1):
            w = edges.pop((u, c))
            if i < len(children) - 1:
                r = new_node()
                edges[(last, r)] = 0.0
                edges[(r, c)] = w
                last = r
            else:
                edges[(last, c)] = w

    # Parents after the outdegree step (a parent may now be one of the new nodes).
    current_parents: dict[Any, list[Any]] = defaultdict(list)
    for (u, v) in edges:
        current_parents[v].append(u)
    for v in sorted(network.nodes, key=str):
        parents = sorted(current_parents[v], key=str)
        if len(parents) <= 2:
            continue
        last = v
        for i, p in enumerate(parents[1:], start=1):
            w = edges.pop((p, v))
            if i < len(parents) - 1:
                r = new_node()
                edges[(r, last)] = 0.0
                edges[(p, r)] = w
                last = r
            else:
                edges[(p, last)] = w

    res = DirectedNetwork()
    for (u, v), w in edges.items():
        res.add_edge(u, v, weight=w)
    res.add_leaves_from([v for v in res.nodes if res.out_degree(v) == 0])
    return res


def resolution_to_string(network: DirectedNetwork) -> str:
    """Encode a network as ``(tail,head,length)`` tokens, sorted for reproducibility."""
    return ",".join(
        f"({u},{v},{network.weight(u, v)!r})"
        for (u, v) in sorted(network.edges, key=lambda e: (str(e[0]), str(e[1])))
    )


def _run_edgewidth_extensions_csv() -> None:
    """Read ``simulated_networks.csv`` and write ``simulated_esw_extensions.csv`` for the sample."""
    base = Path(__file__).resolve().parent
    src = base / "simulated_networks.csv"
    dst = base / "simulated_esw_extensions.csv"

    wanted = set(sample_ids(src))
    n_done = 0
    with src.open(encoding="utf-8") as f_in, dst.open(
        "w", encoding="utf-8", newline=""
    ) as f_out:
        writer = csv.writer(f_out, delimiter=" ")
        writer.writerow(
            ["id", "esw", "time_net", "time_res", "comp_time_xp", "resolution", "extension"]
        )
        for raw in f_in:
            line = raw.strip()
            if not line:
                continue
            parts = line.split(maxsplit=3)
            if parts[0] == "id" or len(parts) < 4 or parts[0] not in wanted:
                continue
            net_id, ntips, _level, newick = parts

            t0 = time.perf_counter()
            network = DirectedNetwork()
            network.load_from_enewick(newick.strip())
            t_net = time.perf_counter() - t0
            assert len(network.leaves) == int(ntips)

            t0 = time.perf_counter()
            resolution = binary_resolution(network)
            t_res = time.perf_counter() - t0
            assert resolution.is_binary()
            assert nx.is_directed_acyclic_graph(resolution)
            assert set(resolution.leaves) == set(network.leaves)

            t0 = time.perf_counter()
            esw, extension = DAG(nx.DiGraph(resolution)).optimal_scanwidth()
            t_xp = time.perf_counter() - t0

            writer.writerow([
                net_id, esw, f"{t_net:.9f}", f"{t_res:.9f}", f"{t_xp:.9f}",
                resolution_to_string(resolution),
                ",".join(str(v) for v in extension.sigma),
            ])
            n_done += 1
            if n_done % 50 == 0:
                print(f"Processed {n_done} networks...", flush=True)

    print(f"Done. Wrote {dst} ({n_done} rows).")


if __name__ == "__main__":
    _run_edgewidth_extensions_csv()
