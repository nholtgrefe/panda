"""
Benchmark PaNDA (phypanda 1.0.3) on unit-cost instances of the sampled networks.

This script must be run with ``phypanda`` version ``1.0.3``, like
``06_compute_simulated_esw.py``.

Reads ``simulated_esw_extensions.csv`` (binary resolutions and their optimal
edge-scanwidth extensions) and ``simulated_networks.csv`` (for the number of leaves and
level). Networks are processed in **level order**, and within a level by ``id``,
as in ``04_run_simulated_experiments.py``.

PaNDA solves MAPPD, i.e. maximizes all-paths diversity with unit taxon costs. The
budgets are ``k = round(0.25 n)``, ``round(0.5 n)`` and ``round(0.9 n)``, with ``n``
the number of leaves (the same fractions of the total taxon cost as in
``04_run_simulated_experiments.py``). Each budget is solved with ``solve_MAPPD`` on the binary
resolution, given the canonical tree extension of the saved extension (so the
edge-scanwidth computation, timed in ``06``, is not repeated here).

Output: ``simulated_esw_results.csv`` (space-delimited) with columns ``id``, ``level``,
``ntips``, ``esw``, ``time_net`` (building the binary resolution from its arcs),
``time_ext`` (building the canonical tree extension), and for each budget
``time_panda_<b>`` (solve time in seconds) and ``val_panda_<b>`` (diversity value),
with ``<b>`` in ``b25``, ``b50``, ``b90``.

The same unit-cost instances are solved with the node-scanwidth algorithm in
``04_run_simulated_experiments.py`` (columns ``time_apu_<b>`` and ``val_apu_<b>`` of
``simulated_results.csv``); the values ``val_panda_<b>`` and ``val_apu_<b>``
should be equal.
"""

from __future__ import annotations

import csv
import re
import time
from pathlib import Path
from typing import Any

import networkx as nx
from phypanda import DirectedNetwork, solve_MAPPD
from phypanda.scanwidth import Extension

_ARC_RE = re.compile(r"\(([^,()]+),([^,()]+),([^,()]+)\)")

# Short suffixes: b25/b50/b90 = 25% / 50% / 90% of the number of leaves (rounded).
_BUDGET_TAGS = ("b25", "b50", "b90")
_BUDGET_FRACTIONS = (0.25, 0.50, 0.90)


def _load_network_info(csv_path: Path) -> dict[str, tuple[int, int]]:
    """Map network id to ``(ntips, level)`` from ``simulated_networks.csv``."""
    out: dict[str, tuple[int, int]] = {}
    with csv_path.open(encoding="utf-8") as f:
        next(f, None)
        for line in f:
            parts = line.strip().split(None, 3)
            if len(parts) < 4:
                continue
            out[parts[0]] = (int(parts[1]), int(parts[2]))
    return out


def _load_extensions(csv_path: Path) -> list[tuple[str, int, str, str]]:
    """Load ``(id, esw, resolution, extension)`` rows from ``simulated_esw_extensions.csv``."""
    rows: list[tuple[str, int, str, str]] = []
    with csv_path.open(encoding="utf-8") as f:
        next(f, None)
        for line in f:
            parts = line.strip().split(None, 6)
            if len(parts) < 7:
                continue
            rows.append((parts[0], int(parts[1]), parts[5], parts[6]))
    return rows


def _network_from_arcs(arcs: str) -> DirectedNetwork:
    """Build a phypanda 1.0.3 network from ``(tail,head,length)`` tokens."""
    network = DirectedNetwork()
    for u, v, w in _ARC_RE.findall(arcs):
        network.add_edge(u, v, weight=float(w))
    network.add_leaves_from([v for v in network.nodes if network.out_degree(v) == 0])
    return network


def _tree_extension_from_labels(network: DirectedNetwork, labels_csv: str) -> Any:
    """Build the canonical tree extension from the saved extension ordering."""
    ordering = [s.strip() for s in labels_csv.split(",")]
    return Extension(nx.DiGraph(network), ordering).canonical_tree_extension()


def _wide_fieldnames() -> list[str]:
    """Header field names for the wide benchmark result row."""
    cols = ["id", "level", "ntips", "esw", "time_net", "time_ext"]
    for tag in _BUDGET_TAGS:
        cols.append(f"time_panda_{tag}")
        cols.append(f"val_panda_{tag}")
    return cols


def main() -> None:
    """Process every sampled network and write ``simulated_esw_results.csv`` beside this file."""
    here = Path(__file__).resolve().parent
    out_path = here / "simulated_esw_results.csv"
    info = _load_network_info(here / "simulated_networks.csv")
    rows = _load_extensions(here / "simulated_esw_extensions.csv")
    rows.sort(key=lambda r: (info[r[0]][1], r[0]))
    fieldnames = _wide_fieldnames()

    n_done = 0
    with out_path.open("w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter=" ")
        writer.writeheader()

        for net_id, esw, arcs, ext_str in rows:
            ntips, level = info[net_id]

            t0 = time.perf_counter()
            network = _network_from_arcs(arcs)
            sec_net = time.perf_counter() - t0
            assert network.is_binary() and len(network.leaves) == ntips

            t0 = time.perf_counter()
            te = _tree_extension_from_labels(network, ext_str)
            sec_te = time.perf_counter() - t0
            assert te.scanwidth() == esw

            row: dict[str, Any] = {
                "id": net_id,
                "level": level,
                "ntips": ntips,
                "esw": esw,
                "time_net": f"{sec_net:.9f}",
                "time_ext": f"{sec_te:.9f}",
            }

            for tag, frac in zip(_BUDGET_TAGS, _BUDGET_FRACTIONS):
                k = max(0, min(ntips, int(round(frac * ntips))))
                t0 = time.perf_counter()
                v, _ = solve_MAPPD(network, k, tree_extension=te)
                dt = time.perf_counter() - t0
                row[f"time_panda_{tag}"] = f"{dt:.9f}"
                row[f"val_panda_{tag}"] = f"{float(v):.17g}"

            writer.writerow(row)
            f_out.flush()
            n_done += 1
            if n_done % 50 == 0:
                print(f"Processed {n_done} networks...", flush=True)

    print(f"Done. Wrote {out_path} ({n_done} rows).")


if __name__ == "__main__":
    main()
