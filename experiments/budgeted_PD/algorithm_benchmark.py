"""
Benchmark PaNDA all-paths (NSW FPT), MaxTree PD, and MinTree PD on ``nonbinary_nets.csv``.

Reads ``extension.csv`` and ``costs.csv`` alongside ``nonbinary_nets.csv``.  Networks
are processed in **level order** (all ``level`` 0 from ``nonbinary_nets.csv``, then
level 1, …); within a level, by ``id``.  One row per network (space-delimited
columns) with prep times plus each algorithm run time and objective value.

Budgets for all-paths and MaxTree: ``25%``, ``50%``, and ``90%`` of total taxon cost
(each ``round(fraction * total_taxon_cost)``, clamped to ``[0, total]``).  MinTree PD
uses the full taxon set.

Output path: ``algorithm_benchmark_results.csv`` in the same directory as this script.
"""

from __future__ import annotations

import csv
import re
import sys
import time
from pathlib import Path
from typing import Any


def _repo_root() -> Path:
    """Walk up from this file until ``src/phypanda`` exists."""
    p = Path(__file__).resolve().parent
    for _ in range(8):
        if (p / "src" / "phypanda").is_dir():
            return p
        if p.parent == p:
            break
        p = p.parent
    raise RuntimeError("Cannot find repository root (expected src/phypanda).")


_ROOT = _repo_root()
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

import phypanda as pp
from phylozoo import DirectedPhyNetwork
from scanwidth import DAG, Extension

_PAIR_RE = re.compile(r"\(([^,]+),(\d+)\)")

# Short suffixes: b25/b50/b90 = 25% / 50% / 90% of total taxon cost (rounded).
_BUDGET_TAGS = ("b25", "b50", "b90")

# CSV columns in order (space-delimited); keep in sync with ``_wide_fieldnames``:
#   id, level, tot_cost, time_net, time_ext,
#   time_ap_b25, val_ap_b25, time_ap_b50, val_ap_b50, time_ap_b90, val_ap_b90,
#   time_mxt_b25, val_mxt_b25, time_mxt_b50, val_mxt_b50, time_mxt_b90, val_mxt_b90,
#   time_mnt_full, val_mnt_full


def _load_rows_sorted_by_level(csv_path: Path) -> list[tuple[str, int, str]]:
    """Load networks as ``(net_id, level, newick)``, sorted by level then id."""
    rows: list[tuple[str, int, str]] = []
    with csv_path.open(encoding="utf-8") as f:
        next(f, None)
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(None, 3)
            if len(parts) < 4:
                continue
            net_id = parts[0].strip()
            level = int(parts[2])
            newick = parts[3].strip()
            rows.append((net_id, level, newick))
    rows.sort(key=lambda r: (r[1], r[0]))
    return rows


def _load_extension_map(csv_path: Path) -> dict[str, str]:
    """Map network id to the comma-separated extension vertex order string."""
    out: dict[str, str] = {}
    with csv_path.open(encoding="utf-8") as f:
        next(f, None)
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(None, 3)
            if len(parts) < 4:
                continue
            out[parts[0].strip()] = parts[3].strip()
    return out


def _load_costs_map(csv_path: Path) -> dict[str, dict[str, int]]:
    """Map network id to taxon label -> integer cost (from ``(tax,cost)`` tokens)."""
    out: dict[str, dict[str, int]] = {}
    with csv_path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            net_id, rest = line.split(None, 1)
            out[net_id] = {tax: int(c) for tax, c in _PAIR_RE.findall(rest)}
    return out


def _label_to_node(network: DirectedPhyNetwork) -> dict[str, Any]:
    """Map each node label string (or fallback node id) to the graph node object."""
    lut: dict[str, Any] = {}
    for node in network._graph.nodes():
        lab = network.get_label(node)
        key = str(lab) if lab is not None else str(node)
        lut[key] = node
    return lut


def _tree_extension_from_labels(network: DirectedPhyNetwork, labels_csv: str) -> Any:
    """Build a canonical tree extension from the NSW ordering string for this network."""
    dag = DAG(network._graph._graph)
    lut = _label_to_node(network)
    ordering = [lut[s.strip()] for s in labels_csv.split(",")]
    return Extension(dag, ordering).to_canonical_tree_extension()


def _wide_fieldnames() -> list[str]:
    """Header field names for the wide benchmark result row."""
    cols = ["id", "level", "tot_cost", "time_net", "time_ext"]
    for tag in _BUDGET_TAGS:
        cols.append(f"time_ap_{tag}")
        cols.append(f"val_ap_{tag}")
    for tag in _BUDGET_TAGS:
        cols.append(f"time_mxt_{tag}")
        cols.append(f"val_mxt_{tag}")
    cols.extend(["time_mnt_full", "val_mnt_full"])
    return cols


def main() -> None:
    """Process every network and write ``algorithm_benchmark_results.csv`` beside this file."""
    here = Path(__file__).resolve().parent
    nets_path = here / "nonbinary_nets.csv"
    out_path = here / "algorithm_benchmark_results.csv"
    extension_map = _load_extension_map(here / "extension.csv")
    costs_map = _load_costs_map(here / "costs.csv")
    fieldnames = _wide_fieldnames()
    networks_sorted = _load_rows_sorted_by_level(nets_path)

    n_done = 0
    with out_path.open("w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter=" ")
        writer.writeheader()

        for net_id, level, newick in networks_sorted:
            ext_str = extension_map[net_id]
            costs_raw = costs_map[net_id]

            t0 = time.perf_counter()
            network = DirectedPhyNetwork.from_string(newick)
            sec_load = time.perf_counter() - t0

            taxa_set = set(network.taxa)
            costs_full = {t: int(costs_raw.get(t, 1)) for t in taxa_set}
            total_cost = sum(costs_full.values())

            t0 = time.perf_counter()
            te = _tree_extension_from_labels(network, ext_str)
            sec_te = time.perf_counter() - t0

            budgets = (
                max(0, min(total_cost, int(round(0.25 * total_cost)))),
                max(0, min(total_cost, int(round(0.50 * total_cost)))),
                max(0, min(total_cost, int(round(0.90 * total_cost)))),
            )

            row: dict[str, Any] = {
                "id": net_id,
                "level": level,
                "tot_cost": total_cost,
                "time_net": f"{sec_load:.9f}",
                "time_ext": f"{sec_te:.9f}",
            }

            for tag, budget in zip(_BUDGET_TAGS, budgets):
                t0 = time.perf_counter()
                v, _ = pp.all_paths.solve_maximization(
                    network,
                    budget=budget,
                    costs=costs_full,
                    algorithm="nsw_fpt_budget",
                    tree_extension=te,
                    numba=True,
                )
                dt = time.perf_counter() - t0
                row[f"time_ap_{tag}"] = f"{dt:.9f}"
                row[f"val_ap_{tag}"] = f"{float(v):.17g}"

            for tag, budget in zip(_BUDGET_TAGS, budgets):
                t0 = time.perf_counter()
                v, _ = pp.max_tree.solve_maximization(
                    network,
                    budget=budget,
                    costs=costs_full,
                    tree_extension=te,
                    numba=True,
                )
                dt = time.perf_counter() - t0
                row[f"time_mxt_{tag}"] = f"{dt:.9f}"
                row[f"val_mxt_{tag}"] = f"{float(v):.17g}"

            t0 = time.perf_counter()
            v = pp.min_tree.compute_diversity(network, taxa_set, tree_extension=te)
            row["time_mnt_full"] = f"{time.perf_counter() - t0:.9f}"
            row["val_mnt_full"] = f"{float(v):.17g}"

            writer.writerow(row)
            n_done += 1
            if n_done % 50 == 0:
                print(f"Processed {n_done} networks...", flush=True)

    print(f"Done. Wrote {out_path} ({n_done} rows).")


if __name__ == "__main__":
    main()
