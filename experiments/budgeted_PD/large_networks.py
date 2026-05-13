"""
Large-network benchmark: merge multiple sub-networks by root identification and
run all-paths (NSW FPT), MaxTree PD, and MinTree PD.

For each of N_TRIALS trials:
  - Randomly sample N_NETWORKS_PER_TRIAL network IDs from [4801, 6400] with
    replacement (IDs formatted as 5-digit strings, e.g. "05312").
  - Load each sub-network from ``nonbinary_nets.csv``.
  - Merge them by identifying all sub-network roots into a single shared root
    vertex (no super-root is added).  Taxa are relabeled with a unique prefix
    to prevent collisions across sub-networks.
  - Assign log-normal costs (mu=2.0, sigma=0.8) to the combined taxon set.
  - Compute a tree extension via node_scanwidth.
  - Run algorithms at budgets 25%, 50%, 90% (all-paths, MaxTree) and on the
    full taxon set (MinTree).

Prints per-trial results plus summary statistics (mean, median, mode, std)
and writes all results to ``large_network_results.csv`` in this directory.
"""

from __future__ import annotations

import csv
import statistics
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
from phylozoo import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.classifications import level as net_level
from scanwidth import DAG, node_scanwidth


def _repo_root() -> Path:
    p = Path(__file__).resolve().parent
    for _ in range(8):
        if (p / "src" / "phypanda").is_dir():
            return p
        if p.parent == p:
            break
        p = p.parent
    raise RuntimeError("Cannot find repository root.")


_ROOT = _repo_root()
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

import phypanda as pp  # noqa: E402

# ── Parameters ────────────────────────────────────────────────────────────────
N_TRIALS: int = 20
N_NETWORKS_PER_TRIAL: int = 5
ID_LOW: int = 4801
ID_HIGH: int = 6400
MU: float = 2.0
SIGMA: float = 0.8
RNG_SEED: int = 42
SHARED_ROOT_ID: str = "__root__"

_BUDGET_TAGS = ("b25", "b50", "b90")
_BUDGET_FRACTIONS = (0.25, 0.50, 0.90)
# ─────────────────────────────────────────────────────────────────────────────


def _read_newick_by_id(net_id: str, csv_path: Path) -> str:
    with csv_path.open(encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter=" ")
        next(reader, None)
        for row in reader:
            parts = [p for p in row if p]
            if len(parts) >= 4 and parts[0] == net_id:
                return parts[3]
    raise ValueError(f"Network id '{net_id}' not found in {csv_path}.")


def _build_combined_network(
    indexed_networks: list[tuple[str, DirectedPhyNetwork]],
) -> DirectedPhyNetwork:
    """Merge sub-networks by identifying all roots into one shared root vertex.

    Each sub-network node is relabeled ``{compound_id}_{original_node}``; the
    root node of each sub-network maps to SHARED_ROOT_ID instead, so all
    sub-network roots become a single vertex.  Taxon labels get the same
    ``{compound_id}_`` prefix to guarantee uniqueness.  The compound id is
    ``{i}_{nid}`` so that duplicate sampled IDs remain distinct.
    """
    all_nodes: list[tuple[Any, dict]] = [(SHARED_ROOT_ID, {})]
    all_edges = []

    for compound_id, network in indexed_networks:
        g = network._graph._graph
        root = network.root_node

        for node, data in g.nodes(data=True):
            if node == root:
                continue
            new_id = f"{compound_id}_{node}"
            label = data.get("label")
            if label is not None:
                all_nodes.append((new_id, {"label": f"{compound_id}_{label}"}))
            else:
                all_nodes.append((new_id, {}))

        for u, v, _key, edata in g.edges(data=True, keys=True):
            new_u = SHARED_ROOT_ID if u == root else f"{compound_id}_{u}"
            new_v = SHARED_ROOT_ID if v == root else f"{compound_id}_{v}"
            edge: dict[str, Any] = {"u": new_u, "v": new_v}
            bl = edata.get("branch_length")
            if bl is not None:
                edge["branch_length"] = bl
            all_edges.append(edge)

    return DirectedPhyNetwork(edges=all_edges, nodes=all_nodes)


def _assign_lognormal_costs(
    taxa: list[str], *, rng: np.random.Generator
) -> dict[str, int]:
    samples = rng.lognormal(mean=MU, sigma=SIGMA, size=len(taxa))
    return {
        taxon: max(1, int(np.rint(sample)))
        for taxon, sample in zip(taxa, samples)
    }


def _compute_tree_extension(
    network: DirectedPhyNetwork,
) -> tuple[Any, int, float, float]:
    """Compute a canonical tree extension for *network*.

    Two steps are timed separately:
      1. ``node_scanwidth`` — the NSW algorithm itself (FPT in node scanwidth).
      2. ``to_canonical_tree_extension`` — converts the NSW ordering into the
         canonical tree extension object consumed by the PD solvers.

    Returns
    -------
    (tree_extension, nsw_value, sec_nsw, sec_te)
        The tree extension, the integer node scanwidth value, and wall-clock
        seconds for each step.
    """
    dag = DAG(network._graph._graph)

    # Step 1: run NSW algorithm to obtain a vertex ordering.
    t0 = time.perf_counter()
    nsw_value, extension = node_scanwidth(dag)
    sec_nsw = time.perf_counter() - t0

    if extension is None:
        raise RuntimeError("node_scanwidth failed to produce a tree extension.")

    # Step 2: convert the NSW ordering into a canonical tree extension.
    t0 = time.perf_counter()
    tree_ext = extension.to_canonical_tree_extension()
    sec_te = time.perf_counter() - t0

    return tree_ext, nsw_value, sec_nsw, sec_te


def _summary_stats(values: list[float]) -> dict[str, float]:
    mean = statistics.mean(values)
    median = statistics.median(values)
    std = statistics.stdev(values) if len(values) > 1 else 0.0
    try:
        mode_val = float(statistics.mode(values))
    except statistics.StatisticsError:
        mode_val = float("nan")
    return {"mean": mean, "median": median, "mode": mode_val, "std": std}


def _fieldnames() -> list[str]:
    # Network identity and topology info.
    # time_nsw   = wall-clock seconds for node_scanwidth itself
    # time_te    = wall-clock seconds for to_canonical_tree_extension conversion
    # nsw_value  = integer node scanwidth of the combined network
    # costs and newick are placed last; newick is the free-form trailing field
    # (same convention as nonbinary_nets.csv).
    cols = ["network_ids", "n_taxa", "tot_cost", "level", "nsw_value",
            "time_nsw", "time_te"]
    for tag in _BUDGET_TAGS:
        cols += [f"time_ap_{tag}", f"val_ap_{tag}"]
    for tag in _BUDGET_TAGS:
        cols += [f"time_mxt_{tag}", f"val_mxt_{tag}"]
    cols += ["time_mnt_full", "val_mnt_full", "costs", "newick"]
    return cols


def _numeric_fieldnames(fieldnames: list[str]) -> list[str]:
    _non_numeric = {"network_ids", "costs", "newick"}
    return [c for c in fieldnames if c not in _non_numeric]


def main() -> None:
    here = Path(__file__).resolve().parent
    nets_path = here / "nonbinary_nets.csv"
    out_path = here / "large_network_results.csv"
    rng = np.random.default_rng(RNG_SEED)
    fieldnames = _fieldnames()
    numeric_cols = _numeric_fieldnames(fieldnames)

    # ── Numba warmup ─────────────────────────────────────────────────────────
    print("Warming up Numba JIT...", flush=True)
    _wid = f"{ID_LOW:05d}"
    _wnewick = _read_newick_by_id(_wid, nets_path)
    _wnet = DirectedPhyNetwork.from_string(_wnewick)
    _wte, _, _, _ = _compute_tree_extension(_wnet)
    _wtaxa = sorted(_wnet.taxa)
    _wcosts = {t: 1 for t in _wtaxa}
    _wtotal = sum(_wcosts.values())
    _wbudget = max(0, min(_wtotal, int(round(0.50 * _wtotal))))
    pp.all_paths.solve_maximization(
        _wnet, budget=_wbudget, costs=_wcosts,
        algorithm="nsw_fpt_budget", tree_extension=_wte, numba=True,
    )
    pp.max_tree.solve_maximization(
        _wnet, budget=_wbudget, costs=_wcosts, tree_extension=_wte, numba=True,
    )
    pp.min_tree.compute_diversity(_wnet, set(_wtaxa), tree_extension=_wte)
    print("Warmup done.\n", flush=True)

    all_rows: list[dict[str, Any]] = []

    with out_path.open("w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter=" ")
        writer.writeheader()

        for trial in range(N_TRIALS):
            sampled_ints = rng.integers(ID_LOW, ID_HIGH + 1, size=N_NETWORKS_PER_TRIAL)
            sampled_ids = [f"{x:05d}" for x in sampled_ints]
            compound_ids = [f"{i}_{nid}" for i, nid in enumerate(sampled_ids)]

            print(f"Trial {trial + 1}/{N_TRIALS}: networks {sampled_ids}", flush=True)

            indexed_networks: list[tuple[str, DirectedPhyNetwork]] = []
            for compound_id, nid in zip(compound_ids, sampled_ids):
                newick = _read_newick_by_id(nid, nets_path)
                net = DirectedPhyNetwork.from_string(newick)
                indexed_networks.append((compound_id, net))

            combined = _build_combined_network(indexed_networks)
            taxa = sorted(combined.taxa)
            costs = _assign_lognormal_costs(taxa, rng=rng)
            total_cost = sum(costs.values())

            # Compute level and newick before the timed section (level can be slow).
            level = net_level(combined)
            newick = combined.to_string()
            # Costs as compact "(taxon,cost)" pairs, same format as costs.csv.
            costs_str = "".join(f"({t},{c})" for t, c in sorted(costs.items()))

            print(
                f"  {combined.number_of_nodes()} nodes, "
                f"{combined.number_of_edges()} edges, "
                f"{len(taxa)} taxa, level {level}, total cost {total_cost}",
                flush=True,
            )

            # Compute tree extension in two separately timed steps.
            te, nsw_value, sec_nsw, sec_te = _compute_tree_extension(combined)
            print(f"  NSW value: {nsw_value}  ({sec_nsw:.3f}s)", flush=True)
            print(f"  Tree ext:  {sec_te:.3f}s", flush=True)

            # Three budget levels: 25%, 50%, 90% of total taxon cost (rounded).
            budgets = tuple(
                max(0, min(total_cost, int(round(f * total_cost))))
                for f in _BUDGET_FRACTIONS
            )

            row: dict[str, Any] = {
                "network_ids": "|".join(sampled_ids),
                "n_taxa": len(taxa),
                "tot_cost": total_cost,
                "level": level,
                "nsw_value": nsw_value,
                "time_nsw": f"{sec_nsw:.9f}",
                "time_te": f"{sec_te:.9f}",
            }

            # ── all-paths NSW FPT solver ──────────────────────────────────
            for tag, budget in zip(_BUDGET_TAGS, budgets):
                t0 = time.perf_counter()
                v, _ = pp.all_paths.solve_maximization(
                    combined, budget=budget, costs=costs,
                    algorithm="nsw_fpt_budget", tree_extension=te, numba=True,
                )
                dt = time.perf_counter() - t0
                row[f"time_ap_{tag}"] = f"{dt:.9f}"
                row[f"val_ap_{tag}"] = f"{float(v):.17g}"
                print(f"  ap_{tag}: val={v:.4f}, time={dt:.3f}s", flush=True)

            # ── MaxTree PD (upper bound via largest spanning tree) ────────
            for tag, budget in zip(_BUDGET_TAGS, budgets):
                t0 = time.perf_counter()
                v, _ = pp.max_tree.solve_maximization(
                    combined, budget=budget, costs=costs,
                    tree_extension=te, numba=True,
                )
                dt = time.perf_counter() - t0
                row[f"time_mxt_{tag}"] = f"{dt:.9f}"
                row[f"val_mxt_{tag}"] = f"{float(v):.17g}"
                print(f"  mxt_{tag}: val={v:.4f}, time={dt:.3f}s", flush=True)

            # ── MinTree PD (lower bound via smallest spanning tree) ───────
            t0 = time.perf_counter()
            v = pp.min_tree.compute_diversity(combined, set(taxa), tree_extension=te)
            dt = time.perf_counter() - t0
            row["time_mnt_full"] = f"{dt:.9f}"
            row["val_mnt_full"] = f"{float(v):.17g}"
            print(f"  mnt_full: val={v:.4f}, time={dt:.3f}s\n", flush=True)

            # Newick and costs go last; newick may contain spaces (tab-delimited
            # file keeps them unambiguous).
            row["costs"] = costs_str
            row["newick"] = newick

            writer.writerow(row)
            all_rows.append(row)

    # ── Print summary table ────────────────────────────────────────────────
    print(f"{'=' * 78}")
    print(f"SUMMARY STATISTICS  ({N_TRIALS} trials)")
    print(f"{'=' * 78}")
    print(f"{'Metric':<22} {'Mean':>13} {'Median':>13} {'Mode':>13} {'Std':>13}")
    print(f"{'-' * 78}")
    for col in numeric_cols:
        values = [float(r[col]) for r in all_rows]
        s = _summary_stats(values)
        print(
            f"{col:<22} {s['mean']:>13.4f} {s['median']:>13.4f} "
            f"{s['mode']:>13.4f} {s['std']:>13.4f}"
        )
    print(f"\nResults written to {out_path}")


if __name__ == "__main__":
    main()
