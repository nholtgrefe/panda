# Budgeted PD benchmarks

Materials for timing and validating **budgeted** phylogenetic diversity on directed networks (MAPD via NSW FPT, max-tree PD, min-tree PD), using the same contracted simulation corpus as in the broader PaNDA experiments.

## Data files (space-separated CSV unless noted)

- **`nonbinary_nets.csv`** — Network id, leaf count, level, eNewick (after edge contraction from the exp1 corpus).
- **`extension.csv`** — Per-network NSW / XP ordering used to build tree extensions.
- **`costs.csv`** — Integer taxon costs for budget experiments.
- **`algorithm_benchmark_results.csv`** — Wide timing table produced by the benchmark script (regenerate after code or data changes).

## Scripts

- **`algorithm_benchmark.py`** — Loads each network, runs MAPPD (three budget fractions), max-tree (same budgets), and full-set min-tree; writes `algorithm_benchmark_results.csv` next to this file.
- **`contract_edges.py`** — One-off pipeline to build `nonbinary_nets.csv` from `../MAPPD/exp1_simulated_networks.csv` (10% shortest-edge contraction).
- **`cost_assignment.py`**, **`node_scanwidth.py`** — Auxiliary utilities for costs and scanwidth on this corpus.

### Example

From the repository root:

```bash
PYTHONPATH=src python experiments/budgeted_PD/algorithm_benchmark.py
```

Plots that consume these tables (e.g. under `sandbox/plot_sims/`) resolve this folder automatically via `get_default_data_dir()`.

## Regenerating outputs

After editing solvers or input CSVs, re-run `algorithm_benchmark.py` so column names and timings stay consistent with the current `phypanda` API.
