# Experiments

This directory holds **reproducibility materials** tied to publications and larger empirical runs: datasets, helper scripts, and benchmark outputs.

## Layout

| Path | Contents |
|------|----------|
| [`MAPPD/`](MAPPD/) | Supplementary material for the first **PaNDA** paper (MAPD simulations, Xiphophorus network, R/Python drivers). See [`MAPPD/README.md`](MAPPD/README.md). |
| [`budgeted_PD/`](budgeted_PD/) | **Budgeted** PD benchmarks: contracted networks, costs, NSW orderings, and timing scripts for all-paths (NSW FPT), max-tree, and min-tree measures. See [`budgeted_PD/README.md`](budgeted_PD/README.md). |

The core library implementation remains under `src/phypanda/`. Run scripts with `PYTHONPATH` pointing at `src` (or an installed `phypanda` package) as described in each subfolder’s README.
