[![PyPI](https://img.shields.io/pypi/v/phypanda)](https://pypi.org/project/phypanda/)
[![License](https://img.shields.io/github/license/nholtgrefe/panda)](https://github.com/nholtgrefe/panda/blob/main/LICENSE)
[![Docs](https://img.shields.io/badge/docs-stable-blue)](https://nholtgrefe.github.io/panda/)

# phypanda

**phypanda** is the Python library for **PaNDA** — phylogenetic diversity optimization on **directed phylogenetic networks**. It builds on [phylozoo](https://pypi.org/project/phylozoo/) for networks and [scanwidth](https://pypi.org/project/scanwidth/) for node- and edge-scanwidth–based dynamic programs, including budgeted **maximum all-paths diversity (MAPD)** and related max-tree / min-tree objectives.

## Key features

- **All-paths (MAPD)**: fixed-set scoring and budgeted maximization via **node-scanwidth FPT** (`nsw_fpt_budget`) and **edge-scanwidth FPT** (`esw_fpt`), with optional precomputed tree extensions.
- **Max-tree and min-tree PD**: budgeted max-tree PD and fixed-set min-tree PD using node-scanwidth DPs (Numba-accelerated merges where applicable).
- **High-level API**: `diversity`, `marginal_diversities`, `greedy_max_diversity`, `solve_max_diversity` dispatching over measure implementations.
- **GUI**: optional desktop interface for interactive exploration — see the [`gui/`](https://github.com/nholtgrefe/panda/tree/main/gui) folder.

## Installation

Install the package from PyPI:

```bash
pip install phypanda
```

From a local clone (editable install):

```bash
pip install -e .
```

Optional extras:

```bash
# development (pytest)
pip install phypanda[dev]

# build the Sphinx documentation
pip install phypanda[docs]
```

Runtime dependencies (installed automatically) include `phylozoo`, `networkx`, `scanwidth`, and `numba`. Several solvers require networks **without parallel edges**; see the documentation.

## Documentation

Installation, quickstart, measure overview, and full API reference:

**[https://nholtgrefe.github.io/panda/](https://nholtgrefe.github.io/panda/)**

*(Site is deployed from this repository; if a page is not live yet, build locally with `pip install -e ".[docs]"` and `sphinx-build -b html docs/source docs/build/html` — see `docs/README.md`.)*

**Maintainers:** publishing to PyPI and deploying docs use separate GitHub Actions workflows on a **`v*`** tag (docs can also be run manually); see [`RELEASING.md`](RELEASING.md).

## Using PaNDA beyond the library

- **GUI**: [gui/README.md](https://github.com/nholtgrefe/panda/tree/main/gui) — screenshots and Tk-based workflow.
- **Example network**: Xiphophorus eNewick in [`experiments/MAPPD/exp2_xiphophorus_network.txt`](https://github.com/nholtgrefe/panda/blob/main/experiments/MAPPD/exp2_xiphophorus_network.txt).
- **Experiments & data**: [experiments/](https://github.com/nholtgrefe/panda/tree/main/experiments) — MAPPD simulation materials under [`experiments/MAPPD/`](https://github.com/nholtgrefe/panda/tree/main/experiments/MAPPD), budgeted PD benchmarks under [`experiments/budgeted_PD/`](https://github.com/nholtgrefe/panda/tree/main/experiments/budgeted_PD).

## Citation

If you use **phypanda** or this repository, please cite the main PaNDA paper:

> **PaNDA: Efficient Optimization of Phylogenetic Diversity in Networks.**  
> *Niels Holtgrefe, Leo van Iersel, Ruben Meuwese, Yukihiro Murakami, Jannik Schestag.*  
> bioRxiv, 2025. doi: [10.1101/2025.11.14.688467](https://www.biorxiv.org/content/10.1101/2025.11.14.688467)

Several **additional budgeted algorithms** (notably node-scanwidth FPT formulations for phylogenetic diversity on networks) are described in:

> **Tractable Optimization of Budgeted Phylogenetic Diversity on Networks using Node-Scanwidth.**  
> *Niels Holtgrefe and Jannik Schestag.*

Cite this second reference when your work builds specifically on those node-scanwidth optimization results.
