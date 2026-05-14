[![PyPI](https://img.shields.io/pypi/v/phypanda)](https://pypi.org/project/phypanda/)
[![License](https://img.shields.io/github/license/nholtgrefe/panda)](https://github.com/nholtgrefe/panda/blob/main/LICENSE)
[![Docs](https://img.shields.io/badge/docs-stable-blue)](https://nholtgrefe.github.io/panda/)

# phypanda

<img src="docs/source/_static/phypanda_full.svg" alt="phypanda logo" width="200" align="right">

phypanda is the Python library for **PaNDA**—an algorithm library for phylogenetic diversity (PD)
optimization on directed phylogenetic networks. It builds on
[phylozoo](https://github.com/nholtgrefe/phylozoo) for network representations and
[scanwidth](https://github.com/nholtgrefe/scanwidth) for node- and edge-scanwidth–based
dynamic programming.

<br>

## Key Features

- **All-paths (MAPPD), Max-tree and Min-tree PD**: compute the diversity score of a fixed taxon set for several PD measures, or maximize PD using exact FPT algorithms parameterized by node scanwidth or edge scanwidth.
- **Budgeted maximization**: assign integer costs to taxa and maximize PD under budget constraints.
- **JIT-compilation**: speed up algorithms by optional `numba` JIT compilation.
- **High-level API** — `compute_diversity`, `marginal_diversities`, `greedy_max_diversity`, and `solve_max_diversity` work with any measure and dispatch to the appropriate solver.

## Installation

```bash
pip install phypanda
```

Runtime dependencies (`phylozoo`, `networkx`, `scanwidth`, `numba`) are installed automatically. For development or documentation extras:

```bash
pip install phypanda[dev]   # testing
pip install phypanda[docs]  # Sphinx documentation
```

## Documentation

For the full manual, API reference, and installation guide, visit the **[phypanda docs](https://nholtgrefe.github.io/panda/)**.

## Citation

If you use phypanda, please cite:

> Niels Holtgrefe, Leo van Iersel, Ruben Meuwese, Yukihiro Murakami, and Jannik Schestag.
> **PaNDA: Efficient Optimization of Phylogenetic Diversity in Networks.**
> *bioRxiv*, 2025. doi: [10.1101/2025.11.14.688467](https://doi.org/10.1101/2025.11.14.688467)

If your work builds specifically on the budgeted node-scanwidth algorithms, please also cite:

> Niels Holtgrefe and Jannik Schestag.
> **Tractable Optimization of Budgeted Phylogenetic Diversity on Networks Utilizing Node-Scanwidth.**
> 2026.


## See also

For the graphical-user interface developed for the first paper, please go to [`gui/`](https://github.com/nholtgrefe/panda/tree/main/gui). 

For the experimental materials corresponding to the above two papers, please go to [`experiments/`](https://github.com/nholtgrefe/panda/tree/main/experiments).