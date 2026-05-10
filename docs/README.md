# phypanda documentation

Sphinx site for the **phypanda** Python package (PaNDA — phylogenetic diversity in networks).

## Structure

| Path | Purpose |
|------|---------|
| `source/index.rst` | Main landing page: scope, paper citation, repo pointers |
| `source/installation.rst` | Pip installs, extras (`dev`, `docs`), dependencies |
| `source/quickstart.rst` | Measures table, how solvers run, examples → API for details |
| `source/api/` | One autodoc page per module (**no duplicate** top-level package page) |

## Build

From the repository root:

```bash
pip install -e ".[docs]"
sphinx-build -b html docs/source docs/build/html
```

Open `docs/build/html/index.html`. Use `-W` locally to treat warnings as errors:

```bash
sphinx-build -b html docs/source docs/build/html -W
```
