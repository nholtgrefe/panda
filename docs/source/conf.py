"""Sphinx configuration for phypanda documentation."""

from __future__ import annotations

import importlib.metadata

try:
    release = importlib.metadata.version("phypanda")
except importlib.metadata.PackageNotFoundError:
    release = "2.1.0"

project = "phypanda"
author = "Niels Holtgrefe"
copyright = "2025–2026, Niels Holtgrefe"
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
]

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}

templates_path = ["_templates"]
exclude_patterns: list[str] = []

html_theme = "pydata_sphinx_theme"
html_title = "PaNDA"
html_theme_options = {
    "logo": {
        "image_light": "_static/phypanda_icon.svg",
        "image_dark": "_static/phypanda_icon.svg",
        "text": "phypanda",
        "alt_text": "phypanda",
    },
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/nholtgrefe/panda",
            "icon": "fa-brands fa-github",
            "type": "fontawesome",
        },
    ],
}
html_static_path = ["_static"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "networkx": ("https://networkx.org/documentation/stable/", None),
}

nitpick_ignore = [
    ("py:class", "phylozoo.core.network.dnetwork.DirectedPhyNetwork"),
    ("py:class", "scanwidth.tree_extension.TreeExtension"),
    ("py:class", "scanwidth.dag.DAG"),
]

# :class:`DiversityMeasure` is documented both on ``phypanda`` and ``phypanda.protocol``.
suppress_warnings = ["ref.python"]
