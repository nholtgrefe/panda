phypanda
========

**phypanda** is the Python library for **PaNDA** — *efficient optimization of phylogenetic diversity in phylogenetic networks*. It works on **directed** phylogenetic networks represented in `phylozoo <https://pypi.org/project/phylozoo/>`_, and uses **scanwidth** (`scanwidth <https://pypi.org/project/scanwidth/>`_) to obtain tree extensions and run fixed-parameter tractable dynamic programs for several budgeted problems.

What you can do here
--------------------

* **Score** a chosen set of taxa under different PD notions (**all-paths** / MAPPD, **max-tree**, **min-tree**, and stubs for tree-only, network, and average-tree measures).
* **Optimize** under an integer **budget** with taxon **costs** (notably MAPPD via node-scanwidth FPT, max-tree PD, and greedy / generic wrappers).
* Provide a **precomputed tree extension** (e.g. from NSW / XP ordering) so repeated solves on the same network avoid rescanning.

The algorithms and objective definitions are those in the PaNDA paper; this site focuses on **how to run the code**. For definitions and proofs, see the reference below.

Where to read next
------------------

.. toctree::
   :maxdepth: 2
   :caption: User guide

   installation
   quickstart
   api/index

Paper and citation
------------------

If you use this software, please cite:

   **PaNDA: Efficient Optimization of Phylogenetic Diversity in Networks.**  
   *Niels Holtgrefe, Leo van Iersel, Ruben Meuwese, Yukihiro Murakami, Jannik Schestag.*  
   bioRxiv, 2025.  
   https://doi.org/10.1101/2025.11.14.688467

The `repository README <https://github.com/nholtgrefe/panda>`_ links the **GUI**, **experiment scripts** (under ``experiments/``), and supplementary datasets for the first MAPPD study (``experiments/MAPPD/``) and budgeted benchmarks (``experiments/budgeted_PD/``).
