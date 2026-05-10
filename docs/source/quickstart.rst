Quickstart
==========

Import pattern
--------------

::

   import phypanda as pp

All concrete **measures** are objects implementing :class:`~phypanda.protocol.DiversityMeasure` (``compute_diversity`` for a fixed taxon set, and usually ``solve_maximization`` for budgeted selection). They are exposed as module-level singletons on :mod:`phypanda`.

Which PD measures are implemented?
----------------------------------

.. list-table::
   :widths: 18 22 20 40
   :header-rows: 1

   * - Measure (``pp.…``)
     - Role
     - Typical algorithm / backend
     - Status
   * - ``all_paths``
     - **All-paths** diversity (MAPD objective on the network).
     - Budgeted MAPPD: **node-scanwidth FPT** (``algorithm="nsw_fpt_budget"``); optional **edge-scanwidth FPT** (``algorithm="esw_fpt"``). Uses ``scanwidth`` for extensions unless you pass a :class:`~scanwidth.tree_extension.TreeExtension`.
     - **Implemented** (see :doc:`api/all_paths`, :doc:`api/all_paths_solvers`).
   * - ``max_tree``
     - **Max-tree PD**: maximize PD over **displayed trees** under a budget.
     - Node-scanwidth DP on a tree extension (Numba-accelerated merges by default).
     - **Implemented** (:doc:`api/max_tree`).
   * - ``min_tree``
     - **Min-tree PD**: minimize PD over displayed trees for a **fixed** taxon set.
     - Node-scanwidth DP; optional extension for the **full** network with automatic adaptation when inducing on a subset (:doc:`api/min_tree`).
     - **Implemented** (:doc:`api/min_tree`).
   * - ``tree``, ``network``, ``average_tree``
     - Placeholders for additional measures.
     - —
     - **Stubs** (raise ``NotImplementedError``; see :doc:`api/tree`, :doc:`api/network`, :doc:`api/average_tree`).

How a typical solve works
-------------------------

1. **Input**: a :class:`~phylozoo.core.network.dnetwork.DirectedPhyNetwork` ``network``, integer **budget** \(B\), and a **cost** per taxon (default: unit costs).
2. **Tree extension** (unless you supply one): for node-scanwidth solvers, ``scanwidth`` builds (or you provide) a **tree extension** Γ of the network DAG—valid for the DP order.
3. **DP**: the solver fills tables along Γ (and merges child states at reticulations). For MAPPD, choose ``algorithm="nsw_fpt_budget"`` or ``"esw_fpt"`` on :meth:`AllPathsDiversity.solve_maximization <phypanda.measure.all_paths.AllPathsDiversity.solve_maximization>`.
4. **Output**: objective value and a set of **selected taxon labels**.

**Parallel edges**: several routines require a network **without** parallel arcs; phylozoo may merge them in preprocessing.

High-level helpers
------------------

:mod:`phypanda.base` adds generic wrappers that dispatch to any measure:

* :func:`~phypanda.base.diversity` — fixed-set score.
* :func:`~phypanda.base.marginal_diversities` — per-taxon gain/loss vs. a current set.
* :func:`~phypanda.base.greedy_max_diversity` — greedy budgeted construction.
* :func:`~phypanda.base.solve_max_diversity` — calls ``measure.solve_maximization`` (exact when the measure implements it).

See :doc:`api/base` for signatures.

Examples
--------

**Fixed taxon set — all-paths diversity**

.. code-block:: python

   import phypanda as pp

   # network: phylozoo DirectedPhyNetwork
   value = pp.diversity(network, {"t1", "t2", "t3"}, measure=pp.all_paths)

**Budgeted MAPPD (default NSW FPT backend)**

.. code-block:: python

   import phypanda as pp

   costs = {t: 1 for t in network.taxa}
   value, taxa = pp.all_paths.solve_maximization(
       network,
       budget=50,
       costs=costs,
       algorithm="nsw_fpt_budget",
   )

**Budgeted max-tree PD**

.. code-block:: python

   import phypanda as pp

   value, taxa = pp.max_tree.solve_maximization(
       network,
       budget=50,
       costs={t: 1 for t in network.taxa},
   )

**Min-tree PD on a subset of taxa** (extension may be computed once on the full network)

.. code-block:: python

   import phypanda as pp
   from scanwidth import DAG, node_scanwidth

   dag = DAG(network._graph._graph)
   _, ext = node_scanwidth(dag)
   te = ext.to_canonical_tree_extension()

   score = pp.min_tree.compute_diversity(
       network,
       {"t1", "t2"},
       tree_extension=te,
   )

Next steps
----------

* Full parameter lists, exceptions, and solver kwargs: :doc:`api/index`.
* Supplementary scripts and datasets: ``experiments/`` in the `GitHub repo <https://github.com/nholtgrefe/panda>`_.
