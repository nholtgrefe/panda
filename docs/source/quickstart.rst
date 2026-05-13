Manual
======

phypanda computes and optimizes phylogenetic diversity (PD) on directed
phylogenetic networks. This manual covers the full workflow: loading a network,
choosing a diversity measure, running the available methods, and understanding
the options for the different algorithms.

The standard import convention used throughout is:

.. code-block:: python

   import phypanda as pp

Loading a network
-----------------

phypanda works with :class:`~phylozoo.core.network.dnetwork.DirectedPhyNetwork`
objects from `phylozoo <https://github.com/nholtgrefe/phylozoo>`_. Networks can be
loaded from an eNewick string, from a file, or constructed from explicit edge and
node lists:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork

   # From an eNewick string — the standard serialization format
   network = DirectedPhyNetwork.from_string("((t1:1,t2:1)#H1:1,(#H1:1,t3:1):1);")

   print(f"Taxa:  {sorted(network.taxa)}")
   print(f"Nodes: {network.number_of_nodes()}")
   print(f"Edges: {network.number_of_edges()}")

The ``network.taxa`` attribute returns the leaf labels; ``network.leaves`` returns
the underlying node objects. See the `phylozoo documentation
<https://github.com/nholtgrefe/phylozoo>`_ for all I/O options, including loading from
files.

Choosing a measure
------------------

phypanda currently provides three diversity measures, each available as a module-level
singleton:

* **pp.all_paths** — all-paths diversity (MAPD). Sums the branch lengths of
  all edges used by at least one root-to-leaf path through the selected taxa.
  Both fixed-set scoring and budgeted optimization are implemented.
* **pp.max_tree** — max-tree PD. The maximum PD value over all displayed trees
  of the network, restricted to the selected taxa. Budgeted optimization is
  implemented; fixed-set scoring is not yet available.
* **pp.min_tree** — min-tree PD. The minimum PD value over all displayed trees.
  Fixed-set scoring is implemented; budgeted optimization is not yet available.

See below for details on the available methods for each measure.

All three implement the :class:`~phypanda.protocol.DiversityMeasure` protocol.
Stubs that currently raise ``NotImplementedError``exist for:

* **pp.tree** - tree PD. The above three measures generalize this standard measure on trees.
* **pp.network** - network PD.
* **pp.average_tree** - average-tree PD. The average PD value over all displayed trees.

Available methods
-----------------

Every measure singleton exposes two core methods, called directly on the singleton:

**Fixed-set scoring** — :meth:`compute_diversity <phypanda.protocol.DiversityMeasure.compute_diversity>`:

This method computes the diversity of a given set of taxa under the measure. The selected set must be a subset of the network's taxa; otherwise, a ``ValueError`` is raised.

.. code-block:: python

   value = pp.all_paths.compute_diversity(network, {"t1", "t2"})
   # or:
   value = pp.min_tree.compute_diversity(network, {"t1", "t2"})

**Budgeted optimization** — :meth:`solve_maximization <phypanda.protocol.DiversityMeasure.solve_maximization>`:

This method finds a subset of taxa that maximizes the diversity under the measure, subject to an integer budget constraint and non-negative integer taxon costs. The budget and all costs must be non-negative integers.

.. code-block:: python

   costs = {t: 1 for t in network.taxa}

   value, selected = pp.all_paths.solve_maximization(network, budget=2, costs=costs)
   # or:
   value, selected = pp.max_tree.solve_maximization(network, budget=2, costs=costs)

When ``costs=None``, unit costs (all ones) are used. The budget and all costs must
be non-negative integers.

:mod:`phypanda.base` additionally provides four generic helpers that accept any
measure as an argument:

* :func:`~phypanda.base.compute_diversity` — validated fixed-set score; checks that all
  taxa are present in the network before calling ``compute_diversity``.
* :func:`~phypanda.base.marginal_diversities` — computes the diversity gain (or
  loss) of adding (or removing) each taxon relative to a current selected set.
* :func:`~phypanda.base.greedy_max_diversity` — greedy budgeted heuristic; at
  each step picks the affordable taxon with the highest normalized gain
  (marginal gain / cost). Works with any measure, including those without an
  exact ``solve_maximization``.
* :func:`~phypanda.base.solve_max_diversity` — calls ``measure.solve_maximization``
  and raises ``NotImplementedError`` if the measure does not support it.

See :doc:`api/base` for signatures and examples.

Measures
--------

All-paths diversity
~~~~~~~~~~~~~~~~~~~

**Fixed-set scoring.**

.. code-block:: python

   value = pp.all_paths.compute_diversity(network, {"t1", "t2"})
   print(f"MAPD: {value:.4f}")

**Budgeted optimization — node-scanwidth FPT** (``"nsw_fpt_budget"``).

The default exact algorithm runs a DP along a tree extension of the network,
parameterized by the node scanwidth. Pass ``numba=True`` to JIT-accelerate
the merge step (slow on first call due to compilation; fast thereafter):

.. code-block:: python

   costs = {t: 1 for t in network.taxa}

   value, selected = pp.all_paths.solve_maximization(
       network,
       budget=2,
       costs=costs,
       algorithm="nsw_fpt_budget",  # default; can be omitted
       numba=True,                   # JIT-accelerate the DP (slow on first call)
   )
   print(f"MAPPD: {value:.4f}, selected: {selected}")

Pass ``numba=False`` to skip JIT compilation (useful for debugging).

**Budgeted optimization — edge-scanwidth FPT** (``"esw_fpt"``).

An alternative exact algorithm parameterized by the *edge* scanwidth. It requires
**unit costs** (all taxon costs equal) and a **binary** network (if the network is
non-binary, phypanda resolves it automatically):

.. code-block:: python

   value, selected = pp.all_paths.solve_maximization(
       network,
       budget=2,
       costs=costs,
       algorithm="esw_fpt",
   )

See :doc:`api/all_paths` for the full class, parameter reference, and low-level
solver functions.

Max-tree PD
~~~~~~~~~~~

Max-tree PD is the maximum PD value over all displayed trees of the network for the
selected taxa. It provides an **upper bound** on the true network PD.

**Budgeted optimization.** An exact algorithm using a node-scanwidth DP.
This solver requires a network **without parallel edges**:

.. code-block:: python

   costs = {t: 1 for t in network.taxa}

   value, selected = pp.max_tree.solve_maximization(
       network,
       budget=2,
       costs=costs,
       numba=True,    # default; set to False for pure-Python merges
   )
   print(f"MaxTPD: {value:.4f}, selected: {selected}")

A precomputed tree extension can be passed as ``tree_extension=te`` to skip
recomputation when solving the same network repeatedly — see :doc:`api/max_tree`
for details.

.. note::

   ``pp.max_tree.compute_diversity`` (fixed-set scoring) is not yet implemented.

See :doc:`api/max_tree` for the full class and parameter reference.

Min-tree PD
~~~~~~~~~~~

Min-tree PD is the minimum PD value over all displayed trees of the network for the
selected taxa. It provides a **lower bound** on the true network PD.

**Fixed-set scoring.** An exact DP along a tree extension of the network
(or its subnetwork induced by the selected taxa):

.. code-block:: python

   score = pp.min_tree.compute_diversity(network, {"t1", "t2"})
   print(f"MinTPD: {score:.4f}")

When the selected set is a strict subset of all taxa, the network is pruned to
those taxa internally before running the DP.

A precomputed tree extension can be passed as ``tree_extension=te``. When the
network is pruned for a subset, phypanda adapts the extension automatically:

.. code-block:: python

   from scanwidth import DAG, node_scanwidth

   dag = DAG(network._graph._graph)
   _, ext = node_scanwidth(dag)
   te = ext.to_canonical_tree_extension()

   # Extension is adapted internally for the pruned subnetwork
   score = pp.min_tree.compute_diversity(network, {"t1", "t2"}, tree_extension=te)

.. note::

   ``pp.min_tree.solve_maximization`` (budgeted optimization) is not yet
   implemented.

See :doc:`api/min_tree` for the full class and utility-function reference,
including :func:`~phypanda.measure.min_tree.prune_to_taxa` and
:func:`~phypanda.measure.min_tree.adapt_tree_extension_for_pruned_network`.

Quick Links
-----------

* :doc:`Installation Guide <installation>` — requirements and install options
* :doc:`API Reference <api/index>` — complete index of all functions and classes
