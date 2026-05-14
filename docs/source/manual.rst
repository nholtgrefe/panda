Manual
======

phypanda computes and optimizes phylogenetic diversity (PD) on directed
phylogenetic networks. This manual covers the full workflow: loading a network,
computing diversity under different measures, and running the available
optimization algorithms.

The standard import convention used throughout is:

.. code-block:: python

   import phypanda as pp

Loading a network
-----------------

phypanda works with :class:`~phylozoo.core.network.dnetwork.DirectedPhyNetwork`
objects from `PhyloZoo <https://github.com/nholtgrefe/phylozoo>`_. Networks can be
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
the underlying node objects. See the `PhyloZoo documentation
<https://github.com/nholtgrefe/phylozoo>`_ for all I/O options, including loading
from files.

PD calculations
---------------

The four main entry points live in :mod:`phypanda.base` and each accept any
measure singleton (``pp.all_paths``, ``pp.max_tree``, ``pp.min_tree``) as the
``measure`` argument. See `PD Measures`_ below for the measure-specific options
and algorithms.

compute_diversity
~~~~~~~~~~~~~~~~~

:func:`~phypanda.base.compute_diversity` scores a fixed, user-supplied set of
taxa. The taxa must all be present in the network; a ``ValueError`` is raised
otherwise:

.. code-block:: python

   costs = {t: 1 for t in network.taxa}  # not needed for scoring, just for context

   # Score {"t1", "t2"} under all-paths diversity
   value = pp.compute_diversity(network, {"t1", "t2"}, measure=pp.all_paths)
   print(f"MAPPD: {value:.4f}")

   # Score the same set under min-tree PD
   value = pp.compute_diversity(network, {"t1", "t2"}, measure=pp.min_tree)
   print(f"MinTPD: {value:.4f}")

Measure-specific keyword arguments (such as ``tree_extension``) are forwarded
via ``**kwargs``. See `PD Measures`_ for details.

marginal_diversities
~~~~~~~~~~~~~~~~~~~~

:func:`~phypanda.base.marginal_diversities` computes, for every taxon in the
network, the diversity gain from adding it to the current set (or the loss from
removing it if it is already selected):

.. code-block:: python

   # Marginal gains/losses relative to the set {"t1"}
   marginals = pp.marginal_diversities(network, {"t1"}, measure=pp.all_paths)
   for taxon, delta in sorted(marginals.items()):
       print(f"  {taxon}: {delta:+.4f}")

Taxa already in the selected set receive a negative value (the loss from
removing them); taxa outside receive a positive value (the gain from adding
them).

greedy_max_diversity
~~~~~~~~~~~~~~~~~~~~

:func:`~phypanda.base.greedy_max_diversity` greedily builds a taxon subset that
maximizes diversity within a budget. At each step it picks the affordable taxon
with the highest normalized gain (marginal gain / cost). This heuristic works
with *any* measure, including those that do not have an exact solver:

.. code-block:: python

   costs = {t: 1 for t in network.taxa}  # unit costs

   value, selected = pp.greedy_max_diversity(
       network,
       budget=2,
       measure=pp.all_paths,
       costs=costs,
   )
   print(f"Greedy MAPPD: {value:.4f}, selected: {selected}")

When ``costs=None``, unit costs are used. The budget and all costs must be
non-negative integers.

solve_max_diversity
~~~~~~~~~~~~~~~~~~~

:func:`~phypanda.base.solve_max_diversity` calls the exact
``measure.solve_maximization`` algorithm and raises ``NotImplementedError`` if
the measure does not support budgeted optimization. Additional keyword arguments
are forwarded to the underlying solver (e.g. ``algorithm``, ``numba``,
``tree_extension``):

.. code-block:: python

   costs = {t: 1 for t in network.taxa}

   value, selected = pp.solve_max_diversity(
       network,
       budget=2,
       measure=pp.all_paths,  # or pp.max_tree
       costs=costs,
   )
   print(f"Optimal value: {value:.4f}, selected: {selected}")

See `PD Measures`_ below for the solver options available for each measure, and
:doc:`api/base` for the full parameter reference.

PD Measures
-----------

phypanda currently provides three diversity measures, each available as a
module-level singleton. All three implement the
:class:`~phypanda.protocol.DiversityMeasure` protocol.

Stubs that raise ``NotImplementedError`` exist for ``pp.tree``,
``pp.network``, and ``pp.average_tree``.

All-paths PD (MAPPD)
~~~~~~~~~~~~~~~~~~~~~~~~~~

All-paths PD (MAPPD) sums the branch lengths of all edges that lie on at
least one directed path from the root to a selected leaf. When the full taxon
set is selected this equals the total branch length of the network.

**Fixed-set scoring.**

.. code-block:: python

   value = pp.compute_diversity(network, {"t1", "t2"}, measure=pp.all_paths)
   print(f"MAPPD: {value:.4f}")

**Budgeted optimization — node-scanwidth FPT** (``"nsw_fpt_budget"``).

The default exact algorithm runs a dynamic program along a tree extension of the
network, parameterized by the node scanwidth. Pass ``numba=True`` to
JIT-accelerate the merge step (slow on the first call due to compilation; fast
thereafter):

.. code-block:: python

   costs = {t: 1 for t in network.taxa}

   value, selected = pp.solve_max_diversity(
       network,
       budget=2,
       measure=pp.all_paths,
       costs=costs,
       algorithm="nsw_fpt_budget",  # default; can be omitted
       numba=True,                   # JIT-accelerate the DP (slow on first call)
   )
   print(f"MAPPD: {value:.4f}, selected: {selected}")

Pass ``numba=False`` to skip JIT compilation (useful for debugging).

**Budgeted optimization — edge-scanwidth FPT** (``"esw_fpt"``).

An alternative exact algorithm parameterized by the *edge* scanwidth. It requires
**unit costs** (all taxon costs equal) and a **binary** network (phypanda
resolves non-binary networks automatically):

.. code-block:: python

   value, selected = pp.solve_max_diversity(
       network,
       budget=2,
       measure=pp.all_paths,
       costs=costs,
       algorithm="esw_fpt",
   )

See :doc:`api/all_paths` for the full parameter reference and low-level solver
functions.

Max-tree PD
~~~~~~~~~~~

Max-tree PD is the maximum of the standard PD over all displayed trees of the network
restricted to the selected taxa.

**Budgeted optimization.** An exact node-scanwidth DP. This solver requires a
network **without parallel edges**:

.. code-block:: python

   costs = {t: 1 for t in network.taxa}

   value, selected = pp.solve_max_diversity(
       network,
       budget=2,
       measure=pp.max_tree,
       costs=costs,
       numba=True,  # default; set to False for pure-Python merges
   )
   print(f"MaxTPD: {value:.4f}, selected: {selected}")

Pass ``tree_extension=te`` (a precomputed ``TreeExtension``) to skip
recomputation when solving the same network repeatedly.

.. note::

   Fixed-set scoring (``pp.compute_diversity(..., measure=pp.max_tree)``) is
   not yet implemented.

See :doc:`api/max_tree` for the full class and parameter reference.

Min-tree PD
~~~~~~~~~~~

Min-tree PD is the minimum of the standard PD over all displayed trees of the network
restricted to the selected taxa. It provides a **lower bound** on the true
network PD.

**Fixed-set scoring.** An exact DP along a tree extension. When the selected
set is a strict subset of all taxa, the network is pruned to those taxa
internally before running the DP:

.. code-block:: python

   score = pp.compute_diversity(network, {"t1", "t2"}, measure=pp.min_tree)
   print(f"MinTPD: {score:.4f}")

Pass ``tree_extension=te`` to reuse a precomputed extension; phypanda adapts
it automatically for the pruned subnetwork:

.. code-block:: python

   from scanwidth import DAG, node_scanwidth

   dag = DAG(network._graph._graph)
   _, ext = node_scanwidth(dag)
   te = ext.to_canonical_tree_extension()

   # Adapt and reuse the extension for a subset
   score = pp.compute_diversity(
       network, {"t1", "t2"}, measure=pp.min_tree, tree_extension=te
   )

.. note::

   Budgeted optimization (``pp.solve_max_diversity(..., measure=pp.min_tree)``)
   is not yet implemented.

See :doc:`api/min_tree` for the full class and utility-function reference.

Quick Links
-----------

* :doc:`Installation Guide <installation>` — requirements and install options
* :doc:`API Reference <api/index>` — complete index of all functions and classes