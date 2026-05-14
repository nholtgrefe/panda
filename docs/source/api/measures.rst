Measures
========

phypanda provides three diversity measures, each available as a module-level
singleton (``pp.all_paths``, ``pp.max_tree``, ``pp.min_tree``). All implement
the :class:`~phypanda.protocol.DiversityMeasure` protocol. Stubs that raise
``NotImplementedError`` exist for ``pp.tree``, ``pp.network``, and
``pp.average_tree``.

.. list-table::
   :widths: 22 38 20 20
   :header-rows: 1

   * - Measure
     - Description
     - Fixed-set scoring
     - Budgeted optimization
   * - :doc:`all_paths`
     - All-paths diversity (MAPD) — edges on any root-to-leaf path
     - yes
     - yes (``nsw_fpt_budget``, ``esw_fpt``)
   * - :doc:`max_tree`
     - Max PD over all displayed trees — upper bound on network PD
     - not yet
     - yes (``nsw_fpt``)
   * - :doc:`min_tree`
     - Min PD over all displayed trees — lower bound on network PD
     - yes
     - not yet
   * - :doc:`tree`, :doc:`network`, :doc:`average_tree`
     - Standard tree PD and other objectives (stubs)
     - not yet
     - not yet

.. toctree::
   :maxdepth: 1

   all_paths
   max_tree
   min_tree
   tree
   network
   average_tree
