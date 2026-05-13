API reference
=============

Docstrings are defined **once** per module below. Import everything from the top level as usual (``import phypanda as pp``); this page maps names to the detailed sections—there is **no** separate duplicate "package" autodoc page.

.. list-table:: Where to find what
   :widths: 35 65
   :header-rows: 1

   * - Topic
     - Page
   * - ``compute_diversity``, ``marginal_diversities``, ``greedy_max_diversity``, ``solve_max_diversity``
     - :doc:`base`
   * - ``DiversityMeasure`` protocol
     - :doc:`protocol`
   * - ``all_paths``, ``AllPathsDiversity``, ``solve_esw_fpt``, ``solve_nsw_fpt_budget``
     - :doc:`all_paths`
   * - ``max_tree``, ``MaxTreeDiversity``
     - :doc:`max_tree`
   * - ``min_tree``, ``MinTreeDiversity``, ``prune_to_taxa``, ``adapt_tree_extension_for_pruned_network``
     - :doc:`min_tree`
   * - ``tree``, ``network``, ``average_tree`` stubs
     - :doc:`tree`, :doc:`network`, :doc:`average_tree`

.. toctree::
   :maxdepth: 2
   :caption: Modules

   base
   protocol
   all_paths
   max_tree
   min_tree
   tree
   network
   average_tree
