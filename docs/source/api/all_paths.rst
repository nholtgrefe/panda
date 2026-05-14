All-paths PD (MAPPD)
====================

.. automodule:: phypanda.measure.all_paths
   :members:
   :exclude-members: _all_paths_diversity_for_taxa

Solvers
-------

Public entry points for **budgeted all-paths maximization**. These are the functions registered on :class:`~phypanda.measure.all_paths.AllPathsDiversity` (``algorithm="nsw_fpt_budget"`` and ``algorithm="esw_fpt"``). Implementation details live in the submodules ``nsw_fpt_budget`` and ``esw_fpt``; they are **not** duplicated here to keep a single documentation target per function.

.. automodule:: phypanda.measure.all_paths_solvers
   :members:
