Changelog
=========

Version History
---------------

Version 2.1.2
~~~~~~~~~~~~~

* Fixed a bug where edges with branch length 0 were treated as edges of length 1
  by the optimization algorithms. Only a missing branch length should default to 1.
  This affected the all-paths solvers (``nsw_fpt_budget`` and ``esw_fpt``), the
  MaxTreePD solver, and MinTreePD (``min_tree.compute_diversity``). For networks
  with zero-length edges (for example, reticulation edges of length 0), these
  returned too high diversity values, and possibly a different optimal taxon set.
  ``all_paths.compute_diversity`` was not affected.
* Added the helper ``phypanda.utils.branch_length``, used by all solvers to read
  branch lengths.
* Updated the regression baselines (checked against brute-force computations) and
  added regression tests for zero-length edges.

Version 2.1.1
~~~~~~~~~~~~~

* Updated the citation of the PaNDA paper to the published version (Journal of
  Computational Biology).

Version 2.1.0
~~~~~~~~~~~~~

* Unified the interface for scoring a set of taxa in ``compute_diversity``.
* Faster child-merge dynamic programming, with optional Numba acceleration.
* Revised documentation.

Version 2.0.0
~~~~~~~~~~~~~

* Complete rewrite on top of `PhyloZoo <https://github.com/nholtgrefe/phylozoo>`_
  and the `scanwidth <https://github.com/nholtgrefe/scanwidth>`_ package.
* Added the budgeted node-scanwidth algorithm for all-paths diversity (taxon costs),
  and the MaxTreePD and MinTreePD measures.

Version 1.0.3
~~~~~~~~~~~~~

* Last release of the original PaNDA implementation (unit-cost MAPPD with the
  edge-scanwidth algorithm).
