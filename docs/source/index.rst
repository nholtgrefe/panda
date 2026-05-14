
.. image:: _static/phypanda_full.svg
   :alt: phypanda logo
   :width: 400px
   :align: left

.. raw:: html

   <div style="clear: both;"></div>
   <br>

**Welcome to the phypanda docs!**

phypanda is the Python package corresponding to the PaNDA algorithm suite for computing and optimizing phylogenetic diversity (PD)
in directed phylogenetic networks. It works on networks represented by `PhyloZoo
<https://github.com/nholtgrefe/phylozoo>`_, uses the `scanwidth <https://github.com/nholtgrefe/scanwidth>`_ package
to compute tree extensions, and provides fixed-parameter tractable algorithms for several
budgeted PD problems. For definitions and proofs, see the references below.

You can score a set of taxa under different PD measures—all-paths (MAPPD), max-tree, and
min-tree—optimize under an integer budget with taxon costs, and supply precomputed tree
extensions to speed up repeated solves on the same network.

Documentation Overview
-----------------------

A good starting point is either the :doc:`Installation <installation>` page or the
:doc:`Manual <manual>` guide.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   manual
   api/index

Indices and tables
------------------

* :ref:`genindex` — index of all functions and classes.
* :ref:`modindex` — index of all modules.
* :ref:`search` — search the documentation.

Citation
--------

If you use this package in research, please cite:

   Niels Holtgrefe, Leo van Iersel, Ruben Meuwese, Yukihiro Murakami, and Jannik Schestag.
   *PaNDA: Efficient Optimization of Phylogenetic Diversity in Networks.*
   bioRxiv, 2025. doi:
   `10.1101/2025.11.14.688467 <https://doi.org/10.1101/2025.11.14.688467>`_.

If your work builds specifically on the budgeted node-scanwidth algorithms, please also cite:

   Niels Holtgrefe and Jannik Schestag.
   *Tractable Optimization of Budgeted Phylogenetic Diversity on Networks Utilizing Node-Scanwidth.*
   2026.
