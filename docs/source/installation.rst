Installation
============

Requirements
------------

* **Python** 3.7 or newer (see ``requires-python`` in ``pyproject.toml``).
* A normal scientific Python environment; heavy lifting for several solvers uses **Numba** (installed with the package).

Install the package
-------------------

**From PyPI** (when a release is published):

.. code-block:: bash

   pip install phypanda

**From a clone of the repository** (editable install recommended for development):

.. code-block:: bash

   cd panda
   python -m venv .venv
   source .venv/bin/activate   # Windows: .venv\Scripts\activate
   pip install -e .

Optional extras
---------------

**Run the test suite** (``pytest``):

.. code-block:: bash

   pip install -e ".[dev]"
   pytest

**Build this documentation** (Sphinx + PyData theme):

.. code-block:: bash

   pip install -e ".[docs]"
   sphinx-build -b html docs/source docs/build/html

Open ``docs/build/html/index.html`` in a browser.

Dependencies (automatic with ``pip install phypanda``)
-------------------------------------------------------

These are declared in ``pyproject.toml`` and pulled in by pip:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Package
     - Role in phypanda
   * - ``phylozoo``
     - Directed phylogenetic networks, I/O, subnetworks.
   * - ``networkx``
     - Graph algorithms shared with phylozoo / scanwidth.
   * - ``scanwidth``
     - Edge- and node-scanwidth, tree extensions, NSW / XP hooks.
   * - ``numba``
     - JIT-accelerated child merges inside several DPs (can be disabled where supported).

Related software
----------------

* **GUI**: optional desktop app in the ``gui/`` folder of the `PaNDA repository <https://github.com/nholtgrefe/panda>`_ (see its README for Tk / install notes).
