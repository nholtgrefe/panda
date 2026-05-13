Installation Guide
==================

Installing phypanda
--------------------

phypanda is a Python package that runs on `Python <https://www.python.org/>`_ (>= 3.7).
Choose one of:

* **From PyPI** — stable release:

  .. code-block:: bash

     pip install phypanda

* **From source** — editable install, recommended for development:

  .. code-block:: bash

     cd panda
     pip install -e .

Optional dependency groups
~~~~~~~~~~~~~~~~~~~~~~~~~~

Development and testing tools:

.. code-block:: bash

   pip install phypanda[dev]

Documentation dependencies:

.. code-block:: bash

   pip install phypanda[docs]

Requirements
^^^^^^^^^^^^

The following are required and installed automatically with pip:

* `phylozoo <https://github.com/nholtgrefe/phylozoo>`_ >= 0.1.2 — directed phylogenetic networks and I/O.
* `networkx <https://networkx.org/>`_ >= 3.0 — graph algorithms shared with phylozoo and scanwidth.
* `scanwidth <https://github.com/nholtgrefe/scanwidth>`_ >= 0.2.5 — tree extensions and node-scanwidth computation.
* `numba <https://numba.pydata.org/>`_ >= 0.56 — JIT-accelerated dynamic programming (can be disabled per call).

Verifying Installation
-----------------------

To verify that phypanda is installed correctly, import it and print the version.
The latest version is |version|.

.. code-block:: python

   >>> import phypanda as pp
   >>> print(pp.__version__)
   x.y.z  # your installed version

Building Documentation
-----------------------

To build the documentation locally, install the optional documentation dependencies:

.. code-block:: bash

   pip install -e ".[docs]"
   sphinx-build -b html docs/source docs/build/html

Open ``docs/build/html/index.html`` in a browser.

Troubleshooting
---------------

**Slow first call with** ``numba=True``: Numba compiles JIT kernels on first invocation.
This is expected — subsequent calls in the same Python session are fast. Pass
``numba=False`` to any solver call to skip JIT compilation entirely.

**Parallel edges**: Some solvers require a network without parallel arcs. phylozoo may
merge parallel edges during preprocessing. If you encounter unexpected solver errors,
check whether your network has parallel arcs.
