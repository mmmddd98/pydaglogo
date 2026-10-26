pydaglogo documentation
=======================

``pydaglogo`` is a Python implementation of the core ``dagLogo`` workflow for
background-aware differential amino-acid logo analysis in proteomics-style
sequence windows.

.. image:: assets/hero_ztest_diff_logo.png
   :alt: pydaglogo ztest differential logo
   :width: 420

The package focuses on:

- differential amino-acid usage against an explicit background
- repaired ``ztest`` and ``fisher`` support
- exact validation against committed R ``dagLogo`` references for the public benchmark
- simple pipeline entrypoints via ``run_daglogo(...)`` and ``make_daglogo(...)``

Quickstart
----------

Install the package plus docs dependencies:

.. code-block:: bash

   pip install -e .[docs]

Build the HTML docs from the ``pydaglogo`` directory:

.. code-block:: bash

   sphinx-build -b html docs docs/_build/html

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Documentation

   usage
   api

