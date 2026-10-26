Usage
=====

Quickest example
----------------

If you already have cleavage-site residues in ``P4..P4'`` columns:

.. code-block:: python

   import pandas as pd

   from pydaglogo import make_daglogo, prepare_proteome

   frame = pd.read_csv("my_cleavage_sites.csv")
   proteome = prepare_proteome("/path/to/proteome.fasta")

   result, ax, matrix = make_daglogo(
       frame,
       proteome=proteome,
       input_mode="positions",
       test_type="ztest",
       labels=["P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'"],
   )

This assumes ``frame`` contains one residue column per position from ``P4``
through ``P4'``. Switch to ``test_type="fisher"`` if you want Fisher testing
instead.

Logo plots use ``Δ residue frequency`` as the default y-axis label. Pass
``y_label="..."`` to customize it, or ``y_label=None`` to suppress it.

Input styles
------------

``pydaglogo`` supports three practical input modes:

- ``input_mode="windows"`` for prebuilt sequence windows
- ``input_mode="positions"`` for one residue column per position
- ``input_mode="flanks"`` for upstream and downstream flank columns

Typical public examples live under ``examples/``:

- ``demo_caspsites_casp3_oneshot_ztest.py``
- ``demo_caspsites_casp3_oneshot_fisher.py``
- ``demo_flanks_oneshot_ztest.py``
- ``demo_flanks_oneshot_fisher.py``

Validation scope
----------------

The shipped public benchmark validates Python outputs against committed R
``dagLogo`` references for:

- motif counts
- background mean
- differential matrix
- statistics
- p-values
- final positive-only logo matrix

For Fisher testing, the package supports both:

- an exact R-backed parity path for the official benchmark
- a SciPy fallback for Python-only workflows

On the shipped public benchmark, the SciPy fallback is already extremely close
to the R reference, with only small residual differences in raw odds-ratio
estimates.

Custom grouping schemes
-----------------------

By default, ``pydaglogo`` operates on the full 20-amino-acid alphabet
(``grouping_scheme="no"``). You can collapse the alphabet into biochemical
groups by passing a :class:`~pydaglogo.GroupingScheme` object instead of a
string.

A :class:`~pydaglogo.GroupingScheme` has four fields:

- ``name`` — a short identifier for the scheme
- ``color`` — maps each group name to a hex color used in the logo
- ``symbols`` — maps each group name to a one-character plot label
- ``groups`` — maps each group name to a tuple of constituent amino-acid letters

Example: collapsing residues into four biochemical classes (hydrophobic, polar,
positively charged, negatively charged):

.. code-block:: python

   from pydaglogo import GroupingScheme, run_daglogo

   my_scheme = GroupingScheme(
       name="biochemical",
       color={
           "hydrophobic": "#FF6600",
           "polar":       "#0066FF",
           "positive":    "#CC0000",
           "negative":    "#009900",
       },
       symbols={
           "hydrophobic": "H",
           "polar":       "P",
           "positive":    "+",
           "negative":    "-",
       },
       groups={
           "hydrophobic": ("A", "F", "G", "I", "L", "M", "P", "V", "W"),
           "polar":       ("C", "N", "Q", "S", "T", "Y"),
           "positive":    ("H", "K", "R"),
           "negative":    ("D", "E"),
       },
   )

   result = run_daglogo(
       input_data,
       proteome=proteome,
       test_type="ztest",
       grouping_scheme=my_scheme,   # pass the object directly
   )

The ``grouping_scheme`` parameter is accepted by :func:`~pydaglogo.run_daglogo`,
:func:`~pydaglogo.make_daglogo`, and :func:`~pydaglogo.test_dau`. Passing the
object directly bypasses the name-lookup and lets you define any grouping
without modifying the package source.
