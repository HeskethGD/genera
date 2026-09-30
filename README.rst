Genera
======

Genera is an early-stage Python package for arbitrary-precision computation
with algebraic curves and Abelian functions, particularly for applications in
integrable systems.

Genera uses `mpmath <https://mpmath.org/>`_ for arbitrary-precision arithmetic.
It is an independent project and does not modify the mpmath namespace.

The migrated public functions are available directly from Genera::

    from genera import kleinian_sigma, rtheta

The numerical implementation began in the mpmath ``algebraic-curve``
development branch. Genera retains the applicable BSD-3-Clause copyright and
licence notice while developing and releasing the higher-genus functionality
independently.

Development
-----------

Create an isolated environment and install the development dependencies::

    python -m venv .venv
    .venv/bin/python -m pip install -e '.[develop]'

Run the local checks with::

    .venv/bin/ruff check .
    .venv/bin/pytest
    .venv/bin/pytest --cov=genera --cov-report=term-missing
    .venv/bin/sphinx-build -W -b html docs build/sphinx/html
    .venv/bin/python -m build

The public numerical API will be added incrementally during migration from the
mpmath algebraic-curve development branch.
