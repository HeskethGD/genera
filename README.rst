Genera
======

Genera is an early-stage Python package for arbitrary-precision computation
with algebraic curves and Abelian functions, particularly for applications in
integrable systems.

Genera uses `mpmath <https://mpmath.org/>`_ for arbitrary-precision arithmetic.
It is an independent project and does not modify the mpmath namespace.

The migrated public interface is available directly from Genera::

    from genera import algebraic_curve, kleinian_sigma, rtheta

    curve = algebraic_curve((0, -1, 0, 1))
    print(curve.genus)

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

Further numerical APIs and documentation will be added incrementally during
migration from the mpmath algebraic-curve development branch.
