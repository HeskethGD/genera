Genera documentation
====================

Genera provides arbitrary-precision numerical algorithms for algebraic curves,
Abelian functions, and integrable systems. It uses mpmath for numerical types,
arithmetic, matrices, and precision contexts.

Genera's implementations receive an mpmath numerical context internally,
allowing them to respect the caller's working precision without modifying
mpmath's public namespace. Pass a custom context with the keyword-only
``ctx`` argument; otherwise Genera uses ``mpmath.mp``.

Algebraic curves
----------------

Construct a curve from ascending coefficients for
:math:`y^2 = \sum_k c_k x^k`::

   >>> from genera import algebraic_curve
   >>> curve = algebraic_curve((0, -1, 0, 1))
   >>> curve.genus
   1

General plane curves accept a sparse mapping from ``(x_power, y_power)`` to
coefficient. Computations are evaluated lazily through the :class:`genera.Curve`
interface.

.. autofunction:: genera.algebraic_curve

.. autoclass:: genera.Curve
   :members:

Riemann theta functions
-----------------------

.. autofunction:: genera.rtheta

.. autofunction:: genera.rtheta_jet

Kleinian functions
------------------

.. autofunction:: genera.kleinian_sigma

.. autofunction:: genera.kleinian_sigma_jet

.. autofunction:: genera.kleinian_zeta

.. autofunction:: genera.kleinian_p

.. autofunction:: genera.kleinian_baker_akhiezer

.. toctree::
   :maxdepth: 1

   references
