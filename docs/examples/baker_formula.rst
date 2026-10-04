Baker's genus-two addition formula
========================================

Baker's addition formula expresses a sigma ratio using only two-index
Kleinian P functions [Baker1907]_, Chapter V, p. 100. In modern notation
it reads as follows; see also [BEH2005]_, equation (3.9):

.. math::

   \frac{\sigma(u+v)\,\sigma(u-v)}{\sigma(u)^2\,\sigma(v)^2}
   = \wp_{22}(u)\wp_{12}(v)-\wp_{12}(u)\wp_{22}(v)
     +\wp_{11}(v)-\wp_{11}(u).

Here :math:`u,v\in\mathbb{C}^2` are arbitrary Abelian argument vectors,
with :math:`\sigma(u)\sigma(v)\ne0`. They do not need to be constructed
from curve points using Abel maps. This particular formula is for genus two;
higher genera have different addition formulae.

A compact calculation
---------------------

We use the smooth hyperelliptic curve
:math:`y^2=4x(x^2-1)(x^2-4)`. Its automatic holomorphic basis is
:math:`(dx/y,x\,dx/y)`. Passing the curve directly to ``ks`` and ``kp``
provides the periods, Riemann characteristic and, for sigma, the canonical
hyperelliptic normalization. That normalization matters because a constant
multiplier of sigma changes the ratio on the left.

Baker numbers the coordinates from one; Genera numbers them from zero.
Thus :math:`\wp_{11},\wp_{12},\wp_{22}` correspond to ``(0, 0)``,
``(0, 1)``, ``(1, 1)`` respectively.

.. literalinclude:: ../../examples/baker_formula/baker_formula_demo.py
   :language: python
   :start-at: from mpmath import mp
   :end-at: return lhs, rhs, residual

The same calculation is checked by the documentation:

.. doctest::

   >>> from examples.baker_formula.baker_formula_demo import documentation_example
   >>> lhs, rhs, residual = documentation_example()
   >>> abs(residual) < 1e-24
   True

At 30 decimal digits, a representative run gives both sides approximately
:math:`121.270036511` and an absolute residual
:math:`|\mathrm{RHS}-\mathrm{LHS}|\approx6.3\times10^{-29}`.
The script enforces an absolute tolerance of :math:`10^{-24}` for this pair
of vectors; this is a numerical check of the identity at the selected arguments.

Run it from the repository root with::

   python -m examples.baker_formula.baker_formula_demo
