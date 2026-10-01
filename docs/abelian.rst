Abelian functions
-----------------

Riemann theta functions generalize the Jacobi theta functions from one
complex variable to several. They arose in Riemann's theory of Abelian
functions and compact Riemann surfaces and are fundamental in complex
analysis, algebraic geometry, and integrable systems.

Theta functions are quasiperiodic. For integer vectors
:math:`m,n\in\mathbb Z^g`, the zero-characteristic function satisfies

.. math::

   \theta(z+m+\tau n\mid\tau)
   = \exp\!\left(-\pi i n^T\tau n-2\pi i n^Tz\right)\theta(z\mid\tau).

Suitable ratios and combinations are multiply-periodic Abelian functions. In
genus two these are meromorphic functions of two complex variables with a
rank-four period lattice. Kleinian sigma, zeta, and P-functions extend the
Weierstrass functions to higher genus.

Riemann theta functions
.......................

.. autofunction:: genera.rtheta

.. autofunction:: genera.rtheta_jet

Curve data for Kleinian functions
.................................

Curve construction, homology markings, period matrices, Riemann constants,
and Abel maps are documented in :doc:`algebraic_curves`. For a recognized
hyperelliptic curve, obtain mutually compatible inputs from one curve object::

   curve = algebraic_curve({(0, 2): 1, (1, 0): 4, (3, 0): -4})
   first = curve.periods_kind_1()
   second = curve.periods_kind_2()
   constant = curve.riemann_constant()
   omega = first.omega
   tau = first.tau
   kappa = second.kappa
   characteristic = constant.characteristic

The same curve supplies first- and second-kind Abelian integrals. Period data
should normally be constructed once and reused for many evaluations.

Kleinian functions
..................

The Kleinian functions use unnormalized Abelian coordinates ``u``. Their
period data are the first-kind half-period matrix ``omega``, the normalized
Riemann matrix ``tau``, and the symmetric matrix ``kappa``. The convention is

.. math::

   \tau=\omega^{-1}\omega',\qquad v=(2\omega)^{-1}u.

The characteristic uses the same literal convention as :func:`genera.rtheta`.
For compatible hyperelliptic curve data, the canonical normalization is
selected with ``normalization="hyperelliptic"``.

.. autofunction:: genera.kleinian_sigma

.. autofunction:: genera.kleinian_sigma_jet

.. autofunction:: genera.kleinian_zeta

.. autofunction:: genera.kleinian_p

The Baker--Akhiezer function remains available as an implementation-level
function in ``genera.kleinian`` while its public top-level export is being
stabilized. It is deliberately not part of the documented public API of this
release.

For the classical development of Abelian, theta, sigma, and multiply periodic
functions, see [BEL1997]_, [CEEK2000]_, and [Onishi2005]_.
