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

The :ref:`Bobenko theta calculation <bobenko-theta-example>` uses ``rtheta``
for characteristic-dependent theta quotients and ``rtheta_jet`` for their
derivatives along the flow.

.. autofunction:: genera.rtheta

The following plots show two real slices and the modulus over two real
variables for genus-two period matrices. Similar slices and surfaces are
illustrated in `DLMF section 21.4 <https://dlmf.nist.gov/21.4>`_.

.. plot::

   import matplotlib.pyplot as plt
   from mpmath import j, plot, re
   from genera import rtheta

   tau = [[j, -0.5], [-0.5, j]]
   curves = [
       lambda x: re(rtheta([x, x/2], tau)),
       lambda x: re(rtheta([x, 2*x], tau)),
   ]
   fig, ax = plt.subplots()
   plot(curves, [-2, 2], axes=ax)
   # mpmath supplies its own styles; distinguish slices without colour.
   for line, style in zip(ax.lines, ("-", "--")):
       line.set_color("#015758")
       line.set_linestyle(style)
   ax.legend([r"$z=(x,x/2)$", r"$z=(x,2x)$"])

.. plot::

   import matplotlib.pyplot as plt
   from mpmath import j, splot
   from genera import rtheta

   tau = [[j, 0.5], [0.5, j]]
   fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
   surface = lambda x, y: abs(rtheta([x, y], tau))
   splot(surface, [-1, 1], [-1, 1], points=35, keep_aspect=False,
         axes=ax, plot3d_kwargs={"color": "#015758"})
   ax.set_zlabel(r"$|\theta(z\mid\tau)|$")

.. autofunction:: genera.rtheta_jet

Curve data for Kleinian functions
.................................

Curve construction, homology markings, period matrices, Riemann constants,
and Abel maps are documented in :doc:`algebraic_curves`. For a recognized
hyperelliptic curve, one curve object provides mutually compatible inputs::

   curve = Curve({(0, 2): 1, (1, 0): 4, (3, 0): -4})
   first = curve.periods_kind_1()
   second = curve.periods_kind_2()
   constant = curve.riemann_constant()
   omega = first.omega
   tau = first.tau
   kappa = second.kappa
   characteristic = constant.characteristic

The same curve supplies first- and second-kind Abelian integrals. Period data
should normally be constructed once and reused for many evaluations.
The optional ``ctx`` argument uses the same convention throughout Genera;
see the :ref:`Curve constructor examples <numerical-contexts>` for an
independent-precision example.

The :ref:`custom-differential-bases` guide covers custom first-kind sequences
and compatible second-kind forms. These are bound to ``Curve`` at construction;
its period, Riemann-constant, and Abel-map methods use the same bases automatically.
The Abel vector, ``omega``, ``kappa``, and P-function indices all refer to
that ordered basis. A change of basis or second-kind convention changes the
meaning of these inputs. In particular, the hyperelliptic sigma normalization
assumes the automatic differential ordering and scaling; it is not a general
normalizer for arbitrary custom bases.

The :ref:`Manakov example <manakov-period-data>` constructs and reuses these
inputs for a hyperelliptic curve; the :ref:`Bernatska example <bernatska-custom-bases>`
uses explicit first- and second-kind differential bases on a trigonal curve.


Kleinian functions
..................

The Kleinian functions use unnormalized Abelian coordinates ``u``. Their
period data are the first-kind half-period matrix ``omega``, the normalized
Riemann matrix ``tau``, and the symmetric matrix ``kappa``. The convention is

.. math::

   \begin{aligned}
   \tau&=\omega^{-1}\omega',\\
   v&=(2\omega)^{-1}u.
   \end{aligned}

The characteristic uses the same literal convention as :func:`genera.rtheta`.
For compatible hyperelliptic curve data, the canonical normalization is
selected with ``normalization="hyperelliptic"``.

The :ref:`Ônishi sigma calculation <onishi-sigma-example>` uses
``kleinian_sigma`` and ``kleinian_sigma_jet`` on Abel images, including
derivatives where sigma itself vanishes. The
:ref:`Manakov solution <manakov-kleinian-example>` uses batched
``kleinian_p`` evaluations to reconstruct amplitudes and momenta.

.. autofunction:: genera.kleinian_sigma

.. autofunction:: genera.kleinian_sigma_jet

.. autofunction:: genera.kleinian_zeta

.. autofunction:: genera.kleinian_p

For the classical development of Abelian, theta, sigma, and multiply periodic
functions, see [BEL1997]_, [CEEK2000]_, and [Onishi2005]_.
