Ônishi's sigma and determinant identities
=========================================

This example follows Theorems 7.2 and 8.3 of Y. Ônishi,
*Determinant Expressions for Hyperelliptic Functions* [Onishi2005]_.
It compares two independent
calculations: a quotient of sigma functions at Abelian images, and a
determinant constructed directly from affine curve coordinates.

The runnable program supports arbitrary genus and determinant size.
This page selects genus two and four points, reaching the first
:math:`y` column of the determinant. It also checks the confluent,
Kiepert-type identity in both Abelian coordinates. No time integration
is involved.

The curve and coordinate conventions
------------------------------------

Ônishi writes the curve as :math:`y^2=f(x)` with monic :math:`f`.
For this demonstration,

.. math::

   f(x)=x(x^2-1)(x^2-4).

Genera uses the equivalent model

.. math::

   \begin{aligned}
   Y^2&=4f(x),\\
   Y&=2y.
   \end{aligned}

Its holomorphic basis is

.. math::

   \begin{aligned}
   du_1&=\frac{dx}{Y},\\
   du_2&=\frac{x\,dx}{Y}.
   \end{aligned}

These are exactly the paper's differentials with denominator :math:`2y`.
The Abel map therefore receives ``(x, Y)``, while the algebraic determinant
uses ``(x, y)``. Keeping that factor of two explicit is essential.

The period and characteristic data are obtained from the same curve::

   from genera import algebraic_curve, kleinian_sigma, kleinian_sigma_jet

   curve = algebraic_curve(coefficients)
   first = curve.periods_kind_1()
   second = curve.periods_kind_2()
   data = (
       first.omega, first.tau, second.kappa,
       curve.riemann_constant().characteristic,
   )

Here ``coefficients`` is the ascending coefficient vector of :math:`4f`.
For each selected point, the example computes its unreduced Abel image::

   image = curve.abel_map_kind_1((x, 2*y))

The images retain the integration paths from infinity. Sums and differences
are formed from those representatives without reducing each sigma argument
independently to a fundamental cell.

Sigma on the one-point stratum
------------------------------

The Abel image :math:`u` of a single curve point lies on :math:`\sigma(u)=0`.
The denominator in Ônishi's identity uses a special sigma derivative,
rather than sigma itself. In genus two,

.. math::

   \begin{aligned}
   \sigma_\sharp&=\sigma_2,\\
   \sigma_\flat&=\sigma.
   \end{aligned}

These correspond to derivative-count keys ``(0, 1)`` and ``(0, 0)`` in
``kleinian_sigma_jet``. A count tuple records how many times each Abelian
coordinate is differentiated; it differs from the index-list notation
used by ``kleinian_p``.

The complete example includes this source helper:

.. literalinclude:: ../../examples/onishi/onishi_determinant_demo.py
   :language: python
   :pyobject: sigma_derivative

Every sigma evaluation uses ``normalization="hyperelliptic"``. This
canonical normalization matters because the overall scale does not cancel
from the quotient.

The four-point identity
-----------------------

For four curve points :math:`P_i=(x_i,y_i)` with Abel images
:math:`u^{(i)}`, define

.. math::

   Q=\frac{
       \sigma\!\left(\sum_{i=1}^{4}u^{(i)}\right)
       \prod_{i<j}\sigma\!\left(u^{(i)}-u^{(j)}\right)
      }{\prod_{i=1}^{4}\sigma_2\!\left(u^{(i)}\right)^4}.

With the convention-adjusted sign derived by the example, Theorem 7.2 is

.. math::

   Q=\det M,

where

.. math::

   M=\begin{pmatrix}
   1&x_1&x_1^2&y_1\\
   1&x_2&x_2^2&y_2\\
   1&x_3&x_3^2&y_3\\
   1&x_4&x_4^2&y_4
   \end{pmatrix}.

The columns are ordered by pole order at infinity. In genus two their
orders are 0, 2, 4 and 5, so the :math:`y` column follows :math:`x^2`.
This exercises more than a Vandermonde determinant of powers of :math:`x`.

The first three x-coordinates are 4, 5 and 6, all on the positive
:math:`y` sheet. The fourth point
varies through :math:`7\leq x_4\leq9` on the same sheet. The dashed
curve is the algebraic determinant; markers are the independent sigma
quotient evaluated at 17 points. Numerical errors are described below.

.. plot::
   :include-source: True
   :context: reset

   from examples.onishi.onishi_determinant_demo import (
       documentation_example, make_figure,
   )

   result = documentation_example()
   make_figure(result)

Fixing the sign independently
------------------------------

The sign in a sigma identity depends on the normalization. The example
derives it before evaluating the full identity: it constructs the leading
Schur–Weierstrass polynomial with exact rational arithmetic and compares
the leading local sigma quotient with the leading determinant near
infinity. The resulting prefactor for this genus-two, four-point case
is :math:`+1`.

This is a local normalization calculation, rather than a sign chosen to
make the numerical values agree. The standalone program reports both the
paper's displayed sign and the sign in Genera's canonical convention.
The residual uses the difference of the two complex values, so an
unresolved sign or phase error would fail the check.

The confluent Kiepert identity
------------------------------

Theorem 8.3 replaces distinct point evaluations with derivatives at one
point. Here the three nonconstant determinant columns are

.. math::

   \begin{aligned}
   \phi_1&=x,\\
   \phi_2&=x^2,\\
   \phi_3&=y.
   \end{aligned}

Define the derivative along Abelian coordinate :math:`u_j` by

.. math::

   D_j=\frac{Y}{x^{j-1}}\frac{d}{dx}.

For :math:`n=4`, the factorial product is :math:`1!2!3!=12`.
The checked identity, with its independently derived prefactor :math:`+1`,
is

.. math::

   12\,\frac{\sigma(4u)}{\sigma_2(u)^{16}}
   =x^{6(j-1)}\det\left(D_j^r\phi_s\right)_{r,s=1}^{3}.

The example checks both :math:`j=1` and :math:`j=2` at the point with
:math:`x=4`. Repeated curve derivatives are evaluated algebraically,
using :math:`Y^2=4f(x)` to reduce expressions to :math:`A(x)+B(x)Y`.
This avoids an additional finite-difference approximation.

Executable validation
----------------------

The plot and doctest share one checked calculation at 30 decimal digits.
The algebraic curve is sampled densely at 129 fourth-point positions;
sigma values and Abel maps are needed at only 17 of them.

.. doctest::

   >>> from examples.onishi.onishi_determinant_demo import documentation_example
   >>> result = documentation_example()
   >>> len(result.sample_x), len(result.scaled_errors), len(result.kiepert_errors)
   (17, 17, 2)
   >>> result.prefactor, result.kiepert_prefactors
   (1, (1, 1))
   >>> max(result.scaled_errors + result.kiepert_errors) < 1e-24
   True
   >>> result.stratum_residual < 1e-24
   True

For each determinant comparison, the enforced scaled error is

.. math::

   \varepsilon=
   \frac{|Q_{\mathrm{sigma}}-Q_{\mathrm{algebraic}}|}
        {\max(1,|Q_{\mathrm{sigma}}|,|Q_{\mathrm{algebraic}}|)}.

This measures relative error for large values and absolute error near
zero. The helper raises an error if any four-point or Kiepert comparison
exceeds :math:`10^{-24}`. It also checks :math:`\sigma(u)=0` for every
one-point Abel image to the same absolute tolerance. The comparison covers
the sampled points and does not bound the identity residual between them.

A representative local run takes about 2.6 seconds before plot rendering.
The maximum four-point scaled error is approximately
:math:`1.3\times10^{-29}`. Both Kiepert comparisons have scaled errors
around :math:`4.1\times10^{-30}`, and the one-point sigma residual is
about :math:`1.4\times10^{-32}`.

Running the standalone example
-------------------------------

From the repository root, the corresponding four-point and confluent cases
can be run with::

   python -m examples.onishi.onishi_determinant_demo \
       --genus 2 --n 4 --kiepert-n 4 --coordinate 1 2 --dps 30

Omitting these options runs the broader genus-three demonstration. Higher
genera and larger determinants remain available, with increasing theta-series
and exact-polynomial costs. Matplotlib is imported only when requesting a
figure.
