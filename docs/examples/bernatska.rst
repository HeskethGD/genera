Bernatska's genus-three trigonal curve
======================================

This example follows Examples 3 and 3a--3b of J. Bernatska,
*Computation of P-Functions on Plane Algebraic Curves* [Bernatska2024]_.
Starting with a non-hyperelliptic
plane curve, Genera computes its branch locus, monodromy, homology and
periods, then evaluates the second- and third-order Kleinian P-functions
printed in the paper.

The comparison illustrates two conventions that matter in practical use:
the paper and Genera can choose different homology bases, and Bernatska's
full periods must be converted to the half-period inputs of
:func:`genera.kleinian_p`.

The curve and its branch locus
-------------------------------

Bernatska's equation (73) is the trigonal :math:`(3,4)` curve

.. math::

   f(x,y)=-y^3+x^4+3x^3+7x^2+16x+9+y(4x^2+5x+11)=0.

Genera accepts a sparse mapping from ``(x_power, y_power)`` to the
coefficient. The mapping below is included directly from the runnable
script:

.. literalinclude:: ../../examples/bernatska/bernatska_trigonal_demo.py
   :language: python
   :pyobject: trigonal_curve

The projection to the :math:`x`-plane is a three-sheeted covering. Write

.. math::

   \begin{aligned}
   P(x)&=x^4+3x^3+7x^2+16x+9,\\
   Q(x)&=4x^2+5x+11.
   \end{aligned}

Its finite branch values satisfy

.. math::

   27P(x)^2-4Q(x)^3=0.

Genera finds these values from the curve equation; the paper's rounded
branch values are used only for comparison. The result consists of two
real points and three complex-conjugate pairs. The figure shows the
computed branch values in the projection plane, rather than a choice of
branch cuts or a drawing of the three sheets.

.. plot::
   :include-source: True
   :context: reset

   from examples.bernatska.bernatska_trigonal_demo import (
       documentation_example, make_figure,
   )

   result = documentation_example()
   make_figure(result)

Infinity is also ramified and is not shown in this finite-plane figure.
Sheet continuation around the finite values gives eight transpositions;
continuation around infinity gives a three-cycle. The resulting homology
has genus three and intersection rank six, giving three canonical pairs
of cycles.

Supplying the differential bases
---------------------------------

For this general plane curve, the caller supplies the holomorphic and
second-kind differentials. A differential :math:`h(x,y)\,dx/f_y` is passed
as a callable returning its coefficient of :math:`dx`, with

.. math::

   f_y=-3y^2+4x^2+5x+11.

The first-kind basis is

.. math::

   \begin{aligned}
   du_1&=\frac{y\,dx}{f_y},\\
   du_2&=\frac{x\,dx}{f_y},\\
   du_3&=\frac{dx}{f_y}.
   \end{aligned}

The supplied second-kind basis is

.. math::

   \begin{aligned}
   dr_1&=\frac{x^2\,dx}{f_y},\\
   dr_2&=\frac{2xy\,dx}{f_y},\\
   dr_3&=\frac{R_5(x,y)\,dx}{f_y},
   \end{aligned}

where

.. math::

   R_5(x,y)=5x^2y+9xy+\frac{32}{3}x^2+7y+\frac{40}{3}x.

The construction uses the mpmath context explicitly::

   from genera import Curve
   from mpmath import mp

   curve = Curve(mp, trigonal_curve())
   first_kind, second_kind = differentials()
   first = curve.periods_kind_1(first_kind)
   second = curve.periods_kind_2(
       differentials=first_kind,
       second_differentials=second_kind,
   )

Here ``trigonal_curve`` and ``differentials`` are helpers in the complete
example. The returned records contain mutually compatible half-periods,
normalized periods and second-kind data. In Genera's notation,

.. math::

   \begin{aligned}
   \tau&=\omega^{-1}\omega',\\
   \varkappa&=\eta\omega^{-1}.
   \end{aligned}

Comparing the homology bases
-----------------------------

Period matrices depend on the chosen cycles. A direct entrywise comparison
with the paper would therefore be misleading before matching the bases.
The script compares the computed first-kind full periods with the rounded
paper matrices and resolves their change of basis to an integer matrix
:math:`C`. It then checks

.. math::

   \begin{aligned}
   |\det C|&=1,\\
   C^TJC&=J,
   \end{aligned}

where :math:`J` is the standard intersection matrix. These checks ensure
that the transformation preserves the full symplectic period lattice.
The same matrix is applied to the independently integrated second-kind
periods; their agreement provides an additional check of the basis change.

Let :math:`\Omega,\Omega'` and :math:`E,E'` denote the transformed full
period matrices in Bernatska's convention. The script forms

.. math::

   \begin{aligned}
   \tau&=\Omega^{-1}\Omega',\\
   \varkappa_{\mathrm{paper}}&=E\Omega^{-1}.
   \end{aligned}

It checks symmetry before removing the small antisymmetric quadrature
residuals, and verifies that :math:`\operatorname{Im}\tau` is positive
definite. These consistency checks are independent of the precision of
the printed period tables.

Evaluating the P-functions
--------------------------

For Example 3a, the paper supplies the Abelian vector

.. math::

   u=\begin{pmatrix}
   -0.270333-1.612257i\\
   -1.116879+0.562199i\\
   0.258194+0.268653i
   \end{pmatrix}.

The literal theta characteristic used by Genera is

.. math::

   \begin{aligned}
   a&=(\tfrac12,0,\tfrac12),\\
   b&=(0,\tfrac12,\tfrac12).
   \end{aligned}

Eight second- and third-order P-functions are evaluated in one batched
call. The runnable source shows both required convention conversions:

.. literalinclude:: ../../examples/bernatska/bernatska_trigonal_demo.py
   :language: python
   :start-at:     indices = (
   :end-before:     expected = (
   :dedent: 4

In this excerpt, ``omega`` denotes Bernatska's full matrix
:math:`\Omega`, and ``kappa`` denotes :math:`\varkappa_{\mathrm{paper}}`.
Dividing ``omega`` by two supplies half-periods. Negating ``kappa`` matches
Genera's sigma convention. The integer indices are zero-based: for example,
``(0, 0)`` denotes the paper's :math:`\wp_{11}`.

The calculation is repeated with the different Abelian vector of Example
3b, checking another eight P-function values from equation (91). In both
cases the Abel vector and characteristic come from the paper and use its
homology basis and base point. They must not be mixed directly with the
characteristic returned for Genera's automatically marked curve. This
example computes the periods from the curve, but does not reconstruct the
paper's divisors or their Abelian integration paths.

Executable validation
----------------------

The branch-locus plot and this doctest use the same checked calculation at
20 decimal digits.

.. doctest::

   >>> from examples.bernatska.bernatska_trigonal_demo import documentation_example
   >>> result = documentation_example()
   >>> result.genus, result.intersection_rank, len(result.branch_values)
   (3, 6, 8)
   >>> len(result.values_3a), len(result.values_3b)
   (8, 8)
   >>> result.residuals["symplectic_change"] == 0
   True
   >>> max(result.residuals["p_3a"], result.residuals["p_3b"]) < 2e-3
   True

The helper raises an error if any of the topology, basis, period, symmetry,
positive-definiteness or P-function checks fail. Representative maximum
absolute residuals and the enforced tolerances are:

.. list-table::
   :header-rows: 1

   * - Check
     - Typical residual
     - Tolerance
   * - Finite branch values against the paper
     - :math:`4.2\times10^{-6}`
     - :math:`8\times10^{-6}`
   * - Distance of the cycle change to integers
     - :math:`8.2\times10^{-7}`
     - :math:`2\times10^{-6}`
   * - First-kind full periods
     - :math:`5.5\times10^{-7}`
     - :math:`2\times10^{-6}`
   * - Second-kind full periods
     - :math:`2.7\times10^{-6}`
     - :math:`2\times10^{-5}`
   * - Symmetry of :math:`\tau`
     - :math:`2.5\times10^{-21}`
     - :math:`10^{-12}`
   * - Symmetry of :math:`\varkappa`
     - :math:`8.4\times10^{-21}`
     - :math:`10^{-11}`
   * - Example 3a P-functions
     - :math:`4.8\times10^{-5}`
     - :math:`2\times10^{-3}`
   * - Example 3b P-functions
     - :math:`1.4\times10^{-4}`
     - :math:`2\times10^{-3}`

A representative local run takes about 4.7 seconds before plot rendering.
The period and Abel-vector tables contain only rounded decimal data, so
working at 20 digits does not imply 20-digit agreement with the paper.
P-function derivatives can amplify perturbations in the supplied Abel
vectors. The table comparisons therefore use looser tolerances than the
internal symmetry checks; increasing working precision alone cannot recover
the missing digits of the reference inputs.

Running the standalone example
-------------------------------

From the repository root::

   python -m examples.bernatska.bernatska_trigonal_demo

This prints the branch permutations, cycle transformation and full
validation diagnostics. Genera and mpmath are sufficient; Matplotlib is
imported only when requesting a figure.
