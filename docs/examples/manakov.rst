Manakov motion from genus-two Kleinian functions
================================================

This example follows Christiansen, Eilbeck, Enolskii and Kostov,
*Quasi-periodic and periodic solutions for systems of coupled nonlinear
Schrödinger equations* [CEEK2000]_. It uses the stationary
real reduction in their equation (2.1): two coupled coordinates are
reconstructed from genus-two Kleinian P-functions and checked against
independent RK4 integration.

The curve and phase are chosen for this demonstration. This calculation
covers the real stationary motion; the complete complex wave fields also
require the Baker--Akhiezer phase integrals.

The stationary equations
-------------------------

For :math:`i=1,2`, the equations are

.. math::

   q_i''=-(q_1^2+q_2^2)q_i+a_iq_i+\frac{C_i^2}{q_i^3}.

The independent variable :math:`x` is the stationary coordinate. Writing
:math:`p_i=q_i'` gives the first-order state
:math:`(q_1,p_1,q_2,p_2)` integrated by RK4. The conserved Hamiltonian is

.. math::

   H=\frac{p_1^2+p_2^2}{2}
     +\frac{(q_1^2+q_2^2)^2}{4}
     -\frac{a_1q_1^2+a_2q_2^2}{2}
     +\frac12\left(\frac{C_1^2}{q_1^2}+\frac{C_2^2}{q_2^2}\right).

.. _manakov-period-data:

The spectral curve and compatible periods
-----------------------------------------

We use the nonsingular genus-two curve

.. math::

   y^2=F(s)=4s(s-\tfrac18)(s-\tfrac14)(s-\tfrac12)(s-1),

with marked parameters :math:`a_1=3/4`, :math:`a_2=3/16`, and

.. math::

   \begin{aligned}
   d&=a_1-a_2,\\
   C_i^2&=-\frac{F(a_i)}{d^2}.
   \end{aligned}

The runnable script expands :math:`F` into ascending polynomial
coefficients. Its period construction uses one curve object throughout::

   from genera import Curve, kleinian_p

   curve = Curve({(0, 2): 1, **{(i, 0): -c for i, c in enumerate(coefficients)}})
   first = curve.periods_kind_1()
   second = curve.periods_kind_2()
   omega, tau = first.omega, first.tau
   kappa = second.kappa
   characteristic = curve.riemann_constant().characteristic

These are mutually compatible half-periods, normalized periods,
second-kind data and a theta characteristic. Reusing them avoids computing
periods again at each point of the trajectory.

.. _manakov-kleinian-example:

The Kleinian solution
----------------------

Genera orders the Abelian differential basis as :math:`(ds/y,s\,ds/y)`.
The paper's one-based indices become zero-based API indices:

.. list-table::
   :header-rows: 1

   * - Paper
     - ``kleinian_p`` indices
   * - :math:`\wp_{22}`
     - ``(1, 1)``
   * - :math:`\wp_{12}`
     - ``(0, 1)``
   * - :math:`\wp_{222}`
     - ``(1, 1, 1)``
   * - :math:`\wp_{122}`
     - ``(0, 1, 1)``

The Abelian point moves along :math:`u(x)=u_0+(0,x)^T`. With all
P-functions evaluated there, the real amplitudes and momenta are

.. math::

   \begin{aligned}
   q_1^2&=\frac{2(a_1^2-a_1\wp_{22}-\wp_{12})}{d},\\
   q_2^2&=-\frac{2(a_2^2-a_2\wp_{22}-\wp_{12})}{d},\\
   p_1&=-\frac{a_1\wp_{222}+\wp_{122}}{dq_1},\\
   p_2&=\frac{a_2\wp_{222}+\wp_{122}}{dq_2}.
   \end{aligned}

One batched call evaluates the four required P-functions:

.. literalinclude:: ../../examples/manakov/manakov_kleinian_demo.py
   :language: python
   :start-at:     wp22, wp12, wp222, wp122 = kleinian_p(
   :end-before:     a1 = data["a1"]
   :dedent: 4

The initial phase used here is

.. math::

   u_0=\begin{pmatrix}0.13\\0.27\end{pmatrix}
       +\omega\begin{pmatrix}1\\1\end{pmatrix}
       +\omega\tau\begin{pmatrix}1\\1\end{pmatrix}.

The half-period translate selects a real component on which both
:math:`q_i^2` remain positive over this interval. The positive square roots
then give continuous amplitudes. The implementation checks reality and
positivity before forming the physical state.

The paper prints a doubled Abelian velocity below equation (2.12), while
equations (3.18)--(3.19) identify :math:`x=u_2`. The example uses unit
velocity in Genera's second coordinate, as checked by the stationary ODE
comparison. Details of this convention and the fourth-order Kleinian
identities are recorded in ``examples/manakov/walkthrough.md``.

Recovering the Abelian phase
----------------------------

The physical state also determines a degree-two divisor on the curve.
Solving the two amplitude and two momentum formulas recovers
:math:`\wp_{22},\wp_{12},\wp_{222},\wp_{122}` without another theta
evaluation. The roots :math:`s_j` of

.. math::

   U(s)=s^2-\wp_{22}s-\wp_{12}

give the divisor's first coordinates; their sheets are fixed by
:math:`y_j=\wp_{222}s_j+\wp_{122}`. The inverse loop uses::

   divisor = divisor_from_state(initial_state, data)
   image = data["curve"].abel_map_kind_1(divisor, reduce=True).value

The example checks that each point lies on the spectral curve, that
``image`` agrees with the original :math:`u_0` modulo full periods, and
that substituting the recovered image into the Kleinian formulas gives
the same physical state. This checks compatibility of the divisor sheets,
Abel-map paths and period basis.

The trajectory
---------------

Dashed curves show the independently integrated RK4 trajectory. Markers
show the more sparsely evaluated Kleinian solution on :math:`[-4,4]`.
RK4 starts from the analytic state at :math:`x=-4` and receives no further
analytic values during integration.

.. plot::
   :include-source: True
   :context: reset

   from examples.manakov.manakov_kleinian_demo import (
       documentation_example, make_figure,
   )

   result = documentation_example()
   make_figure(result)

Executable validation
----------------------

The same calculation runs as a doctest at 30 decimal digits. It saves
257 RK4 plotting samples and evaluates only 33 analytic states. The fine
integration uses 65,536 steps, :math:`h=8/65536`, and the coarse check uses
32,768 steps. Halving the step must reduce the maximum comparison error
by a factor between eight and 32, surrounding RK4's fourth-order
prediction of 16.

.. doctest::

   >>> from examples.manakov.manakov_kleinian_demo import documentation_example
   >>> result = documentation_example()
   >>> len(result.rk4_times), len(result.analytic_times), len(result.errors)
   (257, 33, 4)
   >>> max(result.errors) < 1e-14
   True
   >>> 8 < result.rk4_error_ratio < 32
   True

The helper also enforces fourth-order Kleinian identities, divisor/Abel-map
closure and analytic Hamiltonian conservation to an absolute tolerance of
:math:`10^{-24}`. Numerical Hamiltonian drift must stay below
:math:`10^{-14}`, and the maximum change between fine and coarse RK4
trajectories must stay below :math:`10^{-13}`. Comparison errors are
measured at the 33 analytic samples; they do not bound every point of the
continuous interval.

Working precision controls roundoff, while the RK4 step controls integration
error. The numerical comparison and the analytic identity checks therefore
use separate tolerances.

A representative local run takes about 10.4 seconds before plot rendering.
Its maximum component error is approximately :math:`6.1\times10^{-16}`,
the numerical Hamiltonian drift is :math:`2.2\times10^{-17}`, and halving
the RK4 step reduces the comparison error by a factor of 16.0008. Kleinian
identity and Abel-map residuals are around :math:`10^{-30}` or smaller.

Running the standalone example
-------------------------------

From the repository root::

   python -m examples.manakov.manakov_kleinian_demo

The standalone program offers precision, interval, sampling and phase
options; it defaults to three phases and a coarser RK4 grid. The documentation
calculation selects phase ``0.13`` and uses the finer grid above. Matplotlib
is imported only when requesting a figure.
