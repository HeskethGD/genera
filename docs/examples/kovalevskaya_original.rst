Kovalevskaya's original theta-function solution
===============================================

This example follows S. Kowalevski's *Sur le problème de la rotation d'un
corps solide autour d'un point fixe* [Kowalevski1889]_, sections 2 and 4--7.
It constructs a genus-two spectral curve,
computes its periods and Abel map with Genera, and reconstructs the motion
from Riemann theta quotients. An independent fourth-order Runge--Kutta
integration (RK4) checks the result.

This is the original memoir's construction. The separate
``examples/bobenko_reyman_semenov_tian_shansky/`` example follows the later
Lax-pair treatment. The parameters and initial state here are chosen for
this demonstration, rather than copied from a numerical example in the
memoir.

The equations of motion
-----------------------

Write the angular velocities as :math:`p,q,r` and the direction cosines as
:math:`g_0,g_1,g_2`. In the memoir's normalization the Euler--Poisson system is

.. math::

   \begin{aligned}
   2\dot p&=qr,\\
   2\dot q&=-pr-c_0g_2,\\
   \dot r&=c_0g_1,\\
   \dot g_0&=rg_1-qg_2,\\
   \dot g_1&=pg_2-rg_0,\\
   \dot g_2&=qg_0-pg_1.
   \end{aligned}

The four conserved quantities are

.. math::

   \begin{aligned}
   2(p^2+q^2)+r^2-2c_0g_0&=6l_1,\\
   2(pg_0+qg_1)+rg_2&=2l,\\
   g_0^2+g_1^2+g_2^2&=1,\\
   \left|(p+iq)^2+c_0(g_0+ig_1)\right|^2&=k^2.
   \end{aligned}

We choose :math:`l_1=2`, :math:`k=3/2`, :math:`c_0=1` and :math:`l=1/2`.
The runnable example constructs an initial state satisfying these integrals
and selects one whose separated coordinates are away from branch points.

The genus-two curve
-------------------

The separation of variables leads to

.. math::

   \begin{aligned}
   y^2&=-4(s-e_1)(s-e_2)(s-e_3)(s-k_1)(s-k_2),\\
   k_1&=\frac{l_1+k}{2},\\
   k_2&=\frac{l_1-k}{2}.
   \end{aligned}

where :math:`e_1>e_2>e_3` solve :math:`4s^3-G_2s-G_3=0`, with

.. math::

   \begin{aligned}
   G_2&=k^2-c_0^2+3l_1^2,\\
   G_3&=l_1(k^2-c_0^2-l_1^2)+l^2c_0^2.
   \end{aligned}

Uppercase :math:`G_2,G_3` distinguish the cubic coefficients from the
direction cosines. The five real branch points satisfy
:math:`k_1>e_1>e_2>k_2>e_3`.

The script expands the polynomial and builds a sparse equation mapping for
:ref:`Curve <curve-class>`. The central API calls are::

   curve = Curve({(0, 2): 1, **{(i, 0): -c for i, c in enumerate(curve_coefficients())}})
   periods = curve.periods_kind_1()
   u0 = curve.abel_map_kind_1(divisor_places(INITIAL)).value

Here ``curve_coefficients`` and ``divisor_places`` are helpers in the
complete example. All periods and Abel coordinates come from the same
curve marking. In particular, ``periods.omega`` contains half-periods,
and the theta argument is :math:`v=(2\omega)^{-1}u`.

From the Abel flow to the physical motion
-----------------------------------------

In Genera's ordered basis :math:`(ds/y,s\,ds/y)`, the selected sheets give
the straight-line flow

.. math::

   u(t)=u(0)+\begin{pmatrix}0\\-1\end{pmatrix}t.

The example checks this against Abel maps of RK4 divisors on
:math:`0\leq t\leq0.6`. Fresh Abel maps can differ by period vectors;
``curve.lattice_reduce`` continues their increments across those jumps.

The memoir reconstructs the motion through fifteen quantities
:math:`P_i,P_{ij}`. Each is represented by a theta quotient

.. math::

   P_\alpha(t)=C_\alpha\,
   \frac{\theta[\chi_\alpha](v(t)\mid\tau)}
        {\theta[K](v(t)\mid\tau)}.

Genera's branch-point Abel images determine the numerator characteristics;
``curve.riemann_constant().characteristic`` supplies :math:`K`. The multiplier
:math:`C_\alpha` is normalized using the algebraic value at the initial
divisor, then held fixed and checked at twelve further divisors. These
checks use the actual Abel maps, independently of the prescribed linear
flow. The memoir's characteristic labels cannot simply be transferred to
Genera's automatic homology marking.

The physical variables follow from the memoir's section-5 formulas. For
example, :math:`q=E/\Delta`, where
:math:`\Delta=LP_1+MP_2+NP_3` and
:math:`E=(e_2-e_3)(e_3-e_1)(e_1-e_2)`. The remaining formulas and their
radical conventions are implemented in ``assemble_state`` and documented
in ``examples/kovalevskaya_original/walkthrough.md``. That walkthrough also
records corrections to the available transcription, justified by algebraic
identities and the equations of motion; this example does not establish
whether those errors originated in the printing or the transcription.

The solution and its numerical check
------------------------------------

The plot evaluates the theta solution on :math:`0\leq t\leq2.4`, beyond a
separated-coordinate turning point near :math:`t=0.73`. Dashed lines show
the dense RK4 trajectory; open markers show the more sparsely evaluated theta
solution. Numerical errors are described below rather than plotted.

.. plot::
   :include-source: True
   :context: reset

   from examples.kovalevskaya_original.kovalevskaya_original_demo import (
       documentation_example, make_figure,
   )

   result = documentation_example()
   make_figure(result)

The same computation is an executable doctest. It uses 30 decimal digits,
13 divisors for the algebraic/theta checks, 241 RK4 samples and 21 theta
samples for the plotted motion. RK4 uses 128 substeps between its samples,
giving an integration step of :math:`h=0.000078125`. A second integration
uses 64 substeps (:math:`2h`) to check convergence. It omits the standalone
:math:`[0,20]` RK4 drift experiment.

.. doctest::

   >>> from examples.kovalevskaya_original.kovalevskaya_original_demo import documentation_example
   >>> result = documentation_example()
   >>> len(result.rk4_times), len(result.theta_times), len(result.errors)
   (241, 21, 6)
   >>> max(result.errors) < 1e-14
   True
   >>> 8 < result.rk4_error_ratio < 32
   True

``documentation_example`` also raises an error if branch-point agreement,
Abel linearity, held-out theta quotients, the initial state, reality,
conserved quantities or directly differentiated Euler--Poisson equations
fail their checks. Analytic checks use an absolute tolerance of
:math:`10^{-24}` (the theta-quotient check is scaled), Abel linearity uses
:math:`10^{-8}`, and the RK4 comparison uses :math:`10^{-14}`.
The two RK4 trajectories must differ by less than :math:`10^{-13}`;
halving the step must reduce the maximum comparison error by a factor
between eight and 32, surrounding the fourth-order prediction of 16.

At 30 digits, a representative local run takes about nine seconds before
plot rendering. Curve construction, Abel maps, theta calculations and
analytic checks take about 3.9 seconds. The fine RK4 trajectory takes
about 3.6 seconds, and the coarse convergence check about 1.8 seconds.
Its maximum theta/RK4 difference is approximately
:math:`2.5\times10^{-16}`, while analytic invariant and ODE residuals are
of order :math:`10^{-29}`. The RK4 comparison therefore measures integration
error rather than 30-digit accuracy. Increasing working precision alone
will not improve that comparison; a smaller RK4 step is needed.
Reported comparison errors are measured at the 21 theta samples, rather
than bounding every time in the interval.

Running the complete example
----------------------------

From the repository root::

   python -m examples.kovalevskaya_original.kovalevskaya_original_demo --stage demo

This prints detailed residuals without requiring Matplotlib. Omitting
``--stage demo`` additionally runs the longer RK4 invariant-drift check.
