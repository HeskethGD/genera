Neumann–Moser motion on a genus-three curve
===========================================

This example follows P. G. Baron's
*The Neumann–Moser dynamical system and the Korteweg–de Vries hierarchy*
[Baron2024]_, using the hyperelliptic Kleinian
solution developed in V. M. Buchstaber's
*The Mumford dynamical system and hyperelliptic Kleinian functions*
[Buchstaber2024]_. It reconstructs a ten-coordinate
trajectory from genus-three P-functions, checks the conserved spectral
polynomial and compares the motion against independent RK4 integration.

The curve and Abelian phase are chosen for this demonstration. The
standalone examples also include a genus-two version; this page uses
genus three to exercise a degree-three divisor and additional Kleinian
identities.
The genus-two calculation below also starts from independently specified
physical coefficients and computes their curve and Abelian phase.

The equations of motion
-----------------------

The state consists of the coefficients
:math:`(u_1,u_2,u_3,v_1,v_2,v_3,w_1,w_2,w_3,w_4)`. With
:math:`\Gamma=w_1-u_1`, its equations are

.. math::

   \begin{aligned}
   \dot u_1&=-2v_1,\\
   \dot u_2&=-2v_2,\\
   \dot u_3&=-2v_3,\\
   \dot v_1&=-\Gamma u_1-u_2+w_2,\\
   \dot v_2&=-\Gamma u_2-u_3+w_3,\\
   \dot v_3&=-\Gamma u_3+w_4,\\
   \dot w_1&=2v_1,\\
   \dot w_2&=2v_2+2\Gamma v_1,\\
   \dot w_3&=2v_3+2\Gamma v_2,\\
   \dot w_4&=2\Gamma v_3.
   \end{aligned}

These coefficients define three generating polynomials:

.. math::

   \begin{aligned}
   U(\xi)&=\xi^3+u_1\xi^2+u_2\xi+u_3,\\
   V(\xi)&=v_1\xi^2+v_2\xi+v_3,\\
   W(\xi)&=\xi^4+w_1\xi^3+w_2\xi^2+w_3\xi+w_4.
   \end{aligned}

The conserved relation

.. math::

   U(\xi)W(\xi)+V(\xi)^2=\frac{F(\xi)}{4}

ties the motion to a fixed spectral curve. Checking every coefficient of
this polynomial supplies several independent invariant checks, rather
than checking only one scalar energy.

The spectral curve and period data
-----------------------------------

The demonstration curve is :math:`y^2=F(s)`, with

.. math::

   \begin{aligned}
   F(s)&=4\prod_{a\in B}(s-a),\\
   F(s)&=4s^7+\sum_{j=2}^{7}\lambda_{2j}s^{7-j}.
   \end{aligned}

The seven distinct real branch points are
:math:`B=\{-5,-3,-1.2,-0.3,0.8,2.7,6\}`. Their sum is zero, giving
the canonical form without an :math:`s^6` term. For this curve,
:math:`\lambda_4=-158.92` and :math:`\lambda_6=-109.92`.

The runnable script expands the polynomial into ascending coefficients
and constructs compatible first-kind, second-kind and Riemann-constant
data from one curve object::

   from generapy import Curve, kleinian_p

   curve = Curve({(0, 2): 1, **{(i, 0): -c for i, c in enumerate(coefficients)}})
   first = curve.periods_kind_1()
   second = curve.periods_kind_2()
   omega, tau = first.omega, first.tau
   kappa = second.kappa
   characteristic = curve.riemann_constant().characteristic

The half-period matrix ``omega``, normalized matrix ``tau``, second-kind
combination ``kappa`` and literal theta characteristic are reused for
every evaluation along the trajectory.

From Kleinian functions to the state
------------------------------------

The paper labels Abelian coordinates by weight as :math:`z_1,z_3,z_5`.
Generapy orders its differential basis by ascending powers of :math:`s`,
so these correspond to API coordinates 2, 1 and 0 respectively.
The P-functions used here map as follows:

.. list-table::
   :header-rows: 1

   * - Paper
     - ``kleinian_p`` indices
   * - :math:`\wp_2`
     - ``(2, 2)``
   * - :math:`\wp_4`
     - ``(1, 2)``
   * - :math:`\wp_6`
     - ``(0, 2)``
   * - :math:`\wp_{3,3}`
     - ``(1, 1)``
   * - :math:`\wp_{3,5}`
     - ``(1, 0)``

A prime differentiates along :math:`z_1` and appends an index 2.
Two primes append two index 2s. Thus the following batched call evaluates
nine P-functions of orders two, three and four from one theta jet:

.. literalinclude:: ../../examples/neumann_moser/neumann_moser_genus_three_demo.py
   :language: python
   :start-at:     wp2, wp4, wp6, wp2p, wp4p, wp6p, wp2pp, wp4pp, wp6pp = kleinian_p(
   :end-before:     coordinates = (
   :dedent: 4

The physical coordinates are then reconstructed by

.. math::

   \begin{aligned}
   u_i&=-\wp_{2i},\quad i=1,2,3,\\
   v_i&=\tfrac12\wp'_{2i},\quad i=1,2,3,\\
   w_1&=\wp_2,\\
   w_2&=\tfrac12\wp''_2-\wp_4-2\wp_2^2,\\
   w_3&=\tfrac12\wp''_4-\wp_6-2\wp_2\wp_4,\\
   w_4&=\tfrac12\wp''_6-2\wp_2\wp_6.
   \end{aligned}

The implementation checks that all ten coordinates are real before
discarding their small imaginary roundoff components.

The selected Abelian phase
--------------------------

To distinguish the Abelian vector from the dynamical coefficients
:math:`u_i`, write Generapy's vector as :math:`z_G`. This example uses

.. math::

   z_G(t)=\begin{pmatrix}0.13\\0.27\\0.41+t\end{pmatrix}
          +\omega\tau\begin{pmatrix}1\\1\\1\end{pmatrix}.

Only the last coordinate advances. The imaginary half-period translate
selects a real component on which the sampled motion stays bounded.
The other two coordinates fix the phase; the first is the ``0.13`` phase
selected from the standalone example.

Recovering the divisor and phase
---------------------------------

The roots :math:`s_j` of :math:`U(s)` recover a degree-three divisor.
At each root, the spectral invariant gives :math:`V(s_j)^2=F(s_j)/4`,
so the sheet coordinates are

.. math::

   y_j=2V(s_j).

This factor of two is required by the chosen curve normalization.
The inverse loop uses only the physical state and its polynomial
coefficients to construct the divisor::

   divisor = divisor_from_state(state)
   image = data["curve"].abel_map_kind_1(divisor, reduce=True).value

At :math:`t=0`, the example checks that all three points lie on the
curve, that the image recovers the chosen phase modulo full periods,
and that evaluating at the recovered representative reproduces all
ten state coordinates.

Conventions and identity checks
-------------------------------

The implementation uses the generating-polynomial equations to resolve
inconsistencies in the expanded formulas of the cited paper. The middle
:math:`\dot w_i` range includes :math:`i=3`, and the final equation uses
:math:`v_3`. In the solution for :math:`W`, the factor is the polynomial
:math:`p_I` of degree three; the printed :math:`p_{II}` would give the
wrong degree. The original walkthroughs in ``examples/neumann_moser/``
record these details and the comparison with Buchstaber's solution.

In Generapy's sigma normalization, three additional identities are checked
at the phase origin:

.. math::

   \begin{aligned}
   \wp''_2&=6\wp_2^2+4\wp_4+\tfrac12\lambda_4,\\
   \wp''_4&=6(\wp_2\wp_4+\wp_6)-2\wp_{3,3},\\
   \wp''_6&=6\wp_2\wp_6-2\wp_{3,5}.
   \end{aligned}

The :math:`\lambda_4/2` coefficient is tied to this sigma and curve
normalization. The last two identities exercise quantities that are
absent from the genus-two example.

.. _neumann-moser-initial-conditions:

Starting from physical initial conditions (genus two)
-----------------------------------------------------

The genus-two script also accepts a seven-coordinate initial state,
without choosing an Abelian phase first. For example,

.. math::

   \begin{aligned}
   (u_1,u_2,v_1,v_2,w_1,w_2,w_3)
      &=(-1,-2,1/10,1/5,1,-55/4,-7),\\
   U(s)&=(s+1)(s-2),\\
   V(s)&=s/10+1/5,\\
   W(s)&=(s+4)(s+1/2)(s-7/2),\\
   F(s)&=4\bigl(U(s)W(s)+V(s)^2\bigr).
   \end{aligned}

The invariant constructs the curve. The roots of :math:`U` and the sheets
:math:`y_j=2V(s_j)` give its initial divisor; ``abel_map_kind_1`` then
determines the offset. This variant uses the canonical form
:math:`w_1=-u_1`, giving no :math:`s^4` term in :math:`F`.

.. literalinclude:: ../../examples/neumann_moser/neumann_moser_kleinian_demo.py
   :language: python
   :start-at: def data_from_initial_state(
   :end-before: def period_lattice_residual(

The recovered Kleinian solution is compared with RK4 initialized directly
from the supplied coefficients, rather than from an analytic evaluation:

.. doctest::

   >>> from mpmath import mp
   >>> from examples.neumann_moser.neumann_moser_kleinian_demo import initial_value_example
   >>> with mp.workdps(25):
   ...     rows = initial_value_example()
   ...     coarse = initial_value_example(steps=64)
   >>> max(rows[0][3]) < 1e-20
   True
   >>> fine_error = max(max(row[3]) for row in rows)
   >>> fine_error < 1e-11
   True
   >>> 8 < max(max(row[3]) for row in coarse) / fine_error < 32
   True

On :math:`[0,0.1]` at 25 decimal digits, the initial reconstruction residual
is about :math:`10^{-25}` and the maximum sampled state error with 128 RK4
steps is about :math:`9\times10^{-13}`. Halving the step size reduces the
error by about a factor of sixteen. The initial coefficients are rational
inputs; no phase or trajectory value is fitted.

The same workflow is available from the command line::

   python -m examples.neumann_moser.neumann_moser_kleinian_demo \
       --initial-state -1 -2 0.1 0.2 1 -13.75 -7 \
       --start=0 --stop=0.1 --steps 128 --samples 5

The trajectory
---------------

The three panels show the coefficients of :math:`U`, :math:`V` and
:math:`W`. Dashed curves are the dense RK4 trajectory; markers are the
sparse Kleinian evaluations. RK4 starts from the analytic state at
:math:`t=-0.5` and receives no later analytic values during integration.

.. plot::
   :include-source: True
   :context: reset

   from examples.neumann_moser.neumann_moser_genus_three_demo import (
       documentation_example, make_figure,
   )

   result = documentation_example()
   make_figure(result)

Executable validation
----------------------

The plot and doctest share one checked calculation at 30 decimal digits.
The interval :math:`[-0.5,0.5]` keeps this higher-genus example practical
for documentation builds. It saves 129 RK4 plotting samples and evaluates
17 analytic states. The fine integration uses 32,768 steps, and a second
integration uses 16,384 steps to check fourth-order convergence.

.. doctest::

   >>> from examples.neumann_moser.neumann_moser_genus_three_demo import documentation_example
   >>> result = documentation_example()
   >>> len(result.rk4_times), len(result.analytic_times), len(result.errors)
   (129, 17, 10)
   >>> max(result.errors) < 1e-14
   True
   >>> result.analytic_spectral_error < 1e-24
   True
   >>> 8 < result.rk4_error_ratio < 32
   True

The helper enforces all three Kleinian identities, the divisor/Abel-map
checks and the analytic spectral polynomial at an absolute tolerance of
:math:`10^{-24}`. The numerical state comparison and spectral polynomial
must agree within :math:`10^{-14}`. Fine and coarse RK4 trajectories
must differ by less than :math:`10^{-13}`, with a comparison-error
reduction between eight and 32 when halving the step.

Comparison errors are measured at the 17 analytic samples. Spectral
coefficients are checked at every saved RK4 sample and every analytic
sample. These checks do not bound the continuous trajectory between
samples. Roundoff and integration error are controlled separately by
working precision and RK4 step size.

A representative local run takes about 12.7 seconds before plot rendering.
The maximum state error is approximately :math:`5.7\times10^{-16}`, and
the numerical spectral-polynomial error is :math:`1.2\times10^{-15}`.
Halving the RK4 step reduces the state comparison error by a factor of
16.0013. The analytic spectral error is about :math:`1.5\times10^{-28}`;
the divisor and identity checks remain below :math:`10^{-26}`.

Running the standalone example
-------------------------------

From the repository root::

   python -m examples.neumann_moser.neumann_moser_genus_three_demo

The standalone program defaults to three phases over :math:`[-4,4]`.
That longer run is kept outside the documentation calculation; precision,
interval, sampling and phase options remain available. Matplotlib is
imported only when requesting a figure.
