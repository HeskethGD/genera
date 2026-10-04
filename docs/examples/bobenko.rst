The Bobenko–Reyman–Semenov-Tian-Shansky top
==================================================

This example follows A. I. Bobenko, A. G. Reyman and
M. A. Semenov-Tian-Shansky, *The Kowalewski top 99 years later: a Lax pair,
generalizations and explicit solutions* [BRS1989]_. It evaluates the genus-three
Riemann theta solution in Theorem 7.7, equation (7.42), and compares all
six physical variables with independent Euler–Poisson integration.

This is a different representation of the classical top from the
:doc:`original Kovalevskaya example <kovalevskaya_original>`. Here the
period matrix, flow velocity and marked-point shifts are supplied numerical
data in the runnable script. This page checks the theta reconstruction;
it does not reconstruct those data from a spectral curve.
The :doc:`curve-data companion <bobenko_curve_data>` computes periods,
marked-point Abel shifts and the flow velocity in Generapy's geometric marking.

The physical equations
-----------------------

Write the angular momentum as :math:`\ell=(\ell_1,\ell_2,\ell_3)` and
the gravity direction as :math:`g=(g_1,g_2,g_3)`. The normalization used
by this example gives

.. math::

   \begin{aligned}
   \dot\ell_1&=\ell_2\ell_3,\\
   \dot\ell_2&=-\ell_1\ell_3-g_3,\\
   \dot\ell_3&=g_2,\\
   \dot g_1&=2g_2\ell_3-g_3\ell_2,\\
   \dot g_2&=g_3\ell_1-2g_1\ell_3,\\
   \dot g_3&=g_1\ell_2-g_2\ell_1.
   \end{aligned}

The four conserved quantities checked along the trajectory are

.. math::

   \begin{aligned}
   N&=g_1^2+g_2^2+g_3^2=1,\\
   H&=\frac{\ell_1^2+\ell_2^2+2\ell_3^2}{2}-g_1,\\
   I_1&=(\ell_1g_1+\ell_2g_2+\ell_3g_3)^2,\\
   I_2&=(\ell_1^2-\ell_2^2+2g_1)^2+4(\ell_1\ell_2+g_2)^2.
   \end{aligned}

.. _bobenko-theta-example:

Theta values and derivatives
-----------------------------

Let :math:`\tau` be the supplied symmetric genus-three period matrix,
:math:`v` the supplied velocity, and

.. math::

   \begin{aligned}
   z(t)&=(p_1/2,p_2,p_1/2)^T+vt,\\
   \epsilon&=\left[\begin{smallmatrix}0&0&0\\0&1/2&0\end{smallmatrix}\right].
   \end{aligned}

Write :math:`\theta(z)=\theta(z\mid\tau)` and
:math:`\theta_\epsilon(z)=\theta[\epsilon](z\mid\tau)`.
``rtheta`` evaluates these normalized Riemann theta functions directly.
``rtheta_jet`` supplies their value and first derivatives in one call:

.. literalinclude:: ../../examples/bobenko_reyman_semenov_tian_shansky/kowalewski_genus_three.py
   :language: python
   :start-at: def theta_data(
   :end-before: def solution(

The tuple keys in the jet are multi-indices. Contracting the gradient
with :math:`v` gives the time derivative :math:`D_v\theta`.
No numerical differentiation of a trajectory is needed.

Reconstructing the top
-----------------------

Denote the supplied Abelian shift by :math:`r`, the two marked shifts by
:math:`s_-` and :math:`s_+`, and the supplied scalar by :math:`c_3`.
All theta functions below are evaluated at the moving phase:

.. math::

   \begin{aligned}
   A&=\theta(z+s_-),\\
   B&=\theta(z+s_+),\\
   C&=-\theta_\epsilon(z+s_-),\\
   D&=-\theta_\epsilon(z+s_+),\\
   \Delta&=AD+BC.
   \end{aligned}

The minus signs in :math:`C,D` account for the chosen marked paths:
shifting those paths by an odd second B-cycle supplies a multiplier of
:math:`-1` for the characteristic :math:`\epsilon`.
With :math:`\ell_\pm=\ell_1\pm i\ell_2` and
:math:`g_\pm=g_1\pm ig_2`, the reconstruction is

.. math::

   \begin{aligned}
   \ell_-&=2ic_3\frac{\theta(z-r)}{\theta(z)},\\
   \ell_+&=2ic_3\frac{\theta_\epsilon(z-r)}{\theta_\epsilon(z)},\\
   \ell_3&=-i\left(\frac{D_v\theta_\epsilon(z)}{\theta_\epsilon(z)}
                  -\frac{D_v\theta(z)}{\theta(z)}\right),\\
   g_-&=2\frac{\theta_\epsilon(z)}{\theta(z)}\frac{AB}{\Delta},\\
   g_+&=2\frac{\theta(z)}{\theta_\epsilon(z)}\frac{CD}{\Delta},\\
   g_3&=-\frac{AD-BC}{\Delta}.
   \end{aligned}

The runnable implementation retains complex values through the entire
comparison. It checks their small imaginary parts before the plot shows
the real components.

The trajectory
---------------

The fixed phase selects a numerically real trajectory. RK4 starts from
the theta solution at :math:`t=0` and then integrates the physical equations
without further theta evaluations. Dashed curves show 257 saved RK4
states on :math:`[0,2]`; open markers show 17 theta evaluations.

.. plot::
   :include-source: True
   :context: reset

   from examples.bobenko_reyman_semenov_tian_shansky.kowalewski_genus_three import (
       documentation_example, make_figure,
   )

   result = documentation_example()
   make_figure(result)

Executable validation and precision
------------------------------------

The same checked calculation runs as a doctest at 30 decimal digits:

.. doctest::

   >>> from examples.bobenko_reyman_semenov_tian_shansky.kowalewski_genus_three import documentation_example
   >>> result = documentation_example()
   >>> len(result.rk4_times), len(result.analytic_times), len(result.errors)
   (257, 17, 6)
   >>> max(result.errors) < 1e-10
   True
   >>> 8 < result.rk4_refinement_ratio < 32
   True

Working at 30 digits limits arithmetic roundoff; it cannot restore digits
missing from the supplied period and marked-point data. These data have
roughly 14–16 decimal places. Independently integrated period entries
differ from their transposes by about :math:`2\times10^{-14}`; the script
projects both period matrices onto their symmetric parts before evaluation.
The marked shifts introduce further error into the reconstructed motion.

For this reason the helper permits absolute errors of :math:`10^{-10}`
in the six-component comparison, imaginary parts, analytic invariant drift
and initial gravity norm. It also checks equation (7.60):
:math:`\Delta/[\theta(z)\theta_\epsilon(z)]` must be constant across two
generic complex phases and two times, to a scaled tolerance of
:math:`10^{-10}`. The two cases of the genus-three/Prym/elliptic theta
decomposition in equation (7.64) are checked at both phases with scaled
tolerance :math:`10^{-12}`.

The RK4 check uses 8,192 steps, :math:`h=2/8192`, with additional runs
at 4,096 and 2,048 steps. The maximum difference between consecutive grids
must decrease by a factor between eight and 32. This checks fourth-order
convergence independently of the rounded theta data. Fine/middle trajectory
differences and numerical invariant drift must both be below
:math:`10^{-12}`. All comparisons use the saved samples; they do not bound
every point of the continuous motion.

A representative local run takes about 5.3 seconds before plot rendering.
The largest theta/RK4 component difference is approximately
:math:`7.2\times10^{-12}`, and analytic invariant drift is
:math:`1.2\times10^{-11}`. By comparison, the fine/middle RK4 difference
is :math:`9.6\times10^{-15}`, numerical invariant drift is
:math:`1.2\times10^{-15}`, and the refinement ratio is 15.99.
Further RK4 refinement therefore has little effect on agreement with this
particular theta fixture.

Running the standalone example
-------------------------------

From the repository root::

   python -m examples.bobenko_reyman_semenov_tian_shansky.kowalewski_genus_three

The standalone program offers ``--dps``, ``--stop``, ``--steps``,
``--samples`` and ``--tol``. Its default RK4 grid is coarser than the
documentation calculation. Matplotlib is imported only when making a figure.
