Computing the Bobenko curve data
================================

This companion to :doc:`bobenko` computes first-kind periods, marked-point
Abel shifts and the velocity of the Jacobian flow from the Kowalewski
spectral curve [BRS1989]_. It demonstrates custom differentials and local
charts, including a closed-loop integral for a residue at infinity.

The calculation uses Genera's geometric marking throughout. The theta
solution in the other example uses an involution-adapted marking, so these
vectors need a compatible change of cycles before they can replace its
supplied data. The third-kind normalization scalar and recovery of physical
initial conditions are outside this companion.

The curve and forms
--------------------

For :math:`H=3/2`, :math:`I_1=1/5` and :math:`I_2=27/5`, the model is

.. math::

   F(x,y)=x^2y^4-2(x-3x^2+2x^3)y^2+1-\frac{26}{5}x+\frac{27}{5}x^2=0.

Its compact smooth model has genus three. The chosen holomorphic basis is

.. math::

   \begin{aligned}
   du_1&=\frac{x}{F_y}\,dx,\\
   du_2&=\frac{xy}{F_y}\,dx,\\
   du_3&=\frac{xy^2-1}{F_y}\,dx.
   \end{aligned}

The affine equation has special behaviour at its boundary, so the basis
and local charts are explicit. They describe the smooth model; this is not
a demonstration of automatic singularity resolution.

.. literalinclude:: ../../examples/bobenko_reyman_semenov_tian_shansky/curve_data.py
   :language: python
   :start-at: def make_curve():
   :end-before: def zero_chart(

.. _bobenko-marked-places:

Marked places and Abel shifts
------------------------------

There are two places above :math:`x=0` and two above infinity. At zero,
:math:`x=t^2` and :math:`y=1/(t(1+tw))` resolve the two branches;
``chart_fibre(chart, 0)`` distinguishes their seeds. At infinity the charts
are

.. math::

   \begin{aligned}
   x&=t^{-2},\\
   y_{\mathrm{growing}}&=t^{-1}w,\\
   y_{\mathrm{vanishing}}&=tw.
   \end{aligned}

The growing seed is :math:`w(0)=2`; the vanishing seed is
:math:`w(0)=\sqrt{27/20}`. ``chart_place`` represents each endpoint with a
finite junction and its local tail. ``abel_map_kind_1`` includes that tail
when computing the integral to the marked place.

With the growing infinity place as origin, the normalized shifts are

.. math::

   z(P)=(2\omega)^{-1}\bigl(A(P)-A(P_{\infty,+})\bigr).

The program checks Abel's theorem for the principal divisor of :math:`x`:

.. math::

   \operatorname{div}(x)=2P_{0,0}+2P_{0,1}
                          -2P_{\infty,+}-2P_{\infty,-}.

Its Abel image is reduced with ``lattice_reduce`` and must vanish modulo
full periods. This checks the four marked places together.

.. _bobenko-flow-residues:

Residues and flow velocity
---------------------------

A residue is the coefficient of :math:`dt/t` in a meromorphic differential
near a pole. It can be calculated without a symbolic local series:

.. math::

   \operatorname{Res}_{t=0}\alpha
      =\frac{1}{2\pi i}\oint_{|t|=r}\alpha.

The loop surrounds the pole and does not pass through it. Here the forms
are :math:`\alpha_j=y\,du_j/2`. Equation (7.9) gives the normalized flow
velocity from their residue vector:

.. math::

   \begin{aligned}
   \rho_j&=\operatorname{Res}_{P_{\infty,+}}(y\,du_j/2),\\
   \rho&=(0,0,-1/2)^T,\\
   v&=(2\omega)^{-1}\rho.
   \end{aligned}

``chart_integral`` integrates these forms along a closed polygonal loop
in the local :math:`t` coordinate, following the chosen branch:

.. literalinclude:: ../../examples/bobenko_reyman_semenov_tian_shansky/curve_data.py
   :language: python
   :start-at: def velocity_residues(
   :end-before: def residue_example(

The residue part is independently executable without computing periods:

.. doctest::

   >>> from mpmath import mp
   >>> from examples.bobenko_reyman_semenov_tian_shansky.curve_data import residue_example
   >>> with mp.workdps(20):
   ...     error = mp.norm(residue_example() - mp.matrix([0, 0, -mp.mpf(1)/2]))
   >>> error < 1e-18
   True

At 20 decimal digits, the complete calculation gives a residue residual
of about :math:`4\times10^{-21}`, a change below :math:`10^{-21}` when the
loop radius is halved, and a principal-divisor residual below
:math:`10^{-20}`. These are absolute numerical checks, not rigorous error
bounds. The program enforces a tolerance of :math:`10^{-10}`.

The full program prints the period matrix, marked shifts and velocity::

   python -m examples.bobenko_reyman_semenov_tian_shansky.curve_data

The period and marked-place calculation runs separately from the short
residue doctest; it uses no external period or marked-point reference data.
