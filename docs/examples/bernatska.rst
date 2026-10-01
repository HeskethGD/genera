Bernatska's trigonal curve
==========================

This example follows Examples 3 and 3a--3b of J. Bernatska,
`Computation of P-Functions on Plane Algebraic Curves
<https://arxiv.org/abs/2407.05632>`_. It uses a non-hyperelliptic genus-three
curve, so it is a useful first example of Genera's general plane-curve engine.

The complete script is
``examples/bernatska/bernatska_trigonal_demo.py`` and can be run from the
Genera repository root with::

   .venv/bin/python -m examples.bernatska.bernatska_trigonal_demo

The curve
---------

Bernatska's equation (73) is the trigonal :math:`(3,4)` curve

.. math::

   f(x,y)=-y^3+x^4+3x^3+7x^2+16x+9
          +y(4x^2+5x+11)=0.

Genera accepts the sparse coefficient mapping directly. The keys are
``(x_power, y_power)`` pairs::

   def trigonal_curve():
       return {
           (0, 3): -1,
           (4, 0): 1,
           (3, 0): 3,
           (2, 0): 7,
           (1, 0): 16,
           (0, 0): 9,
           (2, 1): 4,
           (1, 1): 5,
           (0, 1): 11,
       }

This is a general plane curve rather than one of the specialized
hyperelliptic models. Genera computes its branch locus, monodromy, genus, and
homology from the projection to the ``x``-plane. The expected genus is three.

Holomorphic and second-kind differentials
------------------------------------------

For a plane curve, a differential written as ``h(x, y) dx / f_y`` is passed to
Genera as a callable returning the coefficient of ``dx``. Here

.. math::

   f_y=-3y^2+4x^2+5x+11,

and the three holomorphic forms are

.. math::

   du=\begin{pmatrix}y\\x\\1\end{pmatrix}\frac{dx}{f_y}.

The example also supplies three second-kind forms. The last numerator is

.. math::

   R_5=5x^2y+9xy+\frac{32}{3}x^2+7y+\frac{40}{3}x.

The corresponding code is kept in the runnable example and has this shape::

   def differentials():
       def denominator(x, y):
           return -3*y**2 + 4*x**2 + 5*x + 11

       first_kind = (
           lambda x, y: y / denominator(x, y),
           lambda x, y: x / denominator(x, y),
           lambda x, y: 1 / denominator(x, y),
       )

       def r5(x, y):
           return (5*x**2*y + 9*x*y + mp.mpf(32)/3*x**2
                   + 7*y + mp.mpf(40)/3*x)

       second_kind = (
           lambda x, y: x**2 / denominator(x, y),
           lambda x, y: 2*x*y / denominator(x, y),
           lambda x, y: r5(x, y) / denominator(x, y),
       )
       return first_kind, second_kind

Constructing compatible period data
------------------------------------

Because the caller supplies the differential bases, the example constructs a
``Curve`` with the mpmath context and passes the forms explicitly::

   from genera import Curve, kleinian_p
   from mpmath import mp

   mp.dps = 20
   curve = Curve(mp, trigonal_curve())
   first_kind, second_kind = differentials()
   first = curve.periods_kind_1(first_kind)
   second = curve.periods_kind_2(
       differentials=first_kind,
       second_differentials=second_kind,
   )

The returned records contain the first-kind half-periods, second-kind
half-periods, normalized period matrix ``tau``, and the symmetric matrix
``kappa`` used by the Kleinian functions. In particular,

.. math::

   \tau=\omega^{-1}\omega',\qquad
   \varkappa=\eta\omega^{-1}.

The paper and Genera need not choose the same homology basis. The complete
script compares Genera's computed full periods with the rounded matrices in
the paper, recovers an integral symplectic cycle transformation, and then
forms the paper-basis ``omega``, ``omega_prime``, ``eta``, and ``eta_prime``.
This is a change of homology basis, not a fitted numerical scale factor.

Evaluating Kleinian P-functions
-------------------------------

For the Abelian vector and characteristic printed in Example 3a, the script
evaluates eight second- and third-order P-functions in one batched call::

   indices = (
       (0, 0), (0, 1), (0, 2), (1, 1), (1, 2),
       (0, 0, 0), (0, 0, 1), (0, 0, 2),
   )
   values = kleinian_p(
       u, omega / 2, tau, -kappa, indices, characteristic)

The factors in this call are important. Bernatska's tables use full period
matrices, while ``kleinian_p`` uses half-period data. The second-kind sign is
also reversed to match Genera's sigma convention. The same calculation is
then repeated with the modified Abel vector of Example 3b.

The resulting validation chain is::

   curve equation
       -> branch locus, monodromy, and homology
       -> first- and second-kind periods
       -> tau, kappa, and Kleinian P-functions

The complete script reports residuals for the branch points, period matrices,
symmetry conditions, and the P-function values against the six-digit tables
in the paper. The following directive runs that complete script during every
documentation build and displays its captured output. A nonzero exit status
fails the documentation build.

Validation output
-----------------

.. genera-example:: examples.bernatska.bernatska_trigonal_demo
   :caption: Output from the Bernatska validation example:
   :timeout: 120
