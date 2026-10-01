# Bernatska's non-hyperelliptic genus-three example

This example reproduces the trigonal calculations in Example 3 and Example
3a of J. Bernatska, *Computation of P-Functions on Plane Algebraic Curves*
(arXiv:2407.05632v4).  It is the example that was previously out of reach of
Genera's hyperelliptic-only curve-data machinery.

The runnable, no-plot script is `bernatska_trigonal_demo.py`.  From the Genera
repository root, run

```text
.venv/bin/python -m examples.bernatska.bernatska_trigonal_demo
```

Genera and its mpmath dependency are sufficient at runtime.

## 1. The curve

Bernatska's equation (73) is the genus-three `(3,4)` curve

\[
 f(x,y)=-y^3+x^4+3x^3+7x^2+16x+9
        +y(4x^2+5x+11)=0.
\]

Writing

\[
 P(x)=x^4+3x^3+7x^2+16x+9,
 \qquad Q(x)=4x^2+5x+11,
\]

the finite branch values are the eight roots of

\[
 27P(x)^2-4Q(x)^3=0.
\]

The script passes the sparse coefficients of `f(x,y)` directly to
`algebraic_curve_data`.  The routine forms the projection resultant, removes
repeated polynomial factors, and finds the finite branch values internally.
The values printed in the paper are used only as an independent check, not as
input to the computation.  The complete curve-data call is

```python
first_kind, second_kind = differentials()
data = algebraic_curve_data(
    curve=trigonal_curve(),
    differentials_kind_1=first_kind,
    differentials_kind_2=second_kind,
)
```

No branch points, path method, homology cycles, or theta characteristic are
passed to the period calculation.

## 2. Complex monodromy

The branch values consist of two real points and three conjugate pairs.  The
curve-data routine selects the general radial coordinator and chooses a
well-separated exterior base point from a deterministic set of candidates.
Each finite generator is a positive
lollipop loop: a straight spoke approaches one branch value, a small
counter-clockwise circle encloses it, and the spoke is retraced.  Circle
radii account for distances to the other branch values and spokes.

The spokes are ordered clockwise at the base point.  Infinity is represented
by an actual clockwise outer circle in the affine `x`-plane; clockwise is
positive around the local coordinate `1/x`.  Numerical sheet continuation
then finds eight transpositions at finite branch values and a three-cycle at
infinity.

This distinction matters.  The real-axis coordinator developed first for
the Kovalevskaya example can infer an infinity loop by multiplying ordered
real branch loops.  The radial coordinator retains continuation-product and
outward ribbon orders separately and supplies the outer loop explicitly.

The lifted ribbon graph reports

```text
genus = 3
intersection rank = 6
boundary components = 3
```

and its integral reduction produces three canonical `a,b` pairs.

## 3. Differentials

Since

\[
 f_y=-3y^2+Q(x),
\]

the paper's holomorphic differentials are supplied directly:

\[
 du=
 \begin{pmatrix}y\\x\\1\end{pmatrix}\frac{dx}{f_y}.
\]

This follows the project decision that general plane curves initially accept
their differential basis from the caller rather than requiring symbolic
normalization.

The script also supplies Bernatska's second-kind basis

\[
 dr=
 \begin{pmatrix}x^2\\2xy\\R_5\end{pmatrix}\frac{dx}{f_y},
\]

where, for this curve,

\[
 R_5=5x^2y+9xy+\frac{32}{3}x^2+7y+\frac{40}{3}x.
\]

The returned object contains `omega`, `omega_prime`, `tau`, `eta`,
`eta_prime`, and `kappa`.  Its `diagnostics` member retains the automatically
computed branch values, monodromy, graph and raw full-period matrices used for
the paper-basis comparison below.

## 4. Comparing homology bases

A numerical homology basis is not expected to equal the basis drawn in the
paper.  Comparing the realified first-kind period matrices finds the integer
cycle transformation

```text
 1  1 -1  1  1 -1
 1 -1 -1  0  0 -1
 0  0  0  1  0  1
 0  1  0  1  1  0
 0 -1  0  0 -1  0
 0  0 -1  1  0  0
```

for the current deterministic paths.  Its determinant is one and it preserves
the standard intersection matrix exactly, so it is an integral symplectic
lattice change rather than a sublattice.  After the transformation, the
first-kind periods agree with the six-digit printed matrices.  The same exact
cycle transformation is applied to the independently integrated second-kind
periods.

The resulting normalized matrix

\[
 \tau=\omega^{-1}\omega'
\]

is symmetric and has positive-definite imaginary part.  Likewise

\[
 \varkappa=\eta\omega^{-1}
\]

is symmetric.  These are internal checks that do not rely only on agreement
with rounded tables in the paper.

## 5. Examples 3a and 3b P-functions

Bernatska gives the Abel image

\[
 u(D)\approx
 \begin{pmatrix}
 -0.270333-1.612257i\\
 -1.116879+0.562199i\\
  0.258194+0.268653i
 \end{pmatrix}
\]

and the Riemann-constant characteristic

\[
 [K]=\begin{pmatrix}1&0&1\\0&1&1\end{pmatrix}.
\]

The entries of the latter become half-integers in Genera's literal theta
characteristic convention.  Bernatska uses full periods, whereas
`kleinian_p` takes the classical half-period matrix, so the script passes
`omega/2`.  Her displayed `varkappa` has the opposite sign from Genera's
sigma convention, just as in the existing genus-four Bernatska test.

With those explicit conversions, `kleinian_p` evaluates the eight second-
and third-order functions printed in equation (86).  The script repeats the
calculation with the modified divisor and Abel vector of Example 3b, checking
the eight values in equation (91) as an independent second evaluation.  This
closes the chain

```text
curve equation
  -> complex monodromy
  -> canonical homology
  -> first- and second-kind periods
  -> tau and kappa
  -> trigonal Kleinian P-functions.
```

The characteristic used for the printed P-function comparison is taken from
the paper because its Abel vectors, characteristic and periods all use the
paper's point at infinity and homology basis.  The period-data API now also
computes a characteristic automatically, but it belongs to the returned
finite `data.base_place` and the automatically constructed cycle basis.  It
must not be mixed directly with the paper's Abel vectors.  Passing a regular
finite `(x, y)` as `base_place` changes the returned characteristic by the
corresponding normalized Abel shift.

The Abel vector is taken from the paper in this first demonstration.  A
future extension can reconstruct its three marked integration paths using
the local-chart and marked-path machinery already exercised by the
Kovalevskaya example.
