# Verifying Ônishi's arbitrary-genus determinants with Generapy

This example numerically evaluates both sides of Theorems 7.2 and 8.3 in Y.
Ônishi, *Determinant Expressions for Hyperelliptic Functions*
(arXiv:math/0105189).  These are the arbitrary-genus
Frobenius--Stickelberger formula and its confluent Kiepert-type limit for an
odd-degree hyperelliptic curve.  The genus-two and genus-three papers in this
directory provide useful special cases:

- `arXiv-math0105188v4`: genus two;
- `arXiv-math0105187v3`: genus three;
- `arXiv-math0105189v8`: arbitrary genus.

The runnable implementation is `onishi_determinant_demo.py`.  It produces
textual numerical comparisons and no plots.

## 1. The formula

Let

\[
 C:y^2=f(x)
\]

be a genus-\(g\) hyperelliptic curve, where \(f\) is monic of degree
\(2g+1\) with distinct roots.  Ônishi uses the holomorphic differentials

\[
 du_i=\frac{x^{i-1}dx}{2y},\qquad 1\leq i\leq g.
\]

For a point \(P_j=(x_j,y_j)\) on the curve, write

\[
 u^{(j)}=\int_\infty^{P_j}(du_1,\ldots,du_g)^T.
\]

These are one-point Abel images, so they lie on the first theta stratum.
Theorem 7.2 relates sigma functions at these Abel images to an \(n\times n\)
determinant of algebraic functions on the curve.

Ônishi defines

\[
 \natural^n=
 \{i:n+1\leq i\leq g,\ i\equiv n+1\pmod 2\}
\]

and

\[
 \sigma_{\natural^n}(u)=
 \left(\prod_{i\in\natural^n}\frac{\partial}{\partial u_i}\right)
 \sigma(u).
\]

Two derivatives occur throughout the formula:

\[
 \sigma_\sharp=\sigma_{\natural^1},
 \qquad
 \sigma_\flat=\sigma_{\natural^2}.
\]

If \(1\leq n<g\), the sigma side is

\[
 c_{g,n}
 \frac{
   \sigma_{\natural^n}(u^{(1)}+\cdots+u^{(n)})
   \prod_{i<j}\sigma_\flat(u^{(i)}-u^{(j)})
 }{
   \sigma_\sharp(u^{(1)})^n\cdots
   \sigma_\sharp(u^{(n)})^n
 }.
\]

If \(n\geq g\), the derivative in the numerator is replaced by the
undifferentiated sigma function:

\[
 c_{g,n}
 \frac{
   \sigma(u^{(1)}+\cdots+u^{(n)})
   \prod_{i<j}\sigma_\flat(u^{(i)}-u^{(j)})
 }{
   \sigma_\sharp(u^{(1)})^n\cdots
   \sigma_\sharp(u^{(n)})^n
 }.
\]

The other side is the determinant obtained by evaluating the first \(n\)
functions in the sequence

\[
 1,x,x^2,\ldots,x^g,y,x^{g+1},xy,x^{g+2},x^2y,\ldots
\]

at the \(n\) curve points.  This sequence is ordered by pole order at the
unique point at infinity.  Since

\[
 \operatorname{ord}_\infty(x^a)=-2a,
 \qquad
 \operatorname{ord}_\infty(x^b y)=-(2g+1+2b),
\]

the demo can generate the columns for any genus and determinant size rather
than storing separate genus-specific formulas.

## 2. Curve conventions

The Generapy hyperelliptic functions use

\[
 Y^2=P(x),\qquad
 dU_i=\frac{x^{i-1}dx}{Y}.
\]

To match Ônishi without rescaling the Abelian coordinates, the demo takes

\[
 P(x)=4f(x),\qquad Y=2y.
\]

Consequently,

\[
 \frac{x^{i-1}dx}{Y}
 =\frac{x^{i-1}dx}{2y}.
\]

The point passed to `curve.abel_map_kind_1` is therefore `(x, Y)`, while
the determinant uses `(x, y) = (x, Y/2)`.

For a simple reproducible real example, the script chooses

\[
 Y^2=4\prod_{k=-g}^{g}(x-k).
\]

There are \(2g+1\) distinct real branch points.  The evaluation points are
placed to their right and use the positive square root of \(P(x)\).  Complex
curves or points are not excluded by the formula; the real choice merely
makes the demonstration and its output easy to inspect.

## 3. Sigma derivatives in Generapy notation

`kleinian_sigma_jet` uses a zero-based derivative-count tuple.  For example,
in genus four the key

```python
(0, 1, 0, 1)
```

means

\[
 \frac{\partial^2\sigma}{\partial u_2\partial u_4}
\]

in Ônishi's one-based notation.  The helper `natural_index(genus, n)`
performs this translation directly from the definition of \(\natural^n\).
Some low-genus cases are:

| Genus | `sigma_sharp` | `sigma_flat` |
|---:|---|---|
| 1 | `(0,)` | `(0,)` |
| 2 | `(0, 1)` | `(0, 0)` |
| 3 | `(0, 1, 0)` | `(0, 0, 1)` |
| 4 | `(0, 1, 0, 1)` | `(0, 0, 1, 0)` |
| 5 | `(0, 1, 0, 1, 0)` | `(0, 0, 1, 0, 1)` |

A zero tuple selects sigma itself, so `sigma_flat` in genus two is the
undifferentiated sigma function.

Every evaluation uses

```python
normalization="hyperelliptic"
```

which selects the canonical Schur--Weierstrass normalization already used
by the Generapy Kleinian functions.  This is essential here because the
individual normalization constants do not cancel completely from the
formula.

## 4. Fixing the overall sign

The displayed theorem includes a sign \(c_n\).  Its direct numerical use is
convention-sensitive: the sign printed in the paper agrees with the current
Generapy normalization for some \((g,n)\), while for other pairs the two sides
are negatives of each other.  The conversion is not a single fixed sign for
each genus.

The demo does not determine the sign by evaluating the desired identity and
choosing whichever answer works.  Instead, it repeats the local calculation
which fixes the constant in Ônishi's proof.

Near infinity, use a parameter \(t=u_g\).  The leading expansions are

\[
 u_j(t)=\frac{t^{2(g-j)+1}}{2(g-j)+1},
 \qquad
 x(t)=t^{-2},
 \qquad
 y(t)=-t^{-(2g+1)}.
\]

The leading term of the canonically normalized sigma function is the
Schur--Weierstrass polynomial.  The script constructs that polynomial
exactly with rational coefficients from

\[
 S(u)=\det\left(U_{g-2i+j+1}(u)\right)_{1\leq i,j\leq g},
\]

where the coefficients \(U_k\) are generated by

\[
 \exp\left(
   \sum_{j=1}^{g}u_j z^{2(g-j)+1}
 \right)=\sum_{k\geq0}U_k(u)z^k.
\]

Its scale is fixed by the same coefficient condition used internally by
Generapy's hyperelliptic sigma normalization.  Substituting several distinct
rational values of \(t\) into the leading sigma quotient and leading
determinant gives a ratio of exactly \(+1\) or \(-1\).  The inverse of that
ratio is the Generapy prefactor.

Each numerical check reports:

```text
printed c_n             sign displayed in Theorem 7.2
Generapy prefactor        sign from the exact local calculation
convention conversion   ratio between those signs
```

The full sigma identity is evaluated only afterwards.  No absolute values
are used when comparing its two sides, so an unresolved sign error would
remain visible as an order-one residual.

## 5. Structure of the script

The main pieces of `onishi_determinant_demo.py` are:

- `curve_coefficients`: constructs the odd-degree demonstration curve;
- `natural_index`: translates Ônishi's special derivatives;
- `onishi_monomials`: generates determinant columns in pole order;
- `schur_weierstrass_polynomial`: constructs the exact leading sigma term;
- `generapy_prefactor`: derives the convention-adjusted sign locally;
- `evaluate_formula`: independently evaluates the two sides of Theorem 7.2;
- `curve_derivative`: differentiates \(A(x)+B(x)Y\) algebraically;
- `generapy_kiepert_prefactor`: derives the confluent determinant sign;
- `evaluate_kiepert_formula`: evaluates both sides of Theorem 8.3.

The numerical path is:

1. Construct an algebraic curve and obtain `omega`, `tau`, `kappa`, and the
   characteristic from its first- and second-kind period records.
2. Select affine points `(x, Y)` on the curve.
3. Compute every one-point Abel image with `curve.abel_map_kind_1`.
4. Evaluate `sigma_sharp`, `sigma_flat`, and the numerator through
   `kleinian_sigma` or `kleinian_sigma_jet`.
5. Construct and evaluate the algebraic determinant.
6. Print both complex values and their absolute and relative residuals.

The Abel images are not reduced to a fundamental cell.  They retain the
deterministic integration paths from infinity used by `curve.abel_map_kind_1`,
giving consistent representatives on the universal
cover on which the individual sigma factors are defined.

## 6. Running the example

From the Generapy repository root, run the default genus-three cases with:

```bash
.venv/bin/python -m examples.onishi.onishi_determinant_demo
```

For Theorem 7.2, the defaults test:

- \(n=2<g\), using `sigma_natural^2`;
- \(n=g=3\), using sigma at the boundary between the two cases;
- \(n=5\), whose determinant has reached the first \(y\)-column.

The default Kiepert checks use \(n=g\) and \(n=g+2\), differentiating with
respect to both \(u_1\) and \(u_g\).

Specify a genus, determinant sizes, and precision with:

```bash
.venv/bin/python -m examples.onishi.onishi_determinant_demo \
    --genus 4 --n 2 3 4 6 --kiepert-n 4 6 --coordinate 1 4 --dps 25
```

The formula is implemented without a fixed genus limit.  Riemann theta
evaluation grows exponentially with genus, however, and construction of the
exact Schur polynomial also becomes more expensive.  “Arbitrary genus” here
means that the mathematical construction is general, not that its runtime is
independent of genus.

## 7. Example output

At 30 decimal places, the default run includes:

```text
n = 2 (sigma_natural^2 numerator)
  printed c_n:            +1
  Generapy prefactor:       +1
  convention conversion: +1
  sigma side:             (1.0 + 0.0j)
  determinant side:       1.0
  relative residual:      3.9443e-30

n = 3 (sigma numerator)
  printed c_n:            +1
  Generapy prefactor:       -1
  convention conversion: -1
  sigma side:             (2.0 + 0.0j)
  determinant side:       2.0
  relative residual:      1.04524e-29

n = 5 (sigma numerator)
  printed c_n:            -1
  Generapy prefactor:       +1
  convention conversion: -1
  sigma side:             (28.98223816436559 + 0.0j)
  determinant side:       28.98223816436559
  relative residual:      2.19057e-28
```

Separate checks in generapy one, two, and four also agree to their working
precision.  The genus-four \(n=6\) case exercises a second-order
`sigma_sharp`, a nontrivial `sigma_flat`, and a determinant beyond the pure
Vandermonde columns.

The same run now continues with the Kiepert-type output.  For example:

```text
n = 5, derivative coordinate j = 3
  printed c'_n:           -1
  Generapy prefactor:       +1
  convention conversion: -1
  sigma side:             (2.826567837490959e+27 + 0.0j)
  determinant side:       2.826567837490959e+27
  relative residual:      7.34175e-29
```

## 8. Relation to the rest of the branch

This example uses only public functionality already present on the
Kleinian branch:

- `Curve.periods_kind_1` and `periods_kind_2` for coherent
  period and sigma data;
- `Curve.abel_map_kind_1` for points on the one-point stratum;
- `kleinian_sigma` for the \(n\geq g\) numerator;
- `kleinian_sigma_jet` for all special derivatives.

It therefore tests these APIs together in two formulas whose algebraic sides
are computed independently from curve coordinates.  In contrast with
examples based on a dynamical system, no ODE integration or fitted physical
initial condition is involved.

The sign handling remains local to the example.  It does not alter Generapy's
sigma normalization, period basis, Abel-map paths, or public API.

## 9. The Kiepert-type confluent determinant

Ônishi's Theorem 8.3 takes the confluent limit of Theorem 7.2.  Define the
division function

\[
 \psi_n(u)=
 \begin{cases}
 \displaystyle
 \frac{\sigma_{\natural^n}(nu)}{\sigma_\sharp(u)^{n^2}},&n<g,\\[8pt]
 \displaystyle
 \frac{\sigma(nu)}{\sigma_\sharp(u)^{n^2}},&n\geq g.
 \end{cases}
\]

For the determinant case \(n\geq g\), fix an Abelian coordinate
\(u_j\), \(1\leq j\leq g\).  Remove the initial constant from the
pole-ordered monomial sequence and call the next \(n-1\) functions
\(\phi_1,\ldots,\phi_{n-1}\).  The formula is

\[
 c'_{g,n}\left(\prod_{r=1}^{n-1}r!\right)\psi_n(u)
 =x(u)^{(j-1)n(n-1)/2}
 \det\left(
   \frac{d^r\phi_s}{du_j^r}(u)
 \right)_{1\leq r,s\leq n-1}.
\]

Thus distinct rows of point evaluations in Theorem 7.2 become derivative
orders \(1,\ldots,n-1\) at one point.

### Differentiation along the curve

From

\[
 du_j=\frac{x^{j-1}dx}{Y},
\]

the curve derivative is

\[
 \frac{d}{du_j}=\frac{Y}{x^{j-1}}\frac{d}{dx}.
\]

Repeated numerical differentiation would introduce an avoidable second
source of error.  Instead, the script represents every determinant entry as

\[
 A(x)+B(x)Y,
\]

where \(A\) and \(B\) are Laurent polynomials.  Since \(Y^2=P(x)\),

\[
 \frac{d}{du_j}(A+BY)
 =x^{1-j}\left(B'P+\frac12BP'\right)
  +x^{1-j}A'Y.
\]

This recurrence evaluates all derivative orders algebraically.  It handles
both kinds of monomial column, introduces no finite-difference step, and is
valid for every permitted coordinate \(j\).

### Kiepert sign conversion

The Kiepert prefactor is derived independently by the same exact local
Schur calculation used for Theorem 7.2.  Locally,

\[
 \frac{d}{du_j}=t^{1-w_j}\frac{d}{dt},
 \qquad w_j=2(g-j)+1.
\]

Applying this operator exactly to the leading monomials
\(c t^k\), and comparing the resulting rational determinant with the
leading Schur expression for \(\psi_n\), produces the Generapy sign
\(c'_{g,n}\).  This calculation is performed before evaluating the full
theta functions.  The output again shows the printed sign, the Generapy sign,
and their ratio.

The command-line options are:

```text
--kiepert-n N [N ...]       multipliers, each at least max(2, genus)
--coordinate J [J ...]     one-based Abelian derivative coordinates
```

For example:

```bash
.venv/bin/python -m examples.onishi.onishi_determinant_demo \
    --genus 3 --n 2 3 5 --kiepert-n 3 5 --coordinate 1 3 --dps 30
```

In genus three this checks both \(d/du_1\) and \(d/du_3\).  Although the
two determinants have very different individual entries, the compensating
power of \(x\) makes both reproduce the same sigma division function.  For
\(n=3\) and \(n=5\), respectively, the observed relative residuals are
approximately \(3.2\times10^{-30}\) and \(7.4\times10^{-29}\).

## References

- Y. Ônishi, *Determinant Expressions for Abelian Functions in Genus Two*,
  arXiv:math/0105188.
- Y. Ônishi, *Determinant Expressions for Hyperelliptic Functions in Genus
  Three*, arXiv:math/0105187.
- Y. Ônishi, *Determinant Expressions for Hyperelliptic Functions*,
  arXiv:math/0105189, especially Definitions 2.1 and 6.1 and Theorem 7.2.
- V. M. Buchstaber, V. Z. Enolskii and D. V. Leykin, *Hyperelliptic
  Kleinian Functions and Applications*, for the sigma convention used by
  the Generapy implementation.
