# Nine-dimensional Reissner-Nordstrom-de Sitter geodesics: walkthrough

This note explains `rn_desitter_geodesic_demo.py` end to end: the radial geodesic equation, its reduction to a genus-three hyperelliptic curve, the theta-divisor inversion, the numerical continuation and the independent RK4 check. The reference is Enolski, Hackmann, Kagramanova, Kunz and Lammerzahl, [_Inversion of hyperelliptic integrals of arbitrary genus with application to particle motion in General Relativity_](https://arxiv.org/abs/1011.6459), Journal of Geometry and Physics 61 (2011), 899--921 (hereafter EHKL).

## 1. The geodesic problem

The example follows EHKL section 7. It uses a massive test particle ($\delta=1$) in the nine-dimensional Reissner-Nordstrom-de Sitter spacetime, with the dimensionless parameters

$$
\widetilde\Lambda=8.7\times10^{-5},\qquad
\widetilde q=0.4,\qquad E^2=1.045,\qquad
\widetilde L^2=0.25.
$$

Dropping tildes in the implementation, the radial equation is

$$
\left(\frac{dr}{d\varphi}\right)^2=\frac{R_{16}(r)}{r^{10}},
$$

where

$$
L^2R_{16}(r)=\frac{\Lambda r^{16}}{28}
 +\left(E^2-1+\frac{L^2\Lambda}{28}\right)r^{14}
 -L^2r^{12}+r^8+L^2r^6-q^2r^2-q^2L^2.
$$

The three positive zeros of the metric function are

| Horizon | Radius |
| --- | ---: |
| Cauchy | $0.76472441$ |
| Event | $0.96349310$ |
| Cosmological | $567.30863$ |

The cosmological horizon is much farther out than either numerical
trajectory, so the reported radial ranges focus on the two local horizons
and state the cosmological radius separately.

## 2. Reduction to a genus-three curve

Set $u=r^{-2}$. The radial equation becomes

$$
\left(\frac{du}{d\varphi}\right)^2=\frac{4}{L^2}R_8(u),
$$

with

$$
R_8(u)=\frac{\Lambda}{28}
 +\left(E^2-1+\frac{L^2\Lambda}{28}\right)u-L^2u^2
 +u^4+L^2u^5-q^2u^7-L^2q^2u^8.
$$

Choose a root $u_8$ of $R_8$ and make the second substitution

$$
u=\frac1x+u_8.
$$

The selected root is $u_8\simeq-0.6073545616$. It sends one finite branch point to infinity and produces the odd-degree genus-three curve

$$
y^2=P_7(x)=\frac{4}{L^2}x^8R_8\left(\frac1x+u_8\right).
$$

Since

$$
\left(x^2\frac{dx}{d\varphi}\right)^2=P_7(x),
$$

the azimuthal angle is the third holomorphic Abelian integral,

$$
\varphi-\varphi_{\mathrm{in}}=\int\frac{x^2\,dx}{y}.
$$

This is why the physical flow advances the last Abelian coordinate in the demo.

The chosen $u_8$ gives five real finite branch points and one complex conjugate pair. The bound orbit runs between two adjacent real turning points. It starts at $r=0.74276301$, reaches $r=1.7748895$, and returns after the radial period $\Delta\varphi=2.875190201$. The escape orbit starts at $r=2.0643116$ and moves outwards; the comparison segment ends at about $r=6.986389$.

## 3. Inversion on the theta divisor

For the differential basis used by both EHKL and Genera,

$$
(du_1,du_2,du_3)^T=(dx/y,x\,dx/y,x^2\,dx/y)^T,
$$

the genus-three inversion formula is

$$
x=-\frac{\sigma_{13}(\boldsymbol u)}
         {\sigma_{23}(\boldsymbol u)}.
$$

A single point of the curve maps to the first stratum of the theta divisor. Starting from the Abel image $A_i$ of a turning branch point, the required Abelian point is therefore written

$$
\boldsymbol u(\varphi)=A_i+
  \begin{pmatrix}f_1(\varphi)\\f_2(\varphi)\\
                  \varphi-\varphi_{\mathrm{in}}\end{pmatrix},
$$

where $f_1$ and $f_2$ are fixed by

$$
\sigma(\boldsymbol u)=0,\qquad
\sigma_3(\boldsymbol u)=0.
$$

The Genera coordinates are zero based, so the derivatives used in the code are

| EHKL quantity | `kleinian_sigma_jet` key |
| --- | --- |
| $\sigma$ | `(0, 0, 0)` |
| $\sigma_3$ | `(0, 0, 1)` |
| $\sigma_{13}$ | `(1, 0, 1)` |
| $\sigma_{23}$ | `(0, 1, 1)` |

These are ordinary sigma derivatives. Logarithmic Kleinian $\wp$-functions are singular on the theta divisor and are therefore not the right interface for this inversion.

### Why the reciprocal quotient appears in the radius

EHKL equations (4.13) and (6.18) give $x=-\sigma_{13}/\sigma_{23}$, while equation (7.5) contains the reciprocal derivative quotient inside the formula for $r$. This initially looks like a transposition. The substitution resolves it directly:

$$
r^{-2}=u=\frac1x+u_8
 =u_8-\frac{\sigma_{23}}{\sigma_{13}},
$$

and hence

$$
r(\varphi)=
\frac{1}{\sqrt{u_8-\sigma_{23}(\boldsymbol u)/
                         \sigma_{13}(\boldsymbol u)}}.
$$

Thus the two printed quotients serve different roles and are both correct. The demo also evaluates the genuinely transposed radial candidate. Its discrepancy is order one, whereas the implemented expression agrees with RK4 to between $10^{-8}$ and $10^{-7}$.

## 4. Period data and the starting point

The curve data are constructed in one call:

```python
curve = mp.algebraic_curve(p7)
first = curve.periods_kind_1()
second = curve.periods_kind_2()
omega, tau = first.omega, first.tau
kappa = second.kappa
characteristic = curve.riemann_constant().characteristic
```

Here `omega` is the first-kind **half-period** matrix; the complete periods are `2*omega` and `2*omega*tau`. The sigma convention used by `kleinian_sigma_jet` is already consistent with the curve records, so no factor-of-two conversion is needed inside the demo.

For each orbit, `curve.abel_map_kind_1` computes the reduced Abel image of its starting branch point:

```python
start = curve.abel_map_kind_1((x_start, 0), reduce=True)
```

The branch-point identities provide a sensitive convention check:

$$
e_i=-\frac{\sigma_1(A_i)}{\sigma_2(A_i)}
   =-\frac{\sigma_{13}(A_i)}{\sigma_{23}(A_i)}.
$$

Both identities hold at the working-precision scale for the bound and escape starting points. This simultaneously checks the Abel-map path, sheet, coordinate order, period basis and absence of an extra Riemann-constant translation.

## 5. Numerical continuation

At every requested angle the script solves the two stratum constraints for $(f_1,f_2)$ by Newton iteration. One order-two sigma jet supplies the two residuals, the Jacobian and the inversion quotient.

After a point has converged, implicit differentiation gives the tangent used to predict the next point:

$$
\begin{pmatrix}
\sigma_1&\sigma_2\\
\sigma_{13}&\sigma_{23}
\end{pmatrix}
\begin{pmatrix}f_1'\\f_2'\end{pmatrix}
=-
\begin{pmatrix}\sigma_3\\\sigma_{33}\end{pmatrix}.
$$

The first step after a branch point uses an ordinary warm start because the branch-point local parameter is singular. Later points use the tangent prediction followed by Newton correction. This reduces the default run from roughly 25 minutes to about 70 seconds on the development machine.

The default numerical settings are:

| Setting | Value |
| --- | ---: |
| Working precision | 25 decimal places |
| Kleinian samples per orbit | 31 |
| Stratum residual tolerance | $10^{-18}$ |
| RK4 steps per orbit | 3000 |

Near the end of the escape segment, $1/x+u_8$ is a small difference, so it amplifies harmless imaginary roundoff in the reconstructed radius. The code checks that this remains below $10^{-7}$ relative before projecting back to the real slice. A wrong branch or quotient instead produces an order-one imaginary part or discrepancy.

## 6. Independent RK4 comparison

Define

$$
V(r)=\frac{R_{16}(r)}{r^{10}}.
$$

Differentiating $(r')^2=V(r)$ gives the second-order equation

$$
r''=\frac12V'(r),
$$

which the demo integrates as a first-order system using classical RK4, with $r'(0)=0$ at each turning point. This calculation uses the original radial potential and does not reuse the sigma solution.

The RK4 solver stores a dense 3001-point trajectory. The 31 Kleinian values
are compared at corresponding integer grid indices, avoiding equality tests
between independently rounded arbitrary-precision times.

For the default run:

```text
sigma-jet evaluations:                  1411
bound worst stratum residual:           7.58e-19
bound max |r Kleinian - r RK4|:         2.20e-8
escape worst stratum residual:          8.90e-19
escape max |r Kleinian - r RK4|:        2.11e-7
bound transposed-quotient discrepancy:  3.45
escape transposed-quotient discrepancy: 5.98
elapsed time:                           about 69 seconds
```

The script reports the relative radial errors, the two stratum residuals and
the perihelion shift. A many-world orbit crosses horizon regions belonging to
different blocks of the maximally extended spacetime; the example checks its
radial equation and does not attempt a coordinate-plane visualization.

## 7. Running the demo

From the repository root:

```bash
.venv/bin/python -m examples.rn_desitter_9d.rn_desitter_geodesic_demo
```

The useful command-line options are `--dps`, `--samples`, `--steps`,
`--newton-tol` and `--tol`. The RK4 step count must be divisible by
`samples - 1` so every Kleinian sample has an exact comparison index. The
script exits unsuccessfully when the analytic-versus-RK4 radial error exceeds
the selected tolerance.
