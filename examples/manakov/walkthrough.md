# Genus-two Manakov-type solution: walkthrough

This note explains `manakov_kleinian_demo.py` end to end: the stationary
Manakov-type system, its genus-two Kleinian solution, the convention mapping,
the numerical choices and the independent checks. The reference is
Christiansen, Eilbeck, Enolskii and Kostov,
[_Quasi-periodic and periodic solutions for coupled nonlinear Schrodinger
equations of Manakov type_](https://arxiv.org/abs/solv-int/9904017),
Proc. R. Soc. Lond. A 456 (2000), 2263–2281 (hereafter CEEK).

## 1. The problem

The demo solves the two-degree-of-freedom stationary system in CEEK equation
(2.1):

$$
q_i''=-(q_1^2+q_2^2)q_i+a_iq_i+\frac{C_i^2}{q_i^3},
\qquad i=1,2.
$$

Writing $p_i=q_i'$ gives the four-dimensional first-order system integrated
independently by RK4. Its conserved Hamiltonian, CEEK equation (2.2), is

$$
H=\frac{p_1^2+p_2^2}{2}
 +\frac{(q_1^2+q_2^2)^2}{4}
 -\frac{a_1q_1^2+a_2q_2^2}{2}
 +\frac{1}{2}\left(\frac{C_1^2}{q_1^2}
                   +\frac{C_2^2}{q_2^2}\right).
$$

The comparison is deliberately end to end: Genera constructs the spectral
curve's periods, evaluates the analytic Kleinian solution, converts that
solution back to a spectral divisor and integrates the ODE without using the
analytic values again.

## 2. Curve and period data

The demonstration curve is

$$
y^2=4x(x-\tfrac18)(x-\tfrac14)(x-\tfrac12)(x-1).
$$

It is nonsingular, has five real branch points and has genus two. The marked
parameters are $a_1=3/4$ and $a_2=3/16$. The constants $C_i^2$ are obtained
from the curve polynomial exactly as required by CEEK:

```python
coefficients = curve_coefficients()
separation = a1 - a2
c1_squared = -polynomial_value(coefficients, a1) / separation**2
c2_squared = -polynomial_value(coefficients, a2) / separation**2
curve = mp.algebraic_curve(coefficients)
first = curve.periods_kind_1()
second = curve.periods_kind_2()
omega, tau = first.omega, first.tau
kappa = second.kappa
characteristic = curve.riemann_constant().characteristic
```

The two curve records provide first-kind half-periods $\omega$, the normalized
Riemann matrix $\tau$, $\varkappa=\eta\omega^{-1}$ and the characteristic of
the vector of Riemann constants in one coherent homology basis. CEEK also use
the classical half-period convention, so no factor-of-two conversion is
needed at the Kleinian interface.

## 3. Analytic Kleinian solution

The Genera differential basis is ordered as

$$
(du_1,du_2)^T=(dx/y,x\,dx/y)^T.
$$

Consequently CEEK's one-based quantities map to the zero-based Genera API as

| CEEK        | `kleinian_p` indices |
| ----------- | -------------------- |
| $\wp_{22}$  | `(1, 1)`             |
| $\wp_{12}$  | `(0, 1)`             |
| $\wp_{222}$ | `(1, 1, 1)`          |
| $\wp_{122}$ | `(0, 1, 1)`          |

With $d=a_1-a_2$, equation (3.23) becomes

$$
q_1^2=\frac{2(a_1^2-a_1\wp_{22}-\wp_{12})}{d},\qquad
q_2^2=\frac{2(a_2^2-a_2\wp_{22}-\wp_{12})}{-d},
$$

and differentiation gives

$$
p_1=-\frac{a_1\wp_{222}+\wp_{122}}{d q_1},\qquad
p_2= \frac{a_2\wp_{222}+\wp_{122}}{d q_2}.
$$

The positive square-root branches are continuous on the selected real
trajectory. A half-period translate
$\omega(1,1)^T+\omega\tau(1,1)^T$ selects a component on which both
$q_i^2$ stay positive over the comparison interval. The first Abelian coordinate
is then varied through the three default phases `-6`, `0.13` and `6`; the
second advances linearly with the independent variable.

CEEK print $u_2=2x+b$ below equation (2.12), while equations (3.18)–(3.19)
identify $x=u_2$. Direct substitution into equation (2.1), using the paper's
curve and differential normalization, selects $u_2=x+b$; the doubled version
leaves order-one residuals.

## 4. Recovering the divisor and Abelian phase

The physical state contains enough information to reverse the construction.
Define

$$
U(\lambda)=\lambda^2-\wp_{22}\lambda-\wp_{12},\qquad
V(\lambda)=\wp_{222}\lambda+\wp_{122}.
$$

The two formulas for $q_i^2$ give the values of
$a_i\wp_{22}+\wp_{12}$, so a two-by-two linear system recovers
$\wp_{22}$ and $\wp_{12}$. Similarly, $q_ip_i$ recovers $\wp_{222}$ and
$\wp_{122}$. This uses only the physical state and fixed parameters; it does
not perform another theta evaluation.

The roots $x_j$ of $U$ are the divisor's x-coordinates. CEEK equations
(3.10)–(3.12) give their sheets as

$$
y_j=V(x_j)=\wp_{222}x_j+\wp_{122}.
$$

Thus the phase is reconstructed by

```python
divisor = divisor_from_state(initial_state, data)
image = data["curve"].abel_map_kind_1(divisor, reduce=True)
```

The returned representative differs from the original Abelian point by a
full period only. Re-evaluating the Kleinian formulas at that representative
recovers all four physical coordinates. This checks the state-to-divisor
algebra, both sheet signs, the Abel-map integration paths, the period basis
and the fact that no extra Riemann-constant translation is required.

## 5. Independent checks

Four checks accompany the numerical solution:

1. **RK4 comparison.** Classical fourth-order Runge–Kutta evolves equation
   (2.1) from the analytic initial state. Step halving shows the expected
   fourth-order convergence.
2. **Hamiltonian conservation.** The numerical trajectory conserves equation
   (2.2), independently of later analytic evaluations.
3. **Kleinian identities.** CEEK equations (3.16) and (3.17) check the
   fourth-order derivatives. In the sigma normalization used here the curve
   constant in (3.16) enters as `coefficients[3]/2`; this agrees with the
   classical genus-one Weierstrass limit.
4. **Abel-map closure.** The physical state is converted to its divisor and
   mapped back to its original Jacobian point modulo full periods, as
   described above.

At 30 decimal places, with 3200 RK4 steps on $[-4,4]$, the original run gives
position errors below $4.2\times10^{-11}$ and Hamiltonian drift below
$3.7\times10^{-12}$. The phase `0.13` Abel-map check gives residuals at the
working-precision scale:

```text
equation (3.16) residual:             1.8e-31
equation (3.17) residual:             1.2e-32
Abel-map divisor curve residual:      1.4e-32
Abel-map period-lattice residual:     3.5e-30
Abel-map recovered-state residual:    6.9e-31
```

## 6. Running the demo

From the repository root:

```bash
.venv/bin/python -m examples.manakov.manakov_kleinian_demo
```

Command-line options select the precision, integration interval, RK4 step
count, sampling count, phases and acceptance tolerance. The script prints the
maximum component errors and exits unsuccessfully if they exceed the selected
tolerance.


The demo covers the stationary real reduction (2.1), not the complete complex
Manakov fields of equation (3.23). Those fields additionally require the
Baker–Akhiezer phase integrals and continuous square-root phase tracking.
