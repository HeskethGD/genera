# Solving the genus-three Neumann-Moser system with Genera Kleinian functions

This note walks through the genus-three demo in this folder (`neumann_moser_genus_three_demo.py`) end to end. It parallels `walkthrough.md`, the genus-two companion: the printed-formula corrections, the index convention and the verification strategy established there carry over verbatim, so this note spells out only the genuinely new genus-three material. The reference is again P. G. Baron, [arXiv:2402.18079](https://arxiv.org/abs/2402.18079), with the solution theory quoted from V. M. Buchstaber, [arXiv:2402.09218](https://arxiv.org/abs/2402.09218); part 18 records the corrections and part 19 the Abel-map closure.

## 1. The system

In the dimension $n = g = 3$ the system has ten coordinates $(u_1, u_2, u_3, v_1, v_2, v_3, w_1, w_2, w_3, w_4)$ and, with $\Gamma = w_1 - u_1$:

$$
\begin{aligned}
\dot u_1 &= -2v_1, & \dot v_1 &= -\Gamma u_1 - u_2 + w_2, & \dot w_1 &= 2v_1, \\
\dot u_2 &= -2v_2, & \dot v_2 &= -\Gamma u_2 - u_3 + w_3, & \dot w_2 &= 2v_2 + 2\Gamma v_1, \\
\dot u_3 &= -2v_3, & \dot v_3 &= -\Gamma u_3 + w_4, & \dot w_3 &= 2v_3 + 2\Gamma v_2, \\
& & & & \dot w_4 &= 2\Gamma v_3.
\end{aligned}
$$

The repairs to the paper's expanded equations carry over unchanged: the middle $\dot w$ range runs to $n$ and the last equation is $\dot w_{n+1} = 2\Gamma v_n$. The system is the 5-stationary KdV hierarchy, with $\Gamma = 2\wp_2$ its KdV solution.

## 2. The curve and the Genera data

The canonical curve is now of degree seven, with no $x^6$ term so that the seven branch points sum to zero:

$$ y^2 = F(x) = 4x^7 + \lambda_4 x^5 + \lambda_6 x^4 + \lambda_8 x^3 + \lambda_{10} x^2 + \lambda_{12} x + \lambda_{14}. $$

The demo picks seven distinct real roots summing to zero, $(-5, -3, -1.2, -0.3, 0.8, 2.7, 6)$, giving $\lambda_4 = -158.92$ and $\lambda_6 = -109.92$.

```python
from genera import kleinian_p
from mpmath import mp
mp.dps = 30

roots = (mp.mpf(-5), mp.mpf(-3), mp.mpf("-1.2"), mp.mpf("-0.3"),
         mp.mpf("0.8"), mp.mpf("2.7"), mp.mpf(6))
coefficients = [mp.mpf(4)]          # 4 * prod(x - root)
for root in roots:
    coefficients = multiply_by_linear(coefficients, root)

curve = Curve({(0, 2): 1, **{(i, 0): -c for i, c in enumerate(coefficients)}})
first = curve.periods_kind_1()
second = curve.periods_kind_2()
omega, tau = first.omega, first.tau
kappa = second.kappa
characteristic = curve.riemann_constant().characteristic
```

## 3. The Kleinian solution

The solution still flows along $z_1$, with $' = d/dz_1$:

$$
u_i = -\wp_{2i}, \qquad v_i = \tfrac12 \wp'_{2i}, \qquad W_\xi = (\xi + 2\wp_2)\,p_I + \tfrac12 p_{III},
$$

where $p_I = \xi^3 - \wp_2\xi^2 - \wp_4\xi - \wp_6$ and $p_{III} = \wp''_2\xi^2 + \wp''_4\xi + \wp''_6$. Expanding $W_\xi$ gives

$$
w_1 = \wp_2, \qquad w_2 = \tfrac12\wp''_2 - \wp_4 - 2\wp_2^2, \qquad w_3 = \tfrac12\wp''_4 - \wp_6 - 2\wp_2\wp_4, \qquad w_4 = \tfrac12\wp''_6 - 2\wp_2\wp_6.
$$

The corrected $p_I$ factor (the paper prints $p_{II}$, which cannot match the degree-$(g+1)$ generating polynomial) is the same repair as in genus two.

**Index convention.** The weight rule $z_{2k-1} \leftrightarrow$ Genera coordinate $g - k$ now reads $z_1/z_3/z_5 \to$ coordinates $2/1/0$:

| Paper | `kleinian_p` indices |
| --- | --- |
| $\wp_2, \wp_4, \wp_6$ | `(2, 2)`, `(1, 2)`, `(0, 2)` |
| $\wp_{3,3}, \wp_{3,5}$ | `(1, 1)`, `(1, 0)` |
| $\wp'_2, \wp'_4, \wp'_6$ | `(2,2,2)`, `(1,2,2)`, `(0,2,2)` |
| $\wp''_2, \wp''_4, \wp''_6$ | `(2,2,2,2)`, `(1,2,2,2)`, `(0,2,2,2)` |

The flow advances the **last** coordinate, the genus-three analogue of the second coordinate in the genus-two demo; the first two coordinates are constant phases.

```python
def analytic_state(x, data):
    u = [data["u_offset"][0], data["u_offset"][1], data["u_offset"][2] + x]
    wp2, wp4, wp6, wp2p, wp4p, wp6p, wp2pp, wp4pp, wp6pp = kleinian_p(
        u, data["omega"], data["tau"], data["kappa"],
        ((2, 2), (1, 2), (0, 2), (2, 2, 2), (1, 2, 2), (0, 2, 2),
         (2, 2, 2, 2), (1, 2, 2, 2), (0, 2, 2, 2)),
        data["characteristic"])
    return (-wp2, -wp4, -wp6,                      # u1, u2, u3
            wp2p / 2, wp4p / 2, wp6p / 2,          # v1, v2, v3
            wp2,                                   # w1
            wp2pp / 2 - wp4 - 2 * wp2 ** 2,        # w2
            wp4pp / 2 - wp6 - 2 * wp2 * wp4,       # w3
            wp6pp / 2 - 2 * wp2 * wp6)             # w4
```

One `kleinian_p` call evaluates all nine functions, of orders two to four, from a single genus-three theta jet — about 0.2 s per analytic state at 30 decimal places, so a three-phase run takes a few minutes.

## 4. Choosing the initial point

The genus-two story repeats: a generic real offset gives real but unbounded solutions, and the imaginary half-period translate $\omega\tau[1, 1, 1]^T$ moves the trajectory onto the **compact real component** of the Jacobian, where all ten coordinates stay real and bounded (here $\max|\wp| \approx 54$). There are now two constant phases, $z_5$ and $z_3$ (coordinates 0 and 1); the demo varies the first and fixes the second.

```python
lattice_vector = mp.matrix([1, 1, 1])
u_offset = (mp.matrix([mp.zero, mp.mpf("0.27"), mp.mpf("0.41")])
            + omega * tau * lattice_vector)
```

## 5. Verification

The same four checks run as in the genus-two demo, with two genuinely new identity verifications:

1. **RK4 comparison.** The ten expanded ODEs are integrated with classical RK4 from the analytic initial state.
2. **Spectral polynomial.** $H_\xi = U_\xi W_\xi + V_\xi^2$ must be identically $F(\xi)/4$, now a degree-seven polynomial with eight coefficients ($h_1 = 0$ from the canonical form). It held to $\sim 10^{-28}$ before any RK4 was run.
3. **Kleinian identities.** The $i = 1$ case is unchanged, including the $\lambda_4/2$ constant of this sigma normalisation. The $i = 2$ and $i = 3$ cases are new at genus three — the former brings in $\wp_6$, the latter exercises the absent-$\wp_8$ boundary and $\wp_{3,5}$, neither of which exists at genus two:

$$ \wp''_2 = 6\wp_2^2 + 4\wp_4 + \tfrac{\lambda_4}{2}, \qquad \wp''_4 = 6(\wp_2\wp_4 + \wp_6) - 2\wp_{3,3}, \qquad \wp''_6 = 6\wp_2\wp_6 - 2\wp_{3,5}. $$

4. **Abel-map closure.** The roots $x_i$ of $U(x) = x^3 + u_1x^2 + u_2x + u_3$ and the values $y_i = 2V(x_i) = 2(v_1x_i^2 + v_2x_i + v_3)$ form the degree-three divisor. Passing it to `curve.abel_map_kind_1(..., reduce=True)` recovers the chosen $z(0)$ modulo the full period lattice, and evaluating the solution at the recovered representative reproduces all ten state coordinates — the first genus-three end-to-end exercise of the Abel map:

```python
roots = mp.polyroots([u3, u2, u1, mp.one])               # roots of U_xi
divisor = [(x, 2 * (v1 * x**2 + v2 * x + v3)) for x in roots]
image = curve.abel_map_kind_1(divisor, reduce=True).value
```

At 30 decimal places with 3200 RK4 steps over $t \in [-4, 4]$, three phases:

```
phase     max u1 error   max u2 error   max u3 error   max w4 error   spectral drift
  -6.0     2.34e-8        4.59e-8        4.33e-8        3.36e-7        1.49e-7
  0.13     3.20e-8        3.84e-8        5.98e-8        3.49e-7        1.03e-7
   6.0     2.81e-8        5.58e-8        4.67e-8        3.63e-7        1.61e-7
Kleinian wp2'' identity residual: 6.74e-30
Kleinian wp4'' identity residual: 2.52e-29
Kleinian wp6'' identity residual: 3.17e-29
Abel-map divisor curve residual: 8.23e-27
Abel-map period-lattice residual: 1.51e-31
Abel-map recovered-state residual: 2.21e-29
```

The coordinate errors sit at the RK4 truncation scale ($w_4$ is largest, as it carries the highest-order jet), and the identity and Abel residuals at rounding level.

## 6. Output


The full runnable script is `neumann_moser_genus_three_demo.py`; run it from the repository root with

```
.venv/bin/python -m examples.neumann_moser.neumann_moser_genus_three_demo
```
