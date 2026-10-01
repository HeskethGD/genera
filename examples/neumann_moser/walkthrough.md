# Solving the Neumann-Moser system with Genera Kleinian functions

This note walks through the genus-two demo in this folder (`neumann_moser_kleinian_demo.py`) end to end: the equations, the solution formulas, the assumptions and fixes made along the way, and how the curve and initial point were chosen. The reference is P. G. Baron, [arXiv:2402.18079](https://arxiv.org/abs/2402.18079), with the solution theory quoted from V. M. Buchstaber, [arXiv:2402.09218](https://arxiv.org/abs/2402.09218). Corrections and convention details are documented at greater length in `project_update_part_18_neumann_moser_demo_2026-09-18.md`.

## 1. The system

The Neumann-Moser system of Baron's paper lives on $\mathbb{C}^{3n+1}$. In the smallest case $n = g = 2$ it has seven coordinates $(u_1, u_2, v_1, v_2, w_1, w_2, w_3)$, and with $\Gamma = w_1 - u_1$ its expanded equations are

$$
\begin{aligned}
\dot u_1 &= -2v_1, & \dot v_1 &= -\Gamma u_1 - u_2 + w_2, & \dot w_1 &= 2v_1, \\
\dot u_2 &= -2v_2, & \dot v_2 &= -\Gamma u_2 + w_3, & \dot w_2 &= 2v_2 + 2\Gamma v_1, \\
& & & & \dot w_3 &= 2\Gamma v_2.
\end{aligned}
$$

Two of the paper's printed equations needed repair here: the middle $\dot w$ range must run to $n$ (giving $\dot w_2 = 2v_2 + 2\Gamma v_1$ for $n = 2$), and the last equation is $\dot w_{n+1} = 2\Gamma v_n$, not $2\Gamma u_n$. Both follow from the generating form $\dot W_\xi = 2(\xi + \Gamma)V_\xi$ and both are confirmed by the RK4 comparison below.

The system is the Moser image of the classical C. Neumann problem (in imaginary time), coincides with the Mumford system on its $t_1$ flow, and is the 3-stationary KdV hierarchy; $\Gamma = 2\wp_2$ is the KdV solution.

## 2. The curve and the Genera data

Section 7.2 of the paper uses the canonical odd-degree curve

$$ y^2 = F(x) = 4x^5 + \lambda_4 x^3 + \lambda_6 x^2 + \lambda_8 x + \lambda_{10}, $$

with no $x^4$ term, so the five branch points must sum to zero. The demo picks five distinct real roots summing to zero, $(-4, -1.5, -0.3, 1.3, 4.5)$, giving $\lambda_4 = -80.56$ and $\lambda_6 = -34.56$. Real roots keep $\tau$ purely imaginary and admit real trajectories.

```python
from genera import kleinian_p
from mpmath import mp
mp.dps = 30

roots = (mp.mpf(-4), mp.mpf("-1.5"), mp.mpf("-0.3"),
         mp.mpf("1.3"), mp.mpf("4.5"))
coefficients = [mp.mpf(4)]          # 4 * prod(x - root)
for root in roots:
    coefficients = multiply_by_linear(coefficients, root)

curve = mp.algebraic_curve(coefficients)
first = curve.periods_kind_1()
second = curve.periods_kind_2()
omega, tau = first.omega, first.tau
kappa = second.kappa
characteristic = curve.riemann_constant().characteristic
```

The curve records return the first-kind half-period matrix, the normalised period matrix, the second-kind combination $\kappa$ and the Riemann characteristic, ready for the Kleinian function evaluations.

## 3. The Kleinian solution

The solution flows along the $z_1$ direction of the Jacobian, with $z(t) = z(0) + t\,e_1$ and $' = d/dz_1$:

$$
u_i = -\wp_{2i}, \qquad v_i = \tfrac12 \wp'_{2i}, \qquad W_\xi = (\xi + 2\wp_2)\,p_I + \tfrac12 p_{III},
$$

where $p_I = \xi^2 - \wp_2\xi - \wp_4$ and $p_{III} = \wp''_2\xi + \wp''_4$. Expanding $W_\xi$ gives

$$
w_1 = \wp_2, \qquad w_2 = \tfrac12\wp''_2 - \wp_4 - 2\wp_2^2, \qquad w_3 = \tfrac12\wp''_4 - 2\wp_2\wp_4.
$$

The paper prints $p_{II}$ in place of $p_I$ in $W_\xi$; that cannot be right by degree count alone (it gives a degree-$g$ polynomial where $W_\xi$ has degree $g+1$), and Buchstaber's theorem confirms $p_I$. With this correction everything below checks numerically.

**Index convention.** The paper writes $\wp_{2k} = -\partial^2 \log\sigma/\partial z_1 \partial z_{2k-1}$. Its $z_1$ is Genera Abelian coordinate **1** and its $z_3$ is coordinate **0** — found empirically via the Kleinian identities, but in fact forced by weights: with $\mathrm{wt}(x) = 2$ and $\mathrm{wt}(y) = 2g+1$ the coordinate from $x^r dx/y$ has weight $2g-1-2r$, while the paper's KdV times are weight-labelled $z_{2k-1}$, so $z_{2k-1}$ is always Genera coordinate $g-k$. Genera's ascending-power ordering ($dx/y$ first) matches the CEEK and Enolskii references it is validated against, as well as Sage and Abelfunctions; the paper's weight ordering is simply the reverse. Thus:

| Paper | `kleinian_p` indices |
| --- | --- |
| $\wp_2$ | `(1, 1)` |
| $\wp_4$ | `(0, 1)` |
| $\wp_{3,3}$ | `(0, 0)` |
| $\wp'_2, \wp'_4$ | `(1,1,1)`, `(0,1,1)` |
| $\wp''_2, \wp''_4$ | `(1,1,1,1)`, `(0,1,1,1)` |

```python
def analytic_state(x, data):
    u = [data["u_offset"][0], data["u_offset"][1] + x]   # flow in coord 1
    wp2, wp4, wp2p, wp4p, wp2pp, wp4pp = kleinian_p(
        u, data["omega"], data["tau"], data["kappa"],
        ((1, 1), (0, 1), (1, 1, 1), (0, 1, 1),
         (1, 1, 1, 1), (0, 1, 1, 1)),
        data["characteristic"])
    return (-wp2,                                   # u1
            -wp4,                                   # u2
            wp2p / 2, wp4p / 2,                     # v1, v2
            wp2,                                    # w1
            wp2pp / 2 - wp4 - 2 * wp2 ** 2,        # w2
            wp4pp / 2 - 2 * wp2 * wp4)              # w3
```

One `kleinian_p` call evaluates all six functions (orders two to four) from a single theta jet, so the whole analytic state costs one evaluation per time point.

## 4. Choosing the initial point

The paper (like the Manakov reference) does not say how to reach $z(0)$ from physical initial data. The demo still chooses a convenient $z(0)$ directly, as follows, but `curve.abel_map_kind_1` now lets it close the inverse loop afterwards.

* A generic real offset produces real but *unbounded* solutions: the $z_1$-line crosses near the theta divisor and the $\wp$-functions reach $10^{6}$–$10^{10}$ on $[-4, 4]$.
* Translating by the imaginary half-period combination $\omega\tau[1, 1]^T$ moves the line onto the **compact real component** of the Jacobian, where all seven coordinates stay real and bounded (here $\max|\wp| \approx 60$). This was found by scanning the sixteen half-period translates and ranking the maximum magnitude of the trajectory.
* The remaining freedom is the constant $z_3$ coordinate (the demo's `--phases` values) and the small real offset; any generic choice works. The demo uses $z(0) = [0, 0.27] + \omega[1,1] + \omega\tau[1,1]$ plus the phase.

```python
lattice_vector = mp.matrix([1, 1])
u_offset = (mp.matrix([mp.zero, mp.mpf("0.27")])
            + omega * lattice_vector + omega * tau * lattice_vector)
```

## 5. Verification

Four independent checks run alongside the trajectory:

1. **RK4 comparison.** The seven expanded ODEs are integrated with classical RK4 from the analytic initial state; agreement at the level of pure truncation error validates the ODE transcription and the solution formulas together.
2. **Spectral polynomial.** $H_\xi = U_\xi W_\xi + V_\xi^2$ must be constant and equal to the fixed curve $F(\xi)/4$ (the paper's integrals $h_i$ are $\lambda_i/4$, including $h_1 = 0$ forced by the canonical curve form). This held to $10^{-29}$ before any RK4 was run, isolating the solution formulas from the ODEs.
3. **Kleinian identities.** The genus-two specialisations of the paper's relation (55) are checked as residuals, with the curve constant at $\lambda_4/2$ in this sigma normalisation — one quarter of the printed $2\lambda_4$, and exactly the classical Weierstrass constant $-g_2/2$: the genus-one reduction of the same Genera sigma reproduces $(\wp')^2 = 4\wp^3 - g_2\wp - g_3$ to rounding level, so the outlier is the printed constant, not this normalisation (the Manakov demo sees the same constant in Christiansen et al. equation (3.16)):

$$ \wp''_2 = 6\wp_2^2 + 4\wp_4 + \tfrac{\lambda_4}{2}, \qquad \wp''_4 = 6\wp_2\wp_4 - 2\wp_{3,3}. $$

4. **Abel-map closure.** From the analytic initial state, the roots $x_i$ of $U(x)=x^2+u_1x+u_2$ and the values $y_i=2V(x_i)=2(v_1x_i+v_2)$ form a degree-two divisor on the spectral curve. Passing this divisor to `curve.abel_map_kind_1(..., reduce=True)` recovers the chosen $z(0)$ modulo the full period lattice. Evaluating the solution at the recovered representative reproduces all seven state coordinates. This checks the factor of two in $y_i$, both sheet choices, the absence of a separate Riemann-constant shift and the compatibility of the Abel-map paths with the period basis.

At 30 decimal places with 3200 RK4 steps over $t \in [-4, 4]$:

```
phase     max u1 error   max u2 error   max w3 error   spectral drift
  -6.0     1.60e-8        1.59e-8       1.12e-7        1.45e-8
  0.13     1.78e-8        1.18e-8       8.83e-8        2.11e-8
   6.0     1.65e-8        1.67e-8       1.14e-7        1.10e-8
Kleinian wp2'' identity residual: 6.3e-30
Kleinian wp4'' identity residual: 6.3e-30
Abel-map divisor curve residual: 2.0e-28
Abel-map period-lattice residual: 5.1e-32
Abel-map recovered-state residual: 6.3e-30
```

The coordinate errors sit exactly at the RK4 truncation scale and the identity residuals at rounding level.

## 6. Output


The full runnable script is `neumann_moser_kleinian_demo.py`; run it from the repository root with

```
.venv/bin/python -m examples.neumann_moser.neumann_moser_kleinian_demo
```
