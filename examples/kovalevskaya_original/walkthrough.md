# The Kovalevskaya top from the original 1889 memoir

`kovalevskaya_original_demo.py` reconstructs the solution of the
Kovalevskaya top from Sophie Kowalevski's *Sur le
problème de la rotation d'un corps solide autour d'un point fixe*, Acta
Mathematica **12** (1889) 177--232 (English translation in `latex/`),
using Genera's algebraic-curve engine and Riemann theta functions, and
compares it against a high-precision RK4 integration of the equations
of motion.  The companion folder `../kovalevskaya_top/` treats the same
top through the modern Bobenko--Reyman--Semenov-Tian-Shansky Lax-pair
paper; this example follows the original memoir alone.

Nothing is taken from the memoir's numerical examples: the integrals of
motion, the initial data, the curve and the evaluation points are chosen
for this demo. The marking-dependent theta characteristics are derived
from branch-point half-periods, and the corrections to printed formulas
are supported by the memoir's own algebraic identities. RK4 is retained
as an independent, lower-accuracy dynamical comparison.

## What the memoir gives

Section 2 reduces the Kowalevski case ($A = B = 2C$, $z_0 = 0$,
normalised to $C = 1$, $y_0 = 0$, $c_0 = Mgx_0$) to the system for
$(p, q, r, \gamma, \gamma', \gamma'')$ -- the code calls the last three
$(g_0, g_1, g_2)$ --

$$
2\dot p = qr, \qquad
2\dot q = -pr - c_0 \gamma'', \qquad
\dot r = c_0 \gamma',
$$

$$
\dot\gamma = r\gamma' - q\gamma'', \qquad
\dot\gamma' = p\gamma'' - r\gamma, \qquad
\dot\gamma'' = q\gamma - p\gamma',
$$

with the four algebraic integrals

$$
2(p^2+q^2) + r^2 = 2c_0\gamma + 6l_1, \qquad
2(p\gamma + q\gamma') + r\gamma'' = 2l,
$$

$$
\gamma^2 + \gamma'^2 + \gamma''^2 = 1, \qquad
\bigl((p+qi)^2 + c_0(\gamma + i\gamma')\bigr)
\bigl((p-qi)^2 + c_0(\gamma - i\gamma')\bigr) = k^2 .
$$

Section 4 separates the motion: with $x_1 = p + qi$, $x_2 = p - qi$,

$$
R(x) = -x^4 + 6l_1 x^2 + 4lc_0 x + c_0^2 - k^2,
$$

$$
R(x_1 x_2) = -x_1^2 x_2^2 + 6l_1 x_1 x_2 + 2lc_0(x_1 + x_2) + c_0^2 - k^2
$$

(a two-variable function of the memoir, *not* the quartic at the
product), the change of variables

$$
s_{1,2} = \frac{R(x_1 x_2) \mp \sqrt{R(x_1)R(x_2)}}{2(x_1-x_2)^2}
        + \frac{l_1}{2}
$$

puts the flow on the genus-two curve

$$
y^2 = R_1(s) = -4(s-e_1)(s-e_2)(s-e_3)(s-k_1)(s-k_2),
$$

$$
k_{1,2} = \frac{l_1 \pm k}{2}, \qquad
g_2 = k^2 - c_0^2 + 3l_1^2, \qquad
g_3 = l_1(k^2 - c_0^2 - l_1^2) + l^2 c_0^2,
$$

where $e_1 > e_2 > e_3$ are the roots of $4s^3 - g_2 s - g_3 = 0$, with
the time-linearised Abel flow (eq. 15)

$$
dt = \frac{s_1\,ds_1}{\sqrt{R_1(s_1)}} + \frac{s_2\,ds_2}{\sqrt{R_1(s_2)}},
\qquad
0 = \frac{ds_1}{\sqrt{R_1(s_1)}} + \frac{ds_2}{\sqrt{R_1(s_2)}} .
$$

Section 7's regime $l_1 > k > c_0 > 0$, $l^2 < (3l_1-k)/2$ gives five
real roots ordered $k_1 > e_1 > e_2 > k_2 > e_3$ (the demo uses $l_1 = 2$,
$k = 3/2$, $c_0 = 1$, $l = 1/2$), $s_1$ oscillating in $(e_1, k_1)$ and
$s_2$ below $e_3$.  Sections 5--7 then give the physical variables
through fifteen P-quantities, $P_a = \sqrt{(s_1-a)(s_2-a)}$ and $P_{ab}$
(the memoir's general definition with unit $c$-constants), each a
Rosenhain theta quotient, combined as

$$
p = -i\,\frac{L_1 P_1 + M_1 P_2 + N_1 P_3}{\Delta}, \qquad
q = \frac{E}{\Delta}, \qquad \Delta = L P_1 + M P_2 + N P_3,
$$

$$
r = -i\,\frac{L P_{23} + M P_{13} + N P_{12}}{\Delta}, \qquad
c_0 \gamma'' = \frac{L_1 P_{23} + M_1 P_{13} + N_1 P_{12}}{\Delta},
$$

$$
c_0 \gamma' = -\frac{i}{2}\,\frac{\mathrm{bracket}}{\Delta^2}, \qquad
2c_0 \gamma = -4l_1 + \frac{L_2 P_1 + M_2 P_2 + N_2 P_3}{\Delta}
              - \Bigl(\frac{r}{i}\Bigr)^{\!2},
$$

with $E = (e_2-e_3)(e_3-e_1)(e_1-e_2)$ and $L, M, N$ built from sigma
multipliers $\sigma_\lambda(w)/\sigma(w)$ at the point $w$ with
$\wp(w) = -l_1$ (see "Conventions").

## What Genera provides

* `Curve({(0, 2): 1, **{(i, 0): -c for i, c in enumerate(coefficients)}})` -- the genus-2 hyperelliptic curve
  from the ascending coefficients of $y^2 = R_1(s)$: genus, branch
  locus (matches the five memoir roots to $8\cdot10^{-31}$).
* `curve.periods_kind_1()` -- $\omega$, $\omega'$ and the Riemann
  matrix $\tau$ of the automatic (Baker-marked) hyperelliptic engine.
  $\omega$ is the matrix of half-periods: the period lattice of the
  engine's Abel coordinates is $[2\omega,\, 2\omega']$.
* `curve.abel_map_kind_1(places)` -- the Abel map of a divisor given as
  $(x, y)$ pairs, in unnormalised engine coordinates.
* `curve.lattice_reduce(value, periods)` -- period-lattice reduction,
  used to repair the $O(1)$ jumps of fresh `abel_map_kind_1` values when the
  engine's integration path crosses a branch cut.
* `rtheta(z, tau, characteristic=(a, b))` -- the genus-two theta
  function with arbitrary real characteristics; the sixteen
  two-torsion characteristics are the demo's entire theta vocabulary.
* `mp.polyroots`, `mp.findroot` and 30-digit arithmetic for the
  initial-data construction.

## Pipeline

Stages (`--stage rk4|curve|abel|theta|demo|all`):

1. **rk4** -- constructs initial data with exactly the prescribed
   integrals (parametrising the Kowalevski integral by an angle,
   `findroot` for the area integral, scanning the direction cosines for
   the state whose separation variables sit farthest from the branch
   points), then validates RK4: invariant drift $2.5\cdot10^{-9}$ over
   $[0, 20]$.
2. **curve** -- builds the curve, checks genus 2, the branch locus and
   $\tau$ (symmetric, purely imaginary for this real curve).
3. **abel** -- follows the divisor $(s_1(t), s_2(t))$ along the RK4
   trajectory and verifies the memoir's central claim: the Abel map is
   a straight line, $u(t) = u_0 + Vt$. Equation (15) and the selected
   sheets give $V = (0, -1)$ in the engine basis, verified to
   $4.4\cdot10^{-11}$ over $[0, 0.6]$, after continuing the engine's
   period jumps.
4. **theta** -- resolves each branch-point Abel image as a half-period,
   adds its characteristic to the Riemann-constant characteristic, and
   thereby derives each representation
   $P = \kappa\,\theta[\chi](v)/\theta[\chi'](v)$ with
   $v = (2\omega)^{-1}u$. Each $\kappa$ is fixed at one divisor from the
   defining algebraic $P$-value and verified at the remaining divisors.
5. **demo** -- evaluates the memoir's final formulas purely from theta
   functions along $v(t) = v_0 + Vt$ on $[0, 2.4]$ and compares
   against RK4, with the invariants included.

## Conventions, and how they differ from the memoir

* **Theta series.**  The memoir's
  $\vartheta(v) = \sum_{\nu\in\mathbb{Z}^2}
  \exp\{\pi i\,\nu^T\tau\nu + 2\pi i\,\nu^T v\}$ (section 6) is
  *exactly* Genera's `rtheta` convention.  But the memoir's $\tau$
  belongs to its own Rosenhain (real-oval) marking, while Genera's
  $\tau$ belongs to the automatic Baker marking -- different matrices.
  So the memoir's $\vartheta_\lambda$ *labels* are not transplanted.
  Instead, the characteristics are derived afresh in Genera's marking
  from its branch-point half-periods. The memoir's section-7 constants
  $c_\lambda$ are replaced by equivalent quotient multipliers fixed from
  the algebraic normalization.
* **Abel coordinates.**  `abel_map_kind_1` returns unnormalised values with
  quasi-periods $(2\omega, 2\omega')$; the theta argument is
  $v = (2\omega)^{-1}u$ by definition. Replacing $u$ by $-u$ gives the
  equivalent representation obtained from theta parity; it is not a
  distinct numerical normalization. The Riemann constant occurs as the
  common denominator characteristic rather than as an added argument.
* **Characteristics.**  Genera takes literal vectors $(a, b)$; the demo
  prints them as `[a1a2;b1b2]` bit patterns.  All fifteen P's share one
  denominator characteristic, $[11;01]$ -- the analogue of the memoir's
  common $\vartheta_5$ -- and the fifteen numerators exhaust the
  remaining two-torsion points. If $\delta_i$ is the half-characteristic
  of branch point $a_i$ and $K=[11;01]$, the derived rules are
  $\chi(P_i)=K+\delta_i$ and
  $\chi(P_{ij})=K+\delta_i+\delta_j$ modulo one.
* **Sigma branches.**  The memoir's constants are built from
  $\sigma_\lambda(w)/\sigma(w) = \sqrt{-(l_1 + e_\lambda)}$; the
  demo's orientation is the branch
  $\sigma_\lambda(w)/\sigma(w) = -i\sqrt{l_1 + e_\lambda}$, fixed by
  continuation on the selected real oval and the initial sheet.
* **Flow velocity.**  The memoir's $(du_1, du_2) = (dt, 0)$ appears as
  $V = (0, -1)$ in Genera's engine coordinates -- the same statement in
  a different basis.
* **Reality bookkeeping.**  The memoir's section-5 reality analysis
  ($L, M, N$ and $P_1, P_2, P_3$ imaginary; $L_1, M_1, N_1$ real)
  reappears verbatim: the derived quotient multipliers come out real or
  purely imaginary in exactly those slots, and the assembled solution
  has imaginary part $5\cdot10^{-33}$.

## What is derived, normalized and checked

**Taken verbatim from the memoir and verified** (residuals quoted from
the run): the reduced system and four integrals (exact by
construction); the separation formulas, including the two-variable
$R(x_1 x_2)$ -- the identity
$R(x_1)R(x_2) - R(x_1 x_2)^2 = (x_1-x_2)^2 R_1(x_1 x_2)$ holds to
roundoff and pins this reading; the spectral curve and root ordering;
the linearisation (eq. 15), verified as the straight-line Abel flow;
the general $P_a$/$P_{ab}$ definition and the three P-identities, which
hold to $2.6\cdot10^{-29}$ once the radical orientation is fixed; and
the section-5 formulas for $p$, $q$, $r$ and $\gamma''$, which
reproduce RK4 exactly.

The marking-dependent data are obtained as follows:

1. The theta argument is $v=(2\omega)^{-1}u$ by the definition of the
   normalized Riemann matrix.
2. Each branch-point Abel image is resolved against
   $[2\omega,2\omega']$. Its half-integral lattice coordinates give its
   literal `rtheta` characteristic. Adding these to the Riemann-constant
   characteristic gives all fifteen numerator characteristics.
3. Each quotient multiplier $\kappa$ is fixed from the defining algebraic
   value of the corresponding P-function at the initial divisor. It is
   then held fixed; the residual over the other twelve divisors is at most
   $3.7\cdot10^{-29}$.
4. The radical orientation is the continuous orientation on the selected
   real oval for which $q(0)=E/\Delta$. This selects an initial sheet, not
   a fitted coefficient. The P-identities independently check it.
5. The branches $\sigma_\lambda(w)/\sigma(w)=
   -i\sqrt{l_1+e_\lambda}$ use positive real square roots and the same
   selected sheet. Substitution recovers all six initial variables to
   $3.8\cdot10^{-30}$.
6. The initial data are constructed for this demo, not taken from the
   memoir.
7. The demonstration window $[0,2.4]$ exceeds the algebraic verification
   window $[0,0.6]$. The latter is chosen turning-point-free ($s_2$ turns near
   $t \approx 0.73$, where the divisor crosses the branch point to the
   other sheet), but the theta identity continues through it. The analytic
   solution satisfies the four invariants near working precision and is
   also differentiated directly in the six Euler--Poisson equations.

## Potential errors in the printed transcription

The table separates what is printed from what is used or inferred. Some
entries may be errors in the 1889 printing and some may have entered during
modern transcription; the original facsimile must be consulted to decide
between those possibilities.

| Location | Printed | Corrected or intended | Internal evidence |
|---|---|---|---|
| Section 3, equation (14), final term of $h_2'''$ | $(e_1^2-e_3^2)\sigma_3$ | probably $(e_1^2-e_2^2)\sigma_3$ | Restores the cyclic pattern of $h_2,h_2',h_2''$; already noted as TN2 in the translation. It is not used by this demo. |
| Section 5, $L_1,M_1,N_1$ | Positive products of $\sqrt{l_1+e_j}$ | Each product has an additional minus sign | The defining sigma ratios are all $\pm i\sqrt{l_1+e_j}$, so their pairwise products contain $i^2=-1$. |
| Section 5, displayed $P_{23},P_{13},P_{12}$ | The special radical expressions printed immediately after the formula for $r$ | Use the preceding general definition of $P_{\alpha\beta}$ | The special displays have the wrong distribution of $s_1,s_2$ factors and contradict both the reality discussion and the later section-7 table. |
| Section 5, $\gamma''$ | Denominator $LP_{23}+MP_{13}+NP_{12}$ | Denominator $\Delta=LP_1+MP_2+NP_3$ | The printed quotient is imaginary in the stated real regime. Differentiating $q=E/\Delta$ gives the corrected denominator. |
| Section 5, final $\gamma'$ formula | Prefactor $+i/2$ | Prefactor $-i/2$ | It follows directly by differentiating $r=-iX/\Delta$. |
| Same $\gamma'$ bracket | The three $(e_j-e_k)^2P_l$ terms have coefficient 1 | Each has coefficient $R_0=-4$ | The immediately preceding fourth P-identity contains precisely these $R_0$ factors. |
| Section 5, final $\gamma$ formula | $(L_2P_1+M_2P_2+N_2P_3)/\Delta$ | $R_0(L_2P_1+M_2P_2+N_2P_3)/\Delta$ | Substitution into the energy integral requires the factor. Keeping it in the final formula preserves the printed sigma-product definitions of $L_2,M_2,N_2$. |
| Section 7, ordered roots | $k_1=a_0$ and $e_3=a_0$ | $a_0=e_3,\ a_1=k_2,\ a_2=e_2,\ a_3=e_1,\ a_4=k_1$ | Forced by the preceding strict ordering and confirmed by the following P-table. |
| Section 7, second radical of $P_{35}$ | The factor $(s_2-e_1)$ occurs twice | One occurrence should be $(s_2-e_2)$ | Follows from the complementary-root pattern and the general $P_{\alpha\beta}$ definition; already noted as TN3. |

The translation also notes the anomalous equation number “(14)” after “(9)”.
That is editorial numbering rather than a mathematical error.

## Results (dps 30)

    initial data               integrals exact by construction
    RK4 drift [0, 20]          2.5e-9
    branch locus vs roots      7.9e-31
    Abel flow linearity        4.4e-11   (V = (0, -1))
    P-identities               2.6e-29 over the verification window
    derived theta quotients     max held-out residual 3.7e-29
    initial-condition check     3.8e-30
    theta vs RK4 on [0, 2.4]   p 1.2e-10 ... g1 1.3e-9 (RK4 truncation)
    imaginary part             6.9e-33
    invariants, theta solution 2.1e-29
    direct ODE residual         1.3e-29

The theta side is limited only by dps: every theta-identity residual away
from the normalization point sits at the $10^{-29}$ level, and the $10^{-10}$--$10^{-9}$
comparison errors are RK4's own fourth-order truncation, removable with
more substeps.

## Running

From the Genera repository root:

    .venv/bin/python -m examples.kovalevskaya_original.kovalevskaya_original_demo
    .venv/bin/python -m examples.kovalevskaya_original.kovalevskaya_original_demo --stage theta

Stages are `rk4`, `curve`, `abel`, `theta`, `demo`, `all`; every stage
prints its residuals. The `demo` stage exits unsuccessfully if the maximum
theta-versus-RK4 component error exceeds `--tol`.
