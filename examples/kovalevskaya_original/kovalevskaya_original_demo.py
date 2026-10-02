#!/usr/bin/env python3
"""Kovalevskaya top: the original 1889 memoir's solution versus RK4.

This demo follows the English translation of S. Kowalevski, "Sur le
probleme de la rotation d'un corps solide autour d'un point fixe",
Acta Mathematica 12 (1889) 177-232, kept at latex/kovalevskaya_top.tex,
exercising Genera's algebraic-curve and Riemann-theta machinery.

Section 2 of the memoir reduces the Kowalevski case (A = B = 2C, z0 = 0,
normalised to C = 1, y0 = 0, c0 = M g x0) to

    2 dp/dt = q r,            2 dq/dt = -p r - c0 g2,
    dr/dt   = c0 g1,
    dg0/dt  = r g1 - q g2,    dg1/dt = p g2 - r g0,   dg2/dt = q g0 - p g1,

with the four algebraic integrals

    2(p^2+q^2) + r^2      = 2 c0 g0 + 6 l1,
    2(p g0 + q g1) + r g2 = 2 l,
    g0^2 + g1^2 + g2^2    = 1,
    ((p+qi)^2 + c0(g0 + i g1))((p-qi)^2 + c0(g0 - i g1)) = k^2.

Putting x1 = p + qi, x2 = p - qi and R(x) = -x^4 + 6 l1 x^2 + 4 l c0 x
+ c0^2 - k^2, the memoir's change of variables (section 4)

    s1 = (R(x1 x2) - sqrt(R(x1) R(x2))) / (2 (x1-x2)^2) + l1/2,
    s2 = (R(x1 x2) + sqrt(R(x1) R(x2))) / (2 (x1-x2)^2) + l1/2

maps the flow onto the genus-two curve (section 4 eq. 15 / section 7)

    y^2 = R1(s) = -4 (s - e1)(s - e2)(s - e3)(s - k1)(s - k2),

    k1 = (l1 + k)/2,   k2 = (l1 - k)/2,
    4 s^3 - g2 s - g3 = 0 with roots e1 > e2 > e3,
    g2 = k^2 - c0^2 + 3 l1^2,   g3 = l1(k^2 - c0^2 - l1^2) + l^2 c0^2,

with the time-linearised Abel flow (section 4 eq. 15)

    dt = s1 ds1/sqrt(R1(s1)) + s2 ds2/sqrt(R1(s2)),
     0 =   ds1/sqrt(R1(s1)) + ds2/sqrt(R1(s2)).

Section 7's case l1 > k > c0 > 0, l^2 < (3 l1 - k)/2 gives five real
roots ordered k1 > e1 > e2 > k2 > e3, with s1 moving in (e1, k1) and
s2 below e3; the solution is then given by Rosenhain theta quotients.

Implementation notes
--------------------
* All marking-dependent quantities (tau, Abel coordinates, theta
  characteristics) come from the Genera automatic hyperelliptic engine.
  Branch-point Abel images determine their half-characteristics in the
  Baker marking.  Riemann's vanishing theorem then gives the characteristic
  of each of the fifteen P-quotients; no characteristic search is used.
* The P_ab follow the memoir's general definition ("if one sets
  R(s) = R0(s-a0)..."); the printed special forms of P23, P13, P12 in
  section 5 carry transcription anomalies, while the general definition
  gives the real values the memoir's own reality analysis requires.
* Sign conventions of the radicals are fixed by the memoir's
  P-identities (the three relations "one easily finds") and the selected
  real initial sheet. The corrected section-5 formulas reconstruct the
  initial state algebraically and the final theta solution is differentiated
  directly in the Euler--Poisson equations.

Stages (--stage): rk4, curve, abel, theta, demo, all.
"""

import argparse
from dataclasses import dataclass

from genera import algebraic_curve, rtheta
from mpmath import mp

from examples._rk4 import rk4_trajectory as integrate_rk4

mp.dps = 30

# Demo constants: the integrals (l1, k, l) and the parameter c0, chosen
# in section 7's regime l1 > k > c0 > 0 and l^2 < (3 l1 - k)/2.
L1 = mp.mpf(2)                 # l1: energy-type integral / 6
K_INTEGRAL = mp.mpf(3) / 2      # k:  sqrt of the Kowalevski integral
C0 = mp.mpf(1)                 # c0 = M g x0
L = mp.mpf(1) / 2               # l:  area-type integral / 2
R0 = mp.mpf(-4)                # leading factor of R1


# ---------------------------------------------------------------------------
# Euler-Poisson system of section 2 and RK4
# ---------------------------------------------------------------------------

def euler_poisson_rhs(state):
    """Right-hand side of the section 2 system for (p, q, r, g0, g1, g2)."""
    p, q, r, g0, g1, g2 = state
    return (
        q * r / 2,
        (-p * r - C0 * g2) / 2,
        C0 * g1,
        r * g1 - q * g2,
        p * g2 - r * g0,
        q * g0 - p * g1,
    )


def rk4_trajectory(state0, t_values, substeps=16):
    """Integrate from state0, returning states at each t in t_values."""
    return integrate_rk4(
        euler_poisson_rhs, state0, t_values, substeps=substeps)


def invariants(state):
    """The four algebraic integrals: (6 l1, 2 l, |g|^2-1, k^2)."""
    p, q, r, g0, g1, g2 = state
    xi1 = (p + 1j * q) ** 2 + C0 * (g0 + 1j * g1)
    xi2 = (p - 1j * q) ** 2 + C0 * (g0 - 1j * g1)
    return (
        2 * (p * p + q * q) + r * r - 2 * C0 * g0,
        2 * (p * g0 + q * g1) + r * g2,
        g0 * g0 + g1 * g1 + g2 * g2 - 1,
        xi1 * xi2,
    )



def initial_state():
    """Construct well-conditioned initial data with the given integrals.

    The direction cosines (g0, g1, g2) and the sign of r are scanned;
    the modulus structure of the Kowalevski integral,

        u = p^2 - q^2 + c0 g0,   v = 2 p q + c0 g1,   u^2 + v^2 = k^2,

    parametrises (p, q) by an angle, and the area integral fixes that
    angle.  Among all candidates the one whose separation variables
    sit farthest from the branch points is kept, so that the demo
    window contains no turning points of s1 or s2.
    """
    def build(theta, rsign, g0, g1, g2):
        uu = K_INTEGRAL * mp.cos(theta) - C0 * g0     # p^2 - q^2
        vv = K_INTEGRAL * mp.sin(theta) - C0 * g1     # 2 p q
        ssum = mp.sqrt(uu ** 2 + vv ** 2)             # p^2 + q^2
        p = mp.sqrt((ssum + uu) / 2)
        q = vv / (2 * p)
        r2 = 2 * C0 * g0 + 6 * L1 - 2 * ssum
        if r2 <= 0:
            return None
        return (p, q, rsign * mp.sqrt(r2), g0, g1, g2)

    def residual(theta, rsign, g0, g1, g2):
        state = build(theta, rsign, g0, g1, g2)
        if state is None:
            return mp.mpf(10) ** 6
        p, q, r = state[0], state[1], state[2]
        return 2 * (p * g0 + q * g1) + r * g2 - 2 * L

    def conditioning(state):
        """Distance of the separation variables from the branch points."""
        try:
            s1, s2 = separation(state)
        except RuntimeError:
            return None
        if K2 < s1 < E2:
            margin1 = min(s1 - K2, E2 - s1)
        else:
            margin1 = min(s1 - E1, K1 - s1)
        return min(margin1, E3 - s2)

    best = None
    seen = []
    for g0 in (mp.mpf("0.3"), mp.mpf("0.45"), mp.mpf("0.6"), mp.mpf("0.75")):
        for g2 in (mp.mpf("0.2"), mp.mpf("0.4"), mp.mpf("0.6")):
            square = 1 - g0 ** 2 - g2 ** 2
            if square <= 0:
                continue
            g1 = mp.sqrt(square)
            for rsign in (1, -1):
                for step in range(1, 16):
                    seed = step * mp.pi / 8
                    try:
                        root = mp.findroot(
                            lambda th: residual(th, rsign, g0, g1, g2),
                            seed)
                    except Exception:
                        continue
                    state = build(root, rsign, g0, g1, g2)
                    if state is None:
                        continue
                    if abs(residual(root, rsign, g0, g1, g2)) > mp.mpf(10) ** -20:
                        continue
                    if any(abs(root - other) < mp.mpf(10) ** -12
                           for other in seen):
                        continue
                    seen.append(root)
                    score = conditioning(state)
                    if score is not None and (best is None or score > best[0]):
                        best = (score, state)
    if best is None:
        raise RuntimeError("could not construct well-conditioned initial data")
    return best[1]


# ---------------------------------------------------------------------------
# Separation variables and the spectral curve
# ---------------------------------------------------------------------------

def quartic(x):
    """R(x) = -x^4 + 6 l1 x^2 + 4 l c0 x + c0^2 - k^2 (section 4)."""
    return (-x ** 4 + 6 * L1 * x ** 2 + 4 * L * C0 * x
            + C0 ** 2 - K_INTEGRAL ** 2)


def r_cross(x1, x2):
    """The memoir's two-variable function R(x1 x2) (section 4, eq. 8).

        R(x1 x2) = -x1^2 x2^2 + 6 l1 x1 x2 + 2 l c0 (x1 + x2) + c0^2 - k^2,

    which satisfies the identity R(x1) R(x2) - R(x1 x2)^2
    = (x1 - x2)^2 R1(x1 x2) with R1(x1 x2) = A C - B^2.
    """
    return (-(x1 * x2) ** 2 + 6 * L1 * (x1 * x2) + 2 * L * C0 * (x1 + x2)
            + C0 ** 2 - K_INTEGRAL ** 2)


def real_root(value):
    """Collapse a conjugate-pair artefact to its real part."""
    value = mp.mpc(value)
    if abs(mp.im(value)) < mp.mpf(10) ** -25 * (1 + abs(mp.re(value))):
        return mp.re(value)
    return value


def spectral_data():
    """Roots e1 > e2 > e3 and k1 > k2 in section 7's regime."""
    g2 = K_INTEGRAL ** 2 - C0 ** 2 + 3 * L1 ** 2
    g3 = L1 * (K_INTEGRAL ** 2 - C0 ** 2 - L1 ** 2) + L ** 2 * C0 ** 2
    cubic = sorted((real_root(r) for r in
                    mp.polyroots([-g3, -g2, 0, 4], asc=True)),
                   reverse=True)
    e1, e2, e3 = (mp.mpf(r) for r in cubic)
    k1 = (L1 + K_INTEGRAL) / 2
    k2 = (L1 - K_INTEGRAL) / 2
    if not (k1 > e1 > e2 > k2 > e3):
        raise RuntimeError("section 7 ordering k1>e1>e2>k2>e3 violated: %s"
                           % [mp.nstr(x, 6) for x in (k1, e1, e2, k2, e3)])
    return {"g2": g2, "g3": g3, "e1": e1, "e2": e2, "e3": e3,
            "k1": k1, "k2": k2}


SPECTRAL = spectral_data()
E1, E2, E3 = SPECTRAL["e1"], SPECTRAL["e2"], SPECTRAL["e3"]
K1, K2 = SPECTRAL["k1"], SPECTRAL["k2"]
FIVE = (E1, E2, E3, K1, K2)     # memoir labels a1..a5 = e1, e2, e3, k1, k2


def r1(s):
    """R1(s) = -4 (s-e1)(s-e2)(s-e3)(s-k1)(s-k2)."""
    return R0 * mp.fprod((s - a) for a in FIVE)


def curve_coefficients():
    """Ascending coefficients of y^2 = R1(s)."""
    poly = [mp.one]                       # descending, monic
    for a in FIVE:
        new = [mp.zero] * (len(poly) + 1)
        for i, c in enumerate(poly):
            new[i] += c
            new[i + 1] -= a * c
        poly = new
    return [R0 * c for c in reversed(poly)]


def separation(state):
    """Return (s1, s2): s1 in a bounded real oval, s2 below e3 (section 4).

    R1 is positive on (-inf, e3), (k2, e2) and (e1, k1); the trajectory
    of the demo data has s1 on one of the two bounded ovals and s2 on
    the unbounded component.
    """
    p, q = state[0], state[1]
    x1 = p + 1j * q
    x2 = p - 1j * q
    cross = r_cross(x1, x2)
    radical = mp.sqrt(quartic(x1) * quartic(x2))
    scale = 2 * (x1 - x2) ** 2
    pair = [real_root((cross - radical) / scale + L1 / 2),
            real_root((cross + radical) / scale + L1 / 2)]
    for candidate, other in ((pair[0], pair[1]), (pair[1], pair[0])):
        for low, high in ((K2, E2), (E1, K1)):
            if low < candidate < high and other < E3:
                return candidate, other
    raise RuntimeError("separation outside the real-motion ovals: s = %s"
                       % [mp.nstr(v, 8) for v in pair])


INITIAL = initial_state()



# ---------------------------------------------------------------------------
# Stage 1: RK4 validation of the section 2 system
# ---------------------------------------------------------------------------

def stage_rk4(verbose=True):
    """Validate the RK4 integrator and the four algebraic integrals."""
    mp.dps = 30
    state0 = INITIAL
    if verbose:
        print("initial state (p, q, r, g0, g1, g2):")
        print("  ", [mp.nstr(x, 12) for x in state0])
        target = (6 * L1, 2 * L, mp.zero, K_INTEGRAL ** 2)
        values = invariants(state0)
        print("integrals (6 l1, 2 l, |g|^2-1, k^2):")
        print("  target:", [mp.nstr(x, 12) for x in target])
        print("  actual:", [mp.nstr(x, 12) for x in values])
    s1, s2 = separation(state0)
    if verbose:
        print("separation: s1 =", mp.nstr(s1, 12), "  s2 =", mp.nstr(s2, 12))
        print("  (e1, k1) = (%s, %s), e3 = %s"
              % (mp.nstr(E1, 6), mp.nstr(K1, 6), mp.nstr(E3, 6)))
    t_values = [mp.mpf(i) / 20 for i in range(401)]
    states = rk4_trajectory(state0, t_values)
    drifts = [max(abs(a - b) for a, b in zip(invariants(s),
                                             invariants(state0)))
              for s in states]
    if verbose:
        print("max invariant drift over [0, 20]:", mp.nstr(max(drifts), 4))
    return state0


# ---------------------------------------------------------------------------
# Stage 2: the spectral curve, its periods and Riemann matrix
# ---------------------------------------------------------------------------

def stage_curve(verbose=True):
    """Build the genus-two curve and report its period data."""
    mp.dps = 30
    curve = algebraic_curve(curve_coefficients())
    first = curve.periods_kind_1()
    if verbose:
        print("curve: y^2 = -4 (s-e1)(s-e2)(s-e3)(s-k1)(s-k2)")
        print("  roots (e3, k2, e2, e1, k1):",
              [mp.nstr(a, 10) for a in sorted(FIVE)])
        print("  genus:", curve.genus)
        locus = sorted(curve.branch_locus.branch_values)
        print("  branch locus vs roots: max |diff| =",
              mp.nstr(max(abs(a - b) for a, b in
                          zip(locus, sorted(FIVE))), 4))
        print("  tau:")
        tau = first.tau
        for i in range(2):
            print("   ", [mp.nstr(tau[i, j], 10) for j in range(2)])
        tau = first.tau
        asym = max(abs(tau[i, j] - tau[j, i])
                   for i in range(2) for j in range(2))
        print("  tau symmetry residual:", mp.nstr(asym, 4))
        print("  omega:")
        for i in range(2):
            print("   ", [mp.nstr(first.omega[i, j], 10) for j in range(2)])
    return curve, first


# ---------------------------------------------------------------------------
# Stage 3: the divisor Abel map is a straight-line flow (eq. 15)
# ---------------------------------------------------------------------------

def divisor_places(state):
    """The two curve points above the separation variables."""
    s1, s2 = separation(state)
    return [(s1, mp.sqrt(r1(s1))), (s2, mp.sqrt(r1(s2)))]


def continuous_abel(curve, periods, states, verbose=False):
    """Abel coordinates continued across the engine's period jumps.

    Fresh abel_map_kind_1 values jump by period vectors whenever the engine's
    integration path crosses a branch cut; the true Jacobian flow is
    continuous, so increments larger than a half-period are reduced
    modulo the period lattice before being accumulated.
    """
    raw = [curve.abel_map_kind_1(divisor_places(state)) for state in states]
    values = [raw[0]]
    for i in range(1, len(raw)):
        step = raw[i] - raw[i - 1]
        reduced = curve.lattice_reduce(step, periods)
        step = reduced.value
        if verbose and any(reduced.shift):
            print("  period jump reduced at index", i)
        values.append(values[i - 1] + step)
    return values


def stage_abel(curve, periods, verbose=True):
    """Verify u(t) = u0 + V t for the divisor (s1(t), s2(t)) (eq. 15)."""
    mp.dps = 30
    # window chosen free of turning points of s1 and s2
    t_values = [mp.mpf(i) / 20 for i in range(13)]      # 0 .. 0.6
    states = rk4_trajectory(INITIAL, t_values)
    abel = continuous_abel(curve, periods, states, verbose=verbose)
    measured_velocity = (
        (abel[2] - abel[0]) / (t_values[2] - t_values[0]))
    # The engine orders (ds/y, s*ds/y).  Equation (15), together with the
    # selected sheets, therefore gives this velocity exactly; the measured
    # value below is only an independent check using the RK4 trajectory.
    velocity = mp.matrix([mp.zero, -mp.one])
    if verbose:
        print("Abel map of the divisor at t = 0:")
        print("   ", [mp.nstr(v, 12) for v in abel[0]])
        print("measured flow velocity (u(0.1) - u(0)) / 0.1:")
        print("   ", [mp.nstr(v, 12) for v in measured_velocity])
        print("velocity from equation (15) and the selected sheets:")
        print("   ", [mp.nstr(v, 12) for v in velocity])
    residuals = [mp.norm(abel[i] - abel[0] - velocity * t_values[i])
                 for i in range(len(t_values))]
    pairwise = [(abel[i + 1] - abel[i]) / (t_values[i + 1] - t_values[i])
                for i in range(len(t_values) - 1)]
    spread = max(mp.norm(step - velocity) for step in pairwise)
    if verbose:
        print("max |u(t) - u(0) - V t| over [0, 0.6]:",
              mp.nstr(max(residuals), 4))
        print("max pairwise velocity spread:", mp.nstr(spread, 4))
    return {"abel0": abel[0], "velocity": velocity, "abel": abel,
            "t_values": t_values, "states": states}


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--stage", default="all",
                        choices=("rk4", "curve", "abel", "theta", "demo",
                                 "all"))
    parser.add_argument("--tol", default="1e-6")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    if arguments.stage in ("rk4", "all"):
        stage_rk4()
    curve = periods = None
    if arguments.stage in ("curve", "abel", "theta", "demo", "all"):
        curve, periods = stage_curve()
    if arguments.stage in ("abel", "theta", "demo", "all"):
        abel_data = stage_abel(curve, periods)
    if arguments.stage in ("theta", "demo", "all"):
        config = stage_theta(curve, periods, abel_data)
    if arguments.stage in ("demo", "all"):
        result = stage_demo(
            curve, periods, abel_data, config)
        if max(result.errors) > mp.mpf(arguments.tol):
            raise SystemExit("theta solution failed its RK4 tolerance")
    return curve, periods


# ---------------------------------------------------------------------------
# Stage 4: derive the fifteen theta characteristics in the Baker marking
# ---------------------------------------------------------------------------

SLOTS = ("P1", "P2", "P3", "P4", "P5",
         "P12", "P13", "P14", "P15",
         "P23", "P24", "P25", "P34", "P35", "P45")
AB_SLOTS = ("P12", "P13", "P14", "P15",
            "P23", "P24", "P25", "P34", "P35", "P45")

HALVES = ((0, 0), (0, 1), (1, 0), (1, 1))


def characteristics():
    """The sixteen two-torsion characteristics in Genera format."""
    return [((mp.mpf(ax) / 2, mp.mpf(ay) / 2),
             (mp.mpf(bx) / 2, mp.mpf(by) / 2))
            for ax, ay in HALVES for bx, by in HALVES]


def p_values(state, cache=None, signs=None):
    """The fifteen P-quantities of sections 5-7 at the given state.

    P_a = sqrt((s1 - a)(s2 - a));  P_ab follows the memoir's general
    definition with unit c-constants,

        P_ab = P_a P_b / (s1 - s2) * [ W1/((s1-a_i)(s1-a_j))
                                       - W2/((s2-a_i)(s2-a_j)) ],

    where W1, W2 are the square roots of R1 at s1 and s2.  All radical
    branches are continued from ``cache`` (raw radical values); the
    7-bit ``signs`` dict (a0..a4, W1, W2) fixes the orientation, which
    the memoir's P-identities determine.
    """
    if cache is None:
        cache = {}
    s1, s2 = separation(state)
    values = {}

    def tracked(radical, key):
        value = mp.sqrt(radical)
        previous = cache.get(key)
        if previous is not None and abs(-value - previous) < abs(value - previous):
            value = -value
        cache[key] = value
        return value

    raw = {}
    for i in range(5):
        raw["a%d" % i] = tracked((s1 - FIVE[i]) * (s2 - FIVE[i]), "a%d" % i)
    raw["W1"] = tracked(r1(s1), "W1")
    raw["W2"] = tracked(r1(s2), "W2")

    def signed(key):
        if signs is not None and signs.get(key, 1) < 0:
            return -raw[key]
        return raw[key]

    for i in range(5):
        values["P%d" % (i + 1)] = signed("a%d" % i)
    w1 = signed("W1")
    w2 = signed("W2")
    for i in range(5):
        for j in range(i + 1, 5):
            values["P%d%d" % (i + 1, j + 1)] = (
                values["P%d" % (i + 1)] * values["P%d" % (j + 1)] / (s1 - s2)
                * (w1 / ((s1 - FIVE[i]) * (s1 - FIVE[j]))
                   - w2 / ((s2 - FIVE[i]) * (s2 - FIVE[j]))))
    return values


def identity_residuals(values):
    """The memoir's three P-identities ("one easily finds")."""
    def P(a):
        return values["P%d" % a]

    def Q(a, b):
        return values["P%d%d" % (min(a, b), max(a, b))]

    residuals = []
    for alpha, beta, gamma in ((1, 2, 3), (2, 3, 1), (3, 1, 2)):
        lhs = (R0 * P(alpha) * P(beta) * P(gamma)
               - Q(beta, gamma) * (P(gamma) * Q(alpha, gamma)
                                   - P(beta) * Q(alpha, beta))
               / (FIVE[beta - 1] - FIVE[gamma - 1]))
        residuals.append(abs(lhs - Q(alpha, 4) * Q(alpha, 5)))
        lhs = (R0 * P(alpha) * P(beta) * P(gamma)
               - Q(alpha, gamma) * (P(alpha) * Q(alpha, beta)
                                     - P(gamma) * Q(beta, gamma))
               / (FIVE[gamma - 1] - FIVE[alpha - 1]))
        residuals.append(abs(lhs - Q(beta, 4) * Q(beta, 5)))
        lhs = (R0 * P(alpha) * P(beta) * P(gamma)
               - Q(alpha, beta) * (P(beta) * Q(beta, gamma)
                                     - P(alpha) * Q(alpha, gamma))
               / (FIVE[alpha - 1] - FIVE[beta - 1]))
        residuals.append(abs(lhs - Q(gamma, 4) * Q(gamma, 5)))
    return residuals


def characteristic_sum(*items):
    """Add two-torsion characteristics componentwise modulo one."""
    return tuple(tuple(
        mp.mpf(sum(int(mp.nint(2 * item[row][column]))
                   for item in items) % 2) / 2
        for column in range(2)) for row in range(2))


def branch_characteristic(curve, periods, branch):
    """Derive a branch point's theta characteristic from its Abel image.

    If ``A = 2*omega*m + 2*omega_prime*n``, a branch point has half-integral
    ``m,n``.  In the literal rtheta convention its characteristic is
    ``(n,m)`` because theta[a,b](z) is translated by ``tau*a+b``.
    """
    basis = mp.matrix(4, 4)
    for row in range(2):
        for column in range(2):
            basis[row, column] = mp.re(2 * periods.omega[row, column])
            basis[row, 2 + column] = mp.re(
                2 * periods.omega_prime[row, column])
            basis[2 + row, column] = mp.im(
                2 * periods.omega[row, column])
            basis[2 + row, 2 + column] = mp.im(
                2 * periods.omega_prime[row, column])
    image = curve.abel_map_kind_1((branch, mp.zero))
    target = mp.matrix([
        mp.re(image[0]), mp.re(image[1]),
        mp.im(image[0]), mp.im(image[1]),
    ])
    coefficients = mp.lu_solve(basis, target)
    halves = [mp.nint(2 * value) for value in coefficients]
    residual = max(abs(2 * coefficients[i] - halves[i]) for i in range(4))
    if residual > mp.mpf(1000) * mp.eps:
        raise RuntimeError("branch Abel image is not a numerical half-period")

    def half(bit):
        return mp.mpf(int(bit) % 2) / 2

    m = tuple(half(halves[i]) for i in range(2))
    n = tuple(half(halves[2 + i]) for i in range(2))
    return n, m


def characteristic_label(chi):
    """Readable [a; b] label for a two-torsion characteristic."""
    (a1, a2), (b1, b2) = chi

    def bit(value):
        return "1" if value else "0"

    return "[%s%s;%s%s]" % (bit(a1), bit(a2), bit(b1), bit(b2))


def real_oval_signs():
    """Radical orientation fixed by the selected real oval and q(0)."""
    return {"W1": -1, "W2": -1,
            "a0": -1, "a1": -1, "a2": -1, "a3": -1, "a4": -1}


def stage_theta(curve, periods, abel_data, verbose=True):
    """Derive and verify the fifteen Rosenhain theta quotients."""
    mp.dps = 30
    tau = periods.tau
    states = abel_data["states"]
    signs = real_oval_signs()
    initial_values = p_values(states[0], {}, signs)
    score = max(identity_residuals(initial_values))
    algebraic_state = assemble_state(
        initial_values, section5_constants((1, 1, 1)))
    algebraic_error = max(
        abs(actual - expected)
        for actual, expected in zip(algebraic_state, INITIAL))
    if algebraic_error > mp.mpf(1000) * mp.eps:
        raise RuntimeError("real-oval branch convention failed initial data")
    if verbose:
        print("P-identity residual at t = 0:", mp.nstr(score, 4))
        print("direct algebraic-state residual at t = 0:",
              mp.nstr(algebraic_error, 4))
    cache = {}
    series = {slot: [] for slot in SLOTS}
    identities = []
    for state in states:
        values = p_values(state, cache, signs)
        for slot in SLOTS:
            series[slot].append(values[slot])
        identities.append(max(identity_residuals(values)))
    if verbose:
        print("max P-identity residual over the window:",
              mp.nstr(max(identities), 4))
    us = abel_data["abel"]
    inverse = mp.inverse(2 * periods.omega)
    zs = [inverse * u for u in us]
    denominator = curve.riemann_constant().characteristic
    branch_chars = [branch_characteristic(curve, periods, root)
                    for root in FIVE]
    numerators = {}
    for slot in SLOTS:
        indices = [int(value) - 1 for value in slot[1:]]
        numerators[slot] = characteristic_sum(
            denominator, *(branch_chars[index] for index in indices))
    if len(set(numerators.values())) != len(SLOTS):
        raise RuntimeError("derived P-characteristics are not distinct")

    anchor_table = {
        chi: rtheta(zs[0], tau, characteristic=chi)
        for chi in set(numerators.values()) | {denominator}}
    report = {}
    for slot in SLOTS:
        numerator = numerators[slot]
        kappa = (series[slot][0] * anchor_table[denominator]
                 / anchor_table[numerator])
        ratios = [series[slot][i]
                  * rtheta(zs[i], tau, characteristic=denominator)
                  / rtheta(zs[i], tau, characteristic=numerator)
                  for i in range(len(zs))]
        scale = max(mp.one, abs(kappa))
        error = max(abs(value - kappa) for value in ratios) / scale
        report[slot] = (numerator, denominator, kappa, error)
    score = max(report[slot][3] for slot in SLOTS)
    if max(score, max(identities), algebraic_error) > mp.mpf("1e-24"):
        raise RuntimeError("algebraic or held-out theta-quotient check failed")
    if verbose:
        print("derived branch half-characteristics:")
        for index, characteristic in enumerate(branch_chars, 1):
            print("  a%d: %s" % (index, characteristic_label(characteristic)))
        print("Riemann-constant denominator:",
              characteristic_label(denominator))
        print("max held-out theta-quotient residual:", mp.nstr(score, 4))
        for slot in SLOTS:
            chi, chi_den, kappa, error = report[slot]
            print("  %-4s = %-24s * theta%s(v)/theta%s(v)   error %s"
                  % (slot, mp.nstr(kappa, 10),
                     characteristic_label(chi),
                     characteristic_label(chi_den), mp.nstr(error, 3)))
    return {
        "name": "W(u)",
        "direction": 1,
        "offset": mp.matrix(2, 1),
        "transform": inverse,
        "velocity": inverse * abel_data["velocity"],
        "signs": signs,
        "report": report,
    }


def section5_constants(signs):
    """The constants L, M, N, L1, M1, N1, L2, M2, N2, E of section 5.

    The memoir's sigma multipliers are taken on the branch
    sigma_l(w)/sigma(w) = -i sqrt(l1 + e_l), which the P-identity
    orientation and the selected initial sheet determine. Three departures
    from the printed transcription are made, each verified to ~1e-28
    algebraically: L1, M1, N1 carry the minus sign of the sigma-product
    i*i; the gamma formula carries an R0 multiplier on its L2 term; and
    the gamma' formula takes prefactor -i/2 with R0 on the
    (e_j - e_k)^2 P_l terms.
    """
    sa, sb, sc = signs
    root_a = sa * mp.sqrt(L1 + E1)
    root_b = sb * mp.sqrt(L1 + E2)
    root_c = sc * mp.sqrt(L1 + E3)
    return {
        "E": (E2 - E3) * (E3 - E1) * (E1 - E2),
        "L": -1j * (E2 - E3) * root_a,
        "M": -1j * (E3 - E1) * root_b,
        "N": -1j * (E1 - E2) * root_c,
        "L1": -(E2 - E3) * root_b * root_c,
        "M1": -(E3 - E1) * root_c * root_a,
        "N1": -(E1 - E2) * root_a * root_b,
        "L2": 1j * (E2 ** 2 - E3 ** 2) * root_a,
        "M2": 1j * (E3 ** 2 - E1 ** 2) * root_b,
        "N2": 1j * (E1 ** 2 - E2 ** 2) * root_c,
    }


def assemble_state(values, constants):
    """The section 5 formulas: (p, q, r, g0, g1, g2) from the P's."""
    def P(name):
        return values[name]

    c = constants
    den = c["L"] * P("P1") + c["M"] * P("P2") + c["N"] * P("P3")
    x = c["L"] * P("P23") + c["M"] * P("P13") + c["N"] * P("P12")
    p = -1j * (c["L1"] * P("P1") + c["M1"] * P("P2")
               + c["N1"] * P("P3")) / den
    q = c["E"] / den
    r = -1j * x / den
    g2 = (c["L1"] * P("P23") + c["M1"] * P("P13")
          + c["N1"] * P("P12")) / den / C0
    bracket = (c["L"] ** 2 * P("P14") * P("P15")
               + c["M"] ** 2 * P("P24") * P("P25")
               + c["N"] ** 2 * P("P34") * P("P35")
               + c["M"] * c["N"] * (P("P24") * P("P35")
                                     + P("P25") * P("P34")
                                     + R0 * (E2 - E3) ** 2 * P("P1"))
               + c["N"] * c["L"] * (P("P34") * P("P15")
                                     + P("P35") * P("P14")
                                     + R0 * (E3 - E1) ** 2 * P("P2"))
               + c["L"] * c["M"] * (P("P14") * P("P25")
                                     + P("P15") * P("P24")
                                     + R0 * (E1 - E2) ** 2 * P("P3")))
    g1 = -1j * bracket / (2 * den ** 2) / C0
    g0 = (-4 * L1
          + R0 * (c["L2"] * P("P1") + c["M2"] * P("P2")
                  + c["N2"] * P("P3")) / den
          - (x / den) ** 2) / (2 * C0)
    return (p, q, r, g0, g1, g2)


@dataclass
class MotionComparison:
    """Independent RK4 and theta grids, with errors at the theta samples."""

    rk4_times: list
    rk4_states: tuple
    theta_times: list
    theta_states: list
    errors: list
    rk4_refinement: object = None
    rk4_error_ratio: object = None


def stage_demo(curve, periods, abel_data, config, verbose=True, *,
               rk4_samples=241, theta_stride=12, rk4_substeps=128):
    """Theta-function solution versus RK4 (the memoir's final formulas)."""
    mp.dps = 30
    tau = periods.tau
    report = config["report"]
    transform = config["transform"]
    direction = config["direction"]
    offset = config["offset"]
    velocity = config["velocity"]
    v0 = transform * (direction * abel_data["abel"][0] + offset)
    needed = sorted({report[slot][k] for slot in SLOTS for k in (0, 1)},
                    key=characteristics().index)

    def theta_side_values(v):
        table = {chi: rtheta(v, tau, characteristic=chi) for chi in needed}
        return {slot: report[slot][2]
                * table[report[slot][0]] / table[report[slot][1]]
                for slot in SLOTS}

    # The real-oval continuation fixes sigma_l(w)/sigma(w) to
    # -i*sqrt(l1+e_l), with the positive real square roots below.
    signs = (1, 1, 1)
    constants = section5_constants(signs)
    initial_theta = assemble_state(theta_side_values(v0), constants)
    error = max(abs(a - b) for a, b in zip(initial_theta, INITIAL))
    if verbose:
        print("initial-condition check: root signs %s, max error: %s"
              % (signs, mp.nstr(error, 4)))
        if error > mp.mpf(10) ** -15:
            print("WARNING: initial-condition error is large")

    # The theta representation is a global identity on the Jacobian, so
    # the comparison window extends well beyond the turning point of
    # s2 near t = 0.73 (unlike the algebraic verification window, which
    # must stay free of turning points for the continuously tracked P's).
    # The evaluated denominator stays nonzero on this demonstration window,
    # so the analytic side is regular across the sampled turning point.
    t_values = [mp.mpf(i) * mp.mpf("2.4") / (rk4_samples - 1)
                for i in range(rk4_samples)]
    states = rk4_trajectory(INITIAL, t_values, substeps=rk4_substeps)
    sample_indices = list(range(0, rk4_samples, theta_stride))
    if sample_indices[-1] != rk4_samples - 1:
        sample_indices.append(rk4_samples - 1)
    theta_times = [t_values[i] for i in sample_indices]
    sampled_states = [states[i] for i in sample_indices]
    analytic = []
    for t in theta_times:
        v = v0 + velocity * t
        values = theta_side_values(v)
        analytic.append(assemble_state(values, constants))
    names = ("p", "q", "r", "g0", "g1", "g2")
    errors = [max(abs(a[i] - s[i]) for a, s in zip(analytic, sampled_states))
              for i in range(6)]
    imaginary = [max(abs(mp.im(a[i])) for a in analytic) for i in range(6)]
    invariant_errors = []
    for a in analytic:
        target = (6 * L1, 2 * L, mp.zero, K_INTEGRAL ** 2)
        invariant_errors.append(max(abs(x - y)
                                    for x, y in zip(invariants(a), target)))
    ode_errors = []
    for t in (mp.mpf("0.17"), mp.mpf("0.91"), mp.mpf("1.73")):
        def analytic_state(value):
            return assemble_state(
                theta_side_values(v0 + velocity * value), constants)
        state = analytic_state(t)
        # Differentiate the vector once, sharing theta evaluations across
        # all six components instead of recomputing it for each component.
        derivative = mp.diff(lambda value: mp.matrix(analytic_state(value)), t)
        ode_errors.append(max(abs(a - b) for a, b in
                              zip(derivative, euler_poisson_rhs(state))))
    if verbose:
        print("max |theta solution - RK4| over [0, %s]:"
              % mp.nstr(t_values[-1], 3))
        for name, err in zip(names, errors):
            print("   %s: %s" % (name, mp.nstr(err, 4)))
        print("max imaginary part of the theta solution:",
              mp.nstr(max(imaginary), 4))
        print("max invariant error of the theta solution:",
              mp.nstr(max(invariant_errors), 4))
        print("max direct Euler-Poisson residual:",
              mp.nstr(max(ode_errors), 4))
    if max(error, max(imaginary), max(invariant_errors), max(ode_errors)) > mp.mpf("1e-24"):
        raise RuntimeError("theta solution failed its analytic checks")
    return MotionComparison(t_values, states, theta_times, analytic, errors)


def documentation_example(verbose=False):
    """Run the checked 30-digit example without the long RK4-only stage.

    Return independent grids and component errors. The same calculation
    drives the documentation plot and its doctest.
    """
    with mp.workdps(30):
        curve, periods = stage_curve(verbose)
        locus = sorted(curve.branch_locus.branch_values)
        if curve.genus != 2 or len(locus) != 5:
            raise RuntimeError("unexpected spectral curve topology")
        if max(abs(a - b) for a, b in zip(locus, sorted(FIVE))) > mp.mpf("1e-24"):
            raise RuntimeError("spectral branch points failed their check")
        abel_data = stage_abel(curve, periods, verbose)
        residual = max(mp.norm(u - abel_data["abel0"] - abel_data["velocity"] * t)
                       for u, t in zip(abel_data["abel"], abel_data["t_values"]))
        if residual > mp.mpf("1e-8"):
            raise RuntimeError("Abel flow failed its linearity check")
        config = stage_theta(curve, periods, abel_data, verbose)
        result = stage_demo(curve, periods, abel_data, config, verbose,
                            rk4_samples=241, theta_stride=12, rk4_substeps=128)
        coarse = rk4_trajectory(INITIAL, result.rk4_times, substeps=64)
        result.rk4_refinement = max(abs(a - b)
                                    for fine, rough in zip(result.rk4_states, coarse)
                                    for a, b in zip(fine, rough))
        coarse_error = max(abs(a - b)
                           for t, state in zip(result.theta_times, result.theta_states)
                           for a, b in zip(state, coarse[result.rk4_times.index(t)]))
        result.rk4_error_ratio = coarse_error / max(result.errors)
        if max(result.errors) > mp.mpf("1e-14"):
            raise RuntimeError("theta solution failed its RK4 comparison")
        if result.rk4_refinement > mp.mpf("1e-13"):
            raise RuntimeError("RK4 step-refinement check failed")
        if not 8 < result.rk4_error_ratio < 32:
            raise RuntimeError("RK4 comparison did not show fourth-order convergence")
        if verbose:
            print("max RK4 change on halving the step:", mp.nstr(result.rk4_refinement, 4))
            print("coarse/fine RK4 comparison error ratio:", mp.nstr(result.rk4_error_ratio, 5))
        return result


def make_figure(result):
    """Plot the checked motion; return a Matplotlib figure.

    Matplotlib is optional and is imported only when a plot is requested.
    """
    import matplotlib.pyplot as plt
    from examples._plotting import comparison_styles

    rk4_ts = [float(t) for t in result.rk4_times]
    theta_ts = [float(t) for t in result.theta_times]
    states = result.rk4_states
    analytic = result.theta_states
    names = ("p", "q", "r", "g0", "g1", "g2")
    figure, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    for ax, indices, label in ((axes[0], range(3), "Angular velocities"),
                               (axes[1], range(3, 6), "Direction cosines")):
        for style_index, index in enumerate(indices):
            line_style, marker_style = comparison_styles(style_index)
            ax.plot(rk4_ts, [float(s[index]) for s in states], **line_style,
                            label=names[index] + " (RK4)")
            ax.plot(theta_ts, [float(mp.re(a[index])) for a in analytic],
                    **marker_style, label=names[index] + " (theta)")
        ax.set_ylabel(label)
        ax.legend(fontsize=8, ncol=3)
    axes[1].set_xlabel("Time t")
    figure.suptitle("Kovalevskaya top: the original memoir's solution")
    figure.tight_layout()
    return figure


if __name__ == "__main__":
    main()
