#!/usr/bin/env python3
"""Solve 9-dimensional Reissner-Nordstrom-de Sitter geodesics by Kleinian
inversion, and compare with direct RK4 integration.

The example reproduces the application section of V. Z. Enolskii, E.
Hackmann, V. Kagramanova, J. Kunz and C. Laemmerzahl, "Inversion of
hyperelliptic integrals of arbitrary genus with application to particle
motion in General Relativity", arXiv:1011.6459 (project update 16). A
massive test particle (delta = 1) moves in the nine-dimensional
Reissner-Nordstrom-de Sitter metric with dimensionless parameters

    Lambda-tilde = 8.7e-5,   q-tilde = 0.4,   E^2 = 1.045,
    L-tilde^2 = 0.25,

for which the radial equation reads

    (dr/dphi)^2 = R16(r)/r^10,   L^2 R16(r) = r^16 Lambda/28
        + (E^2 - 1 + L^2 Lambda/28) r^14 - L^2 r^12 + r^8
        + L^2 r^6 - q^2 r^2 - q^2 L^2.

The chain of substitutions r = 1/sqrt(u), then u = 1/x + u8 with u8 any
root of the degree-eight polynomial R8, produces the genus-three spectral
curve

    y^2 = P7(x) = (4/L^2) x^8 R8(1/x + u8) = b7 prod(x - e_i)

whose eighth branch point is sent to infinity, and the azimuthal angle
becomes the third holomorphic integral,

    phi - phi_in = integral x^2 dx / y.

The inversion formula of the paper (its equations (4.13) and (6.18))
then gives the solution on the first stratum of the theta divisor:

    x = -sigma_13(u) / sigma_23(u),
    r(phi) = 1/sqrt(u8 - sigma_23(u)/sigma_13(u)),

where u = A_i + (f1(phi), f2(phi), phi - phi_in)^T, A_i is the half
period equal to the Abel image of the starting branch point, and f1, f2
are fixed at each phi by the stratum conditions

    sigma(u) = 0,   sigma_3(u) = 0,

solved here by Newton iteration warm-started along the trajectory.

Conventions and findings established in this folder:

- The paper's differential basis du_i = x^(i-1) dx/y is the ascending
  Genera ordering, so its u_1, u_2, u_3 are Genera coordinates 0, 1, 2:
  the phi flow advances the last coordinate, sigma_3 is the sigma-jet
  key (0, 0, 1), and sigma_13, sigma_23 are the keys (1, 0, 1), (0, 1, 1).
- Project update 16 flagged the printed quotient -sigma_23/sigma_13
  inside the square root as an apparent transposition of x = -sigma_13/
  sigma_23. Reading the substitution chain resolves it: u = 1/x + u8,
  so r^-2 = u = u8 - sigma_23/sigma_13 exactly as printed; there is no
  typo. Both candidate quotients are still checked numerically below.
- The solution uses ordinary sigma derivatives only, evaluated through
  the public kleinian_sigma_jet, which remains defined on the theta
  divisor where the logarithmic Kleinian functions are singular.
- The starting half period A_i is the raw reduced Abel image of the
  branch point computed by Curve.abel_map_kind_1; the paper's branch
  point identities e_i = -sigma_1/sigma_2(A_i) and x = -sigma_13/
  sigma_23(A_i) = e_i both hold to rounding level, confirming the
  convention alignment with no further shift.
- With the root u8 = -0.6074 chosen below, the five real roots of P7
  are e_1 < e_2 < e_3 < e_4 < e_7 with one conjugate pair, as in the
  paper; the labels of the paper's orbits shift because its u8 choice
  differs. The bound ("many-world") orbit runs between e_2 and e_3 and
  starts at e_2; the escape orbit starts at e_4 and approaches e_7 as
  r tends to infinity.

Run from the Genera repository root with

    .venv/bin/python -m examples.rn_desitter_9d.rn_desitter_geodesic_demo
"""

import argparse
import time

from genera import algebraic_curve, kleinian_sigma_jet
from mpmath import mp, mpf, polyroots, quad, sqrt

from examples._rk4 import rk4_step


LAMBDA = mpf("8.7e-5")
CHARGE_SQUARED = mpf("0.4") ** 2
ENERGY_SQUARED = mpf("1.045")
L_SQUARED = mpf("0.25")


def radial_potential_prime(r):
    """Return dV/dr for (dr/dphi)^2 = V(r) = R16(r)/r^10."""
    return (6 * LAMBDA * r ** 5 / 28
            + 4 * (ENERGY_SQUARED - 1 + L_SQUARED * LAMBDA / 28) * r ** 3
            - 2 * L_SQUARED * r - 2 * r ** -3 - 4 * L_SQUARED * r ** -5
            + 8 * CHARGE_SQUARED * r ** -9
            + 10 * L_SQUARED * CHARGE_SQUARED * r ** -11) / L_SQUARED


def r8_coefficients():
    """Ascending coefficients of R8 in (du/dphi)^2 = 4 R8(u)/L^2."""
    return [
        LAMBDA / 28,
        ENERGY_SQUARED - 1 + L_SQUARED * LAMBDA / 28,
        -L_SQUARED,
        mp.zero,
        mp.one,
        L_SQUARED,
        mp.zero,
        -CHARGE_SQUARED,
        -L_SQUARED * CHARGE_SQUARED,
    ]


def p7_coefficients(r8, u8):
    """Ascending coefficients of P7(x) = (4/L^2) x^8 R8(1/x + u8)."""
    coefficients = [mp.zero] * 8
    for i, a in enumerate(r8):
        for j in range(i + 1):
            power = 8 - i + j
            if power <= 7:
                coefficients[power] += a * mp.binomial(i, j) * u8 ** j
    return [(4 / L_SQUARED) * value for value in coefficients]


def real_part(value, name, tolerance):
    """Discard numerical noise after checking that a value is real.

    The check is relative to the magnitude of the value. Convention
    errors (wrong quotient, wrong branch) produce imaginary parts
    comparable to the value itself, so they trip the guard at any
    tolerance; genuine arithmetic noise instead grows along the escape
    orbit, where the cancellation u = 1/x + u8 tends to zero as r grows
    and amplifies the theta-truncation error.
    """
    if abs(mp.im(value)) > tolerance * max(mp.one, abs(mp.re(value))):
        raise ValueError(f"{name} unexpectedly left the real slice: {value}")
    return mp.re(value)


def problem_data():
    """Build the curve, periods, branch data and both starting points."""
    r8 = r8_coefficients()
    u_roots = sorted(polyroots(r8, maxsteps=200, asc=True),
                     key=lambda z: (mp.re(z), mp.im(z)))
    real_roots = [u for u in u_roots if abs(mp.im(u)) < mpf("1e-15")]
    if len(real_roots) != 6:
        raise ValueError(f"expected six real R8 roots, found {len(real_roots)}")
    positive = sorted(u for u in real_roots if u > 0)
    if len(positive) != 3:
        raise ValueError("expected three positive R8 roots")
    # Choose a moderate negative root as u8: the eighth branch point maps
    # to x = infinity and the physical turning points stay at moderate x.
    u8 = min((u for u in real_roots if u < -mpf("0.01")), key=abs)
    p7 = p7_coefficients(r8, u8)
    x_roots = polyroots(p7, maxsteps=200, asc=True)
    curve = algebraic_curve(p7)
    first = curve.periods_kind_1()
    second = curve.periods_kind_2()
    omega, tau, kappa, characteristic = (
        first.omega, first.tau, second.kappa,
        curve.riemann_constant().characteristic)

    # Turning radii: the bound orbit oscillates in the u interval between
    # the two larger positive roots; the escape orbit starts at the
    # smallest positive root and runs to u = 0.
    u_bound_inner, u_bound_outer, u_escape = positive[2], positive[1], positive[0]
    x_bound = 1 / (u_bound_inner - u8)
    x_escape = 1 / (u_escape - u8)
    r_bound = 1 / sqrt(u_bound_inner)
    r_escape = 1 / sqrt(u_escape)

    # Azimuthal ranges: one full bound radial period, and the escape
    # range up to a large fraction of the asymptotic angle. The values
    # are of order one and the quadrature of the real integrand carries
    # a few ulps of imaginary noise, so a 1e-12 relative guard is ample.
    phi_half = real_part(
        quad(lambda u: sqrt(L_SQUARED / (4 * mp.polyval(
            r8, u, asc=True))),
             [u_bound_outer, u_bound_inner]), "phi_half", mpf(10) ** -12)
    phi_inf = real_part(
        quad(lambda u: sqrt(L_SQUARED / (4 * mp.polyval(
            r8, u, asc=True))),
             [mp.zero, u_escape]), "phi_inf", mpf(10) ** -12)

    # Cauchy, event and cosmological horizon radii: positive real roots
    # of g_tt = 0, that is
    # q^2 - r^6 + r^12 - Lambda r^14/28 = 0.
    horizon_polynomial = [CHARGE_SQUARED, mp.zero, mp.zero, mp.zero,
                          mp.zero, mp.zero, -mp.one, mp.zero, mp.zero,
                          mp.zero, mp.zero, mp.zero, mp.one, mp.zero,
                          -LAMBDA / 28]
    horizons = sorted(r for r in polyroots(
        horizon_polynomial, maxsteps=200, asc=True)
                      if abs(mp.im(r)) < mpf("1e-15") and mp.re(r) > 0)

    return {
        "r8": r8,
        "u8": u8,
        "p7": p7,
        "curve": curve,
        "x_roots": x_roots,
        "omega": omega,
        "tau": tau,
        "kappa": kappa,
        "characteristic": characteristic,
        "u_bound_inner": u_bound_inner,
        "u_bound_outer": u_bound_outer,
        "x_bound": x_bound,
        "x_escape": x_escape,
        "r_bound": r_bound,
        "r_escape": r_escape,
        "phi_half": phi_half,
        "phi_inf": phi_inf,
        "horizons": horizons,
    }


def sigma_jet(data, u):
    """Evaluate the order-two Kleinian sigma jet at Abelian point u."""
    return kleinian_sigma_jet(
        u, data["omega"], data["tau"], data["kappa"], 2,
        data["characteristic"])


def stratum_point(data, start, t, f_init, tolerance):
    """Solve sigma(u) = sigma_3(u) = 0 for the stratum point at time t.

    The point is u = start + (f1, f2, t): the third Abelian coordinate
    carries the azimuthal angle and the first two are unknowns fixed by
    the two stratum conditions. Warm-started Newton iterations follow
    the trajectory; the order-two jet supplies the residuals, the
    Jacobian and the solution quotient in one evaluation.
    """
    f1, f2 = f_init
    jet = None
    evaluations = 0
    for _ in range(80):
        u = [start[0] + f1, start[1] + f2, start[2] + t]
        jet = sigma_jet(data, u)
        evaluations += 1
        residual = jet[(0, 0, 0)]
        constraint = jet[(0, 0, 1)]
        if abs(residual) + abs(constraint) < tolerance:
            break
        jacobian = mp.matrix([[jet[(1, 0, 0)], jet[(0, 1, 0)]],
                              [jet[(1, 0, 1)], jet[(0, 1, 1)]]])
        delta = mp.lu_solve(
            jacobian, -mp.matrix([residual, constraint]))
        f1 += delta[0]
        f2 += delta[1]
    else:
        raise ValueError(f"stratum Newton did not converge at t = {t}")
    # Differentiate sigma = sigma_3 = 0 along the trajectory. This
    # tangent predicts the next pair (f1, f2), leaving Newton to make a
    # small correction rather than restarting at the preceding point.
    jacobian = mp.matrix([[jet[(1, 0, 0)], jet[(0, 1, 0)]],
                          [jet[(1, 0, 1)], jet[(0, 1, 1)]]])
    tangent = mp.lu_solve(
        jacobian, -mp.matrix([jet[(0, 0, 1)], jet[(0, 0, 2)]]))
    return f1, f2, jet, tangent, evaluations


def analytic_orbit(data, x_start, phi_values, tolerance):
    """Return r(phi) from the stratum inversion plus residual maxima."""
    start = data["curve"].abel_map_kind_1((x_start, 0), reduce=True)
    f1 = f2 = mp.zero
    rows = []
    residuals = []
    worst_stratum = mp.zero
    worst_wrong = mp.zero
    tangent = None
    previous_t = None
    evaluations = 0
    for position, t in enumerate(phi_values):
        # At the branch point the local parameter is singular, so use
        # one ordinary warm start before enabling tangent prediction.
        if position >= 2:
            delta_t = t - previous_t
            f1 += delta_t * tangent[0]
            f2 += delta_t * tangent[1]
        f1, f2, jet, tangent, used = stratum_point(
            data, start, t, (f1, f2), tolerance)
        evaluations += used
        stratum_residual = (
            abs(jet[(0, 0, 0)]) + abs(jet[(0, 0, 1)]))
        worst_stratum = max(worst_stratum, stratum_residual)
        x = -jet[(1, 0, 1)] / jet[(0, 1, 1)]
        radius = 1 / sqrt(1 / x + data["u8"])
        if position:
            # The transposed reading, r = 1/sqrt(x + u8), is a genuinely
            # different function; its discrepancy from the solution shows
            # the two printed quotients are not interchangeable.
            worst_wrong = max(
                worst_wrong, abs(1 / sqrt(x + data["u8"]) - radius))
        # Relative guard (see real_part): the escape end amplifies the
        # theta-truncation noise through the 1/x + u8 cancellation, so the
        # imaginary residue can reach ~1e-9 relative to r at this working
        # precision while a convention error would be O(1) relative.
        rows.append((t, real_part(radius, f"r({t})", mpf(10) ** -6)))
        residuals.append((t, stratum_residual))
        previous_t = t
    return (start, rows, residuals, worst_stratum, worst_wrong,
            evaluations)


def rk4_orbit(r_start, phi_range, steps):
    """Integrate r'' = V'(r)/2, r'(0) = 0, by classical RK4."""
    def rhs(state):
        radius, velocity = state
        return velocity, radial_potential_prime(radius) / 2

    step = phi_range / steps
    state = (r_start, mp.zero)
    rows = []
    for index in range(steps + 1):
        rows.append((index * step, state[0]))
        if index < steps:
            state = rk4_step(rhs, state, step)
    return rows


def starting_residuals(data, start, x_start):
    """Check the branch-point identities at the stratum starting point."""
    jet = sigma_jet(data, [start[0], start[1], start[2]])
    divisor = abs(jet[(0, 0, 0)]) + abs(jet[(0, 0, 1)])
    branch = abs(-jet[(1, 0, 0)] / jet[(0, 1, 0)] - x_start)
    inversion = abs(-jet[(1, 0, 1)] / jet[(0, 1, 1)] - x_start)
    return divisor, branch, inversion


def parse_arguments():
    """Parse command-line options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=25)
    parser.add_argument("--samples", type=int, default=31)
    parser.add_argument("--steps", type=int, default=3000)
    parser.add_argument(
        "--newton-tol", default="1e-18",
        help="absolute tolerance for sigma=0 and sigma_3=0")
    parser.add_argument("--tol", default="1e-5",
                        help="maximum analytic-versus-RK4 radial error")
    return parser.parse_args()


def main():
    """Invert both geodesics and compare them with direct RK4."""
    arguments = parse_arguments()
    if arguments.samples < 3 or arguments.steps < 1:
        raise ValueError("samples must be at least 3 and steps positive")
    if arguments.steps % (arguments.samples - 1):
        raise ValueError("steps must be divisible by samples - 1")
    mp.dps = arguments.dps
    tolerance = mpf(arguments.newton_tol)
    if not 0 < tolerance < 1:
        raise ValueError("newton-tol must lie between zero and one")
    started = time.perf_counter()
    data = problem_data()

    phi_period = 2 * data["phi_half"]
    escape_range = mpf("0.85") * data["phi_inf"]
    bound_numeric = rk4_orbit(
        data["r_bound"], phi_period, arguments.steps)
    escape_numeric = rk4_orbit(
        data["r_escape"], escape_range, arguments.steps)
    stride = arguments.steps // (arguments.samples - 1)
    sample_indices = [k * stride for k in range(arguments.samples)]
    bound_grid = [bound_numeric[index][0] for index in sample_indices]
    escape_grid = [escape_numeric[index][0] for index in sample_indices]

    (bound_start, bound_rows, bound_residuals, bound_stratum, bound_wrong,
     bound_evaluations) = analytic_orbit(
        data, data["x_bound"], bound_grid, tolerance)
    (escape_start, escape_rows, escape_residuals, escape_stratum,
     escape_wrong, escape_evaluations) = analytic_orbit(
        data, data["x_escape"], escape_grid, tolerance)

    bound_error = max(
        abs(r - bound_numeric[index][1])
        for (t, r), index in zip(bound_rows, sample_indices))
    escape_error = max(
        abs(r - escape_numeric[index][1])
        for (t, r), index in zip(escape_rows, sample_indices))

    bound_identity = starting_residuals(
        data, bound_start, data["x_bound"])
    escape_identity = starting_residuals(
        data, escape_start, data["x_escape"])

    print(f"precision: {mp.dps} decimal places")
    print(f"RK4 steps: {arguments.steps}, Kleinian samples: "
          f"{arguments.samples}")
    print("Newton tolerance:", mp.nstr(tolerance, 4),
          " sigma-jet evaluations:",
          bound_evaluations + escape_evaluations)
    print("parameters: Lambda=8.7e-5, q=0.4, E^2=1.045, L^2=0.25, massive")
    print("u8:", mp.nstr(data["u8"], 10))
    print("horizons:", [mp.nstr(h, 8) for h in data["horizons"]])
    print("bound: start x =", mp.nstr(data["x_bound"], 10),
          " r =", mp.nstr(data["r_bound"], 8),
          " phi period =", mp.nstr(phi_period, 10))
    print("bound: perihelion shift per radial period:",
          mp.nstr(phi_period - 2 * mp.pi, 10))
    print("escape: start x =", mp.nstr(data["x_escape"], 10),
          " r =", mp.nstr(data["r_escape"], 8),
          " phi_inf =", mp.nstr(data["phi_inf"], 10))
    print("bound starting residuals: sigma-divisor",
          mp.nstr(bound_identity[0], 6), " branch",
          mp.nstr(bound_identity[1], 6), " inversion",
          mp.nstr(bound_identity[2], 6))
    print("escape starting residuals: sigma-divisor",
          mp.nstr(escape_identity[0], 6), " branch",
          mp.nstr(escape_identity[1], 6), " inversion",
          mp.nstr(escape_identity[2], 6))
    print("bound: worst stratum residual:", mp.nstr(bound_stratum, 6))
    print("bound: max |r analytic - r RK4|:", mp.nstr(bound_error, 6))
    print("bound: r range:",
          mp.nstr(min(r for _, r in bound_rows), 8), "to",
          mp.nstr(max(r for _, r in bound_rows), 8))
    print("bound: transposed-quotient discrepancy:",
          mp.nstr(bound_wrong, 6))
    print("escape: worst stratum residual:", mp.nstr(escape_stratum, 6))
    print("escape: max |r analytic - r RK4|:", mp.nstr(escape_error, 6))
    print("escape: r reaches",
          mp.nstr(escape_rows[-1][1], 8))
    print("escape: transposed-quotient discrepancy:",
          mp.nstr(escape_wrong, 6))
    print(f"elapsed: {time.perf_counter() - started:.1f} seconds")
    if max(bound_error, escape_error) > mpf(arguments.tol):
        raise SystemExit("Kleinian inversion failed its RK4 tolerance")


if __name__ == "__main__":
    main()
