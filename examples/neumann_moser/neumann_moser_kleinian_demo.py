#!/usr/bin/env python3
"""Compare a genus-two Kleinian solution with direct RK4 integration.

The example solves the Neumann-Moser dynamical system of P. Baron,
"The Neumann-Moser dynamical system and the Korteweg-de Vries hierarchy",
arXiv:2402.18079, in the smallest physical dimension n = g = 2. The system
lives on C^7 with coordinates (u1, u2, v1, v2, w1, w2, w3) and, with
Gamma = w1 - u1, its expanded equations (the coordinate form of the
paper's equations (20)-(22)) are

    u1' = -2 v1                v1' = -Gamma u1 - u2 + w2
    u2' = -2 v2                v2' = -Gamma u2 + w3
    w1' = 2 v1                 w2' = 2 v2 + 2 Gamma v1
                                w3' = 2 Gamma v2

Section 7.3 of the paper (via Theorem 5.1 and Buchstaber, "The Mumford
dynamical system and hyperelliptic Kleinian functions", arXiv:2402.09218)
solves this system by genus-two hyperelliptic Kleinian functions on the
spectral curve

    y^2 = F(x) = 4 x^5 + lambda_4 x^3 + lambda_6 x^2
                       + lambda_8 x + lambda_10,

with the solution flowing along the z_1 direction of the Jacobian:

    u_i = -wp_{2i},    v_i = wp_{2i}'/2,    W_xi = (xi + 2 wp_2) p_I + p_III/2

where p_I = xi^2 - wp_2 xi - wp_4, ' denotes d/dz_1 and

    w1 = wp_2,  w2 = wp_2''/2 - wp_4 - 2 wp_2^2,  w3 = wp_4''/2 - 2 wp_2 wp_4.

Corrections to the printed formulas, established numerically here and
cross-checked against RK4, the conserved spectral polynomial and a closed
state-to-divisor-to-Abel-map calculation:

1. Baron's section 7.3 prints W_xi = (xi + 2 wp_2) p_II + p_III/2, which
   has degree g and cannot match the degree-(g+1) generating polynomial
   W_xi. The factor must be p_I, exactly as in Buchstaber's theorem (his
   w_xi = z_xi + 2 x_xi (xi + 2 wp_2), halved into the monic
   normalisation of the coordinates above).
2. The paper's expanded system omits w_n' = 2 v_n + 2 Gamma v_{n-1} (the
   range "2 <= i <= n-1" should read "2 <= i <= n") and prints
   w_{n+1}' = 2 Gamma u_n for what the generating form gives as
   w_{n+1}' = 2 Gamma v_n.
3. In this sigma normalisation the genus-two identity (the paper's
   Corollary 7.7(a), also checked below) reads
   wp_2'' = 6 wp_2^2 + 4 wp_4 + lambda_4/2, with the curve coefficient
   lambda_4 entering at half the printed weight (the same normalisation
   the Manakov demo observes in Christiansen et al. equation (3.16)).

The curve below has five real branch points summing to zero, so it is in
the canonical form with vanishing x^4 coefficient. As in the Manakov
demo, the imaginary half-period translate omega tau [1, 1]^T places the
trajectory on the compact real component of the Jacobian, where the
Kleinian functions stay real and bounded; the constant first Abelian
coordinate is a free phase.

Run from the Genera repository root with

    .venv/bin/python -m examples.neumann_moser.neumann_moser_kleinian_demo
"""

import argparse

from genera import algebraic_curve, kleinian_p
from mpmath import mp

from examples._rk4 import rk4_step


def multiply_by_linear(coefficients, root):
    """Multiply an ascending coefficient vector by x - root."""
    result = [mp.zero] * (len(coefficients) + 1)
    for degree, coefficient in enumerate(coefficients):
        result[degree] -= root * coefficient
        result[degree + 1] += coefficient
    return result


def curve_coefficients():
    """Return coefficients of the real genus-two demonstration curve.

    The roots sum to zero, giving the canonical form 4x^5 + lambda_4 x^3
    + lambda_6 x^2 + lambda_8 x + lambda_10 with no x^4 term.
    """
    roots = (mp.mpf(-4), mp.mpf("-1.5"), mp.mpf("-0.3"),
             mp.mpf("1.3"), mp.mpf("4.5"))
    coefficients = [mp.mpf(4)]
    for root in roots:
        coefficients = multiply_by_linear(coefficients, root)
    return coefficients


def real_part(value, name):
    """Discard numerical noise after checking that a value is real."""
    tolerance = 1000 * mp.eps * max(1, abs(mp.re(value)))
    if abs(mp.im(value)) > tolerance:
        raise ValueError(f"{name} unexpectedly left the real slice: {value}")
    return mp.re(value)


def problem_data():
    """Construct the curve, periods and the Abelian starting point."""
    coefficients = curve_coefficients()
    curve = algebraic_curve(coefficients)
    first = curve.periods_kind_1()
    second = curve.periods_kind_2()
    omega, tau, kappa, characteristic = (
        first.omega, first.tau, second.kappa,
        curve.riemann_constant().characteristic)
    # This imaginary half-period translate selects the compact real
    # component on which the seven coordinates remain bounded and real
    # throughout the comparison interval.
    lattice_vector = mp.matrix([1, 1])
    u_offset = (
        mp.matrix([mp.zero, mp.mpf("0.27")])
        + omega * lattice_vector + omega * tau * lattice_vector
    )
    return {
        "coefficients": coefficients,
        "curve": curve,
        "omega": omega,
        "tau": tau,
        "kappa": kappa,
        "characteristic": characteristic,
        "u_offset": u_offset,
    }


def data_for_phase(data, first_coordinate):
    """Return a shallow copy with the selected constant first coordinate."""
    phased_data = dict(data)
    phased_data["u_offset"] = data["u_offset"].copy()
    phased_data["u_offset"][0] += first_coordinate
    return phased_data


def analytic_state(x, data):
    """Return u1, u2, v1, v2, w1, w2, w3 from the Kleinian solution.

    The paper's z_1 is the second Genera Abelian coordinate (the one the
    flow advances, exactly as in the Manakov demo) and z_3 is the first,
    so wp_2 = -d^2 log sigma / dz_1^2, wp_4 = -d^2 log sigma / dz_1 dz_3
    and wp_3,3 = -d^2 log sigma / dz_3^2 map to the kleinian_p coordinate
    indices (1, 1), (0, 1) and (0, 0) respectively.
    """
    u = [data["u_offset"][0], data["u_offset"][1] + x]
    wp2, wp4, wp2p, wp4p, wp2pp, wp4pp = kleinian_p(
        u, data["omega"], data["tau"], data["kappa"],
        ((1, 1), (0, 1), (1, 1, 1), (0, 1, 1), (1, 1, 1, 1), (0, 1, 1, 1)),
        data["characteristic"],
    )
    coordinates = (
        -wp2,                              # u1
        -wp4,                              # u2
        wp2p / 2,                          # v1
        wp4p / 2,                          # v2
        wp2,                               # w1
        wp2pp / 2 - wp4 - 2 * wp2 ** 2,    # w2
        wp4pp / 2 - 2 * wp2 * wp4,         # w3
    )
    return tuple(real_part(value, f"coordinate {index + 1}")
                 for index, value in enumerate(coordinates))


def divisor_from_state(state):
    """Recover the degree-two spectral divisor from a Mumford state.

    At a root x_i of U(x)=x^2+u1*x+u2, the conserved relation
    U*W+V^2=F/4 gives the affine curve point y_i=2*V(x_i).
    """
    u1, u2, v1, v2, unused_w1, unused_w2, unused_w3 = state
    roots = mp.polyroots([u2, u1, mp.one], maxsteps=200, asc=True)
    return tuple((root, 2 * (v1 * root + v2)) for root in roots)


def period_lattice_residual(left, right, data):
    """Return the residual after resolving left-right in the full lattice."""
    genus = data["omega"].rows
    periods = mp.matrix(genus, 2 * genus)
    periods[:, :genus] = 2 * data["omega"]
    periods[:, genus:] = 2 * data["omega"] * data["tau"]
    real_periods = mp.matrix(2 * genus)
    difference = left - right
    real_difference = mp.matrix(2 * genus, 1)
    for row in range(genus):
        real_difference[row] = mp.re(difference[row])
        real_difference[genus + row] = mp.im(difference[row])
        for column in range(2 * genus):
            real_periods[row, column] = mp.re(periods[row, column])
            real_periods[genus + row, column] = mp.im(
                periods[row, column])
    coordinates = mp.lu_solve(real_periods, real_difference)
    lattice_vector = mp.matrix([mp.nint(value) for value in coordinates])
    return mp.norm(difference - periods * lattice_vector)


def abel_map_residuals(data):
    """Close the state-to-divisor-to-Jacobian loop at the initial point."""
    initial_state = analytic_state(mp.zero, data)
    divisor = divisor_from_state(initial_state)
    curve_residual = max(abs(
        y ** 2 - mp.polyval(data["coefficients"], x, asc=True))
        for x, y in divisor)
    image = data["curve"].abel_map_kind_1(divisor, reduce=True)
    lattice_residual = period_lattice_residual(
        image, data["u_offset"], data)
    recovered_data = dict(data)
    recovered_data["u_offset"] = image
    recovered_state = analytic_state(mp.zero, recovered_data)
    state_residual = max(abs(recovered - expected)
                         for recovered, expected
                         in zip(recovered_state, initial_state))
    return curve_residual, lattice_residual, state_residual


def ode_rhs(state, data):
    """Return the first-order form of the expanded Neumann-Moser system."""
    u1, u2, v1, v2, w1, w2, w3 = state
    gamma = w1 - u1
    return (
        -2 * v1,
        -2 * v2,
        -gamma * u1 - u2 + w2,
        -gamma * u2 + w3,
        2 * v1,
        2 * v2 + 2 * gamma * v1,
        2 * gamma * v2,
    )


def spectral_polynomial(state, data):
    """Return the coefficients of H_xi = U_xi W_xi + V_xi^2 in xi.

    For the solution above this is identically the fixed spectral curve
    F(xi)/4, the conserved Hamiltonian family of the paper's section 4
    (the coefficients h_i of H_xi are the curve's lambda constants / 4,
    including the vanishing h_1 forced by the canonical curve form).
    """
    u1, u2, v1, v2, w1, w2, w3 = state
    upper = [u2, u1, mp.one]
    middle = [v2, v1]
    lower = [w3, w2, w1, mp.one]
    product = [mp.zero] * (len(upper) + len(lower) - 1)
    for i, a in enumerate(upper):
        for j, b in enumerate(lower):
            product[i + j] += a * b
    square = [mp.zero] * (2 * len(middle) - 1)
    for i, a in enumerate(middle):
        for j, b in enumerate(middle):
            square[i + j] += a * b
    for i in range(len(square)):
        product[i] += square[i]
    return product


def spectral_error(state, data):
    """Return the maximal deviation of H_xi from the fixed curve F/4."""
    expected = [value / 4 for value in data["coefficients"][:4]]
    expected += [mp.zero, mp.one]
    actual = spectral_polynomial(state, data)
    return max(abs(a - b) for a, b in zip(actual, expected))


def kleinian_residuals(data):
    """Check the genus-two Kleinian identities at the starting point.

    These are the paper's relations (55) specialised to g = 2 (its
    Corollary 7.7(a) and the i = 2 case), with the curve constants in the
    sigma normalisation used by kleinian_sigma.
    """
    u = list(data["u_offset"])
    wp2, wp4, wp33, wp2pp, wp4pp = kleinian_p(
        u, data["omega"], data["tau"], data["kappa"],
        ((1, 1), (0, 1), (0, 0), (1, 1, 1, 1), (0, 1, 1, 1)),
        data["characteristic"],
    )
    lambda4 = data["coefficients"][3]
    residual_second = wp2pp - (6 * wp2 ** 2 + 4 * wp4 + lambda4 / 2)
    residual_cross = wp4pp - (6 * wp2 * wp4 - 2 * wp33)
    return abs(residual_second), abs(residual_cross)


def compute_comparison(data, start, stop, steps, samples):
    """Return sampled analytic and RK4 trajectories."""
    if steps <= 0 or samples <= 1:
        raise ValueError("steps must be positive and samples must exceed one")
    if steps % (samples - 1):
        raise ValueError("steps must be divisible by samples - 1")
    step = (stop - start) / steps
    sample_stride = steps // (samples - 1)
    state = analytic_state(start, data)
    initial_polynomial = spectral_polynomial(state, data)
    rows = []
    for index in range(steps + 1):
        if not index % sample_stride:
            x = start + index * step
            analytic = analytic_state(x, data)
            errors = tuple(abs(numeric - exact)
                           for numeric, exact in zip(state, analytic))
            polynomial = spectral_polynomial(state, data)
            drift = max(abs(a - b) for a, b
                        in zip(polynomial, initial_polynomial))
            rows.append((x, analytic, state, errors, drift))
        if index < steps:
            state = rk4_step(ode_rhs, state, step, data)
    return rows


def parse_arguments():
    """Parse command-line options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=30)
    parser.add_argument("--steps", type=int, default=3200)
    parser.add_argument("--samples", type=int, default=101)
    parser.add_argument("--start", default="-4")
    parser.add_argument("--stop", default="4")
    parser.add_argument(
        "--phases", nargs="+", default=("-6", "0.13", "6"),
        help="constant first Abelian coordinates to compare")
    parser.add_argument("--tol", default="1e-5")
    return parser.parse_args()


def main():
    """Run the numerical comparison and print its residuals."""
    arguments = parse_arguments()
    mp.dps = arguments.dps
    common_data = problem_data()
    summaries = []
    for phase_text in arguments.phases:
        phase = mp.mpf(phase_text)
        data = data_for_phase(common_data, phase)
        rows = compute_comparison(
            data, mp.mpf(arguments.start), mp.mpf(arguments.stop),
            arguments.steps, arguments.samples)
        summaries.append((
            phase,
            max(row[3][0] for row in rows),
            max(row[3][1] for row in rows),
            max(row[3][6] for row in rows),
            max(row[4] for row in rows),
        ))
    reference_data = data_for_phase(common_data, mp.mpf("0.13"))
    residuals = kleinian_residuals(reference_data)
    abel_residuals = abel_map_residuals(reference_data)
    print(f"precision: {mp.dps} decimal places")
    print(f"RK4 steps: {arguments.steps}")
    print("lambda_4:", mp.nstr(common_data["coefficients"][3], 12))
    print("lambda_6:", mp.nstr(common_data["coefficients"][2], 12))
    print("phase     max u1 error   max u2 error   max w3 error"
          "   spectral drift")
    for phase, u1_error, u2_error, w3_error, drift in summaries:
        print(
            f"{mp.nstr(phase, 5):>6}   {mp.nstr(u1_error, 8):>12}   "
            f"{mp.nstr(u2_error, 8):>12}   {mp.nstr(w3_error, 8):>12}   "
            f"{mp.nstr(drift, 8):>14}"
        )
    print("Kleinian wp2'' identity residual:",
          mp.nstr(residuals[0], 8))
    print("Kleinian wp4'' identity residual:",
          mp.nstr(residuals[1], 8))
    print("Abel-map divisor curve residual:",
          mp.nstr(abel_residuals[0], 8))
    print("Abel-map period-lattice residual:",
          mp.nstr(abel_residuals[1], 8))
    print("Abel-map recovered-state residual:",
          mp.nstr(abel_residuals[2], 8))
    if max(max(summary[1:]) for summary in summaries) > mp.mpf(arguments.tol):
        raise SystemExit("analytic solution failed its RK4 tolerance")


if __name__ == "__main__":
    main()
