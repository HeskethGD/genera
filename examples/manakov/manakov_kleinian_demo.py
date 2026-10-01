#!/usr/bin/env python3
"""Compare a genus-two Kleinian solution with direct RK4 integration.

The example follows Christiansen, Eilbeck, Enolskii and Kostov,
"Quasi-periodic and periodic solutions for coupled nonlinear Schrodinger
equations of Manakov type", Proc. R. Soc. Lond. A 456 (2000), equations
(2.1), (3.9), (3.16), (3.17), and (3.23).

Besides the RK4 and identity checks, the demo reconstructs the spectral
divisor from one physical state and uses ``Curve.abel_map_kind_1`` to
recover the original Abelian point modulo its full period lattice.

Run from the Genera repository root with

    .venv/bin/python -m examples.manakov.manakov_kleinian_demo
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
    """Return coefficients of the real genus-two demonstration curve."""
    roots = (mp.zero, mp.mpf(1) / 8, mp.mpf(1) / 4,
             mp.mpf(1) / 2, mp.one)
    coefficients = [mp.mpf(4)]
    for root in roots:
        coefficients = multiply_by_linear(coefficients, root)
    return coefficients


def polynomial_value(coefficients, x):
    """Evaluate a polynomial whose coefficients are in ascending order."""
    return mp.fsum(coefficient * x ** degree
                   for degree, coefficient in enumerate(coefficients))


def real_part(value, name):
    """Discard numerical noise after checking that a value is real."""
    tolerance = 1000 * mp.eps * max(1, abs(mp.re(value)))
    if abs(mp.im(value)) > tolerance:
        raise ValueError(f"{name} unexpectedly left the real slice: {value}")
    return mp.re(value)


def problem_data():
    """Construct the curve, marked points, periods, and Abelian offset."""
    coefficients = curve_coefficients()
    a1 = mp.mpf(3) / 4
    a2 = mp.mpf(3) / 16
    separation = a1 - a2
    c1_squared = -polynomial_value(coefficients, a1) / separation ** 2
    c2_squared = -polynomial_value(coefficients, a2) / separation ** 2
    curve = algebraic_curve(coefficients)
    first = curve.periods_kind_1()
    second = curve.periods_kind_2()
    omega, tau, kappa, characteristic = (
        first.omega, first.tau, second.kappa,
        curve.riemann_constant().characteristic)

    # This half-period translate selects a real component on which q1**2 and
    # q2**2 remain positive throughout the comparison interval.
    lattice_vector = mp.matrix([1, 1])
    u_offset = (
        mp.matrix([mp.zero, mp.mpf("0.27")])
        + omega * lattice_vector + omega * tau * lattice_vector
    )
    return {
        "coefficients": coefficients,
        "curve": curve,
        "a1": a1,
        "a2": a2,
        "separation": separation,
        "c1_squared": c1_squared,
        "c2_squared": c2_squared,
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
    """Return q1, p1, q2, p2 from the Kleinian solution."""
    u = [data["u_offset"][0], data["u_offset"][1] + x]
    wp22, wp12, wp222, wp122 = kleinian_p(
        u, data["omega"], data["tau"], data["kappa"],
        ((1, 1), (0, 1), (1, 1, 1), (0, 1, 1)),
        data["characteristic"],
    )
    a1 = data["a1"]
    a2 = data["a2"]
    separation = data["separation"]
    q1_squared = real_part(
        2 * (a1 ** 2 - a1 * wp22 - wp12) / separation, "q1**2")
    q2_squared = real_part(
        2 * (a2 ** 2 - a2 * wp22 - wp12) / -separation, "q2**2")
    if q1_squared <= 0 or q2_squared <= 0:
        raise ValueError("the selected analytic trajectory crossed q_i**2 = 0")
    q1 = mp.sqrt(q1_squared)
    q2 = mp.sqrt(q2_squared)
    # The paper prints u2 = 2*x + b below (2.12), but later identifies x = u2
    # in (3.18)-(3.19). With its curve and differential normalization, direct
    # substitution in (2.1) selects u2 = x + b. Differentiating (3.23) then
    # gives the momenta directly in terms of the third-order P-functions.
    p1 = real_part(
        -(a1 * wp222 + wp122) / (separation * q1), "p1")
    p2 = real_part(
        (a2 * wp222 + wp122) / (separation * q2), "p2")
    return (q1, p1, q2, p2)


def divisor_from_state(state, data):
    """Recover the degree-two spectral divisor from the physical state.

    Equation (3.23) expresses q_i**2 and q_i*p_i as values of the Jacobi
    inversion polynomials at a_i. Solving the resulting two linear systems
    recovers wp22, wp12, wp222 and wp122 without another theta evaluation.
    """
    q1, p1, q2, p2 = state
    a1 = data["a1"]
    a2 = data["a2"]
    separation = data["separation"]
    even_at_a1 = a1 ** 2 - separation * q1 ** 2 / 2
    even_at_a2 = a2 ** 2 + separation * q2 ** 2 / 2
    wp22 = (even_at_a1 - even_at_a2) / separation
    wp12 = even_at_a1 - a1 * wp22
    odd_at_a1 = -separation * q1 * p1
    odd_at_a2 = separation * q2 * p2
    wp222 = (odd_at_a1 - odd_at_a2) / separation
    wp122 = odd_at_a1 - a1 * wp222
    roots = mp.polyroots(
        [-wp12, -wp22, mp.one], maxsteps=200, asc=True)
    return tuple((root, wp222 * root + wp122) for root in roots)


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
    """Close the physical-state-to-divisor-to-Jacobian loop."""
    initial_state = analytic_state(mp.zero, data)
    divisor = divisor_from_state(initial_state, data)
    curve_residual = max(abs(
        y ** 2 - polynomial_value(data["coefficients"], x))
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
    """Return the first-order form of equation (2.1)."""
    q1, p1, q2, p2 = state
    squared_norm = q1 ** 2 + q2 ** 2
    return (
        p1,
        -squared_norm * q1 + data["a1"] * q1
        + data["c1_squared"] / q1 ** 3,
        p2,
        -squared_norm * q2 + data["a2"] * q2
        + data["c2_squared"] / q2 ** 3,
    )


def hamiltonian(state, data):
    """Evaluate the conserved Hamiltonian in equation (2.2)."""
    q1, p1, q2, p2 = state
    return (
        (p1 ** 2 + p2 ** 2) / 2
        + (q1 ** 2 + q2 ** 2) ** 2 / 4
        - (data["a1"] * q1 ** 2 + data["a2"] * q2 ** 2) / 2
        + (data["c1_squared"] / q1 ** 2
           + data["c2_squared"] / q2 ** 2) / 2
    )


def fourth_order_residuals(data):
    """Check the paper's fourth-order Kleinian identities at the centre."""
    u = list(data["u_offset"])
    values = kleinian_p(
        u, data["omega"], data["tau"], data["kappa"],
        ((1, 1), (0, 1), (0, 0), (1, 1, 1, 1), (0, 1, 1, 1)),
        data["characteristic"],
    )
    wp22, wp12, wp11, wp2222, wp1222 = values
    coefficients = data["coefficients"]
    residual_316 = (
        wp2222
        - (6 * wp22 ** 2 + coefficients[3] / 2
           + coefficients[4] * wp22 + 4 * wp12)
    )
    residual_317 = (
        wp1222
        - (6 * wp22 * wp12 + coefficients[4] * wp12 - 2 * wp11)
    )
    return abs(residual_316), abs(residual_317)


def compute_comparison(data, start, stop, steps, samples):
    """Return sampled analytic and RK4 trajectories."""
    if steps <= 0 or samples <= 1:
        raise ValueError("steps must be positive and samples must exceed one")
    if steps % (samples - 1):
        raise ValueError("steps must be divisible by samples - 1")
    step = (stop - start) / steps
    sample_stride = steps // (samples - 1)
    state = analytic_state(start, data)
    initial_hamiltonian = hamiltonian(state, data)
    rows = []
    for index in range(steps + 1):
        if not index % sample_stride:
            x = start + index * step
            analytic = analytic_state(x, data)
            errors = tuple(abs(numeric - exact)
                           for numeric, exact in zip(state, analytic))
            rows.append((x, analytic, state, errors,
                         abs(hamiltonian(state, data)
                             - initial_hamiltonian)))
        if index < steps:
            state = rk4_step(ode_rhs, state, step, data)
    return rows


def parse_arguments():
    """Parse command-line options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=30)
    parser.add_argument("--steps", type=int, default=3200)
    parser.add_argument("--samples", type=int, default=201)
    parser.add_argument("--start", default="-4")
    parser.add_argument("--stop", default="4")
    parser.add_argument(
        "--phases", nargs="+", default=("-6", "0.13", "6"),
        help="constant first Abelian coordinates to compare")
    parser.add_argument("--tol", default="1e-6")
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
            max(row[3][2] for row in rows),
            max(row[4] for row in rows),
        ))
    reference_data = data_for_phase(common_data, mp.mpf("0.13"))
    residuals = fourth_order_residuals(reference_data)
    abel_residuals = abel_map_residuals(reference_data)
    print(f"precision: {mp.dps} decimal places")
    print(f"RK4 steps: {arguments.steps}")
    print(f"C1^2: {mp.nstr(common_data['c1_squared'], 12)}")
    print(f"C2^2: {mp.nstr(common_data['c2_squared'], 12)}")
    print("phase       max q1 error   max q2 error   Hamiltonian drift")
    for phase, q1_error, q2_error, energy_error in summaries:
        print(
            f"{mp.nstr(phase, 5):>7}   {mp.nstr(q1_error, 8):>12}   "
            f"{mp.nstr(q2_error, 8):>12}   {mp.nstr(energy_error, 8):>17}"
        )
    print("equation (3.16) residual:", mp.nstr(residuals[0], 8))
    print("equation (3.17) residual:", mp.nstr(residuals[1], 8))
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
