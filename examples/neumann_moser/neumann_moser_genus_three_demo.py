#!/usr/bin/env python3
"""Compare a genus-three Kleinian solution with direct RK4 integration.

The example solves the Neumann-Moser dynamical system of P. Baron,
"The Neumann-Moser dynamical system and the Korteweg-de Vries hierarchy",
arXiv:2402.18079, in the dimension n = g = 3. The system lives on C^10 with
coordinates (u1, u2, u3, v1, v2, v3, w1, w2, w3, w4) and, with
Gamma = w1 - u1, its expanded equations (the coordinate form of the paper's
equations (20)-(22)) are

    u1' = -2 v1                v1' = -Gamma u1 - u2 + w2
    u2' = -2 v2                v2' = -Gamma u2 - u3 + w3
    u3' = -2 v3                v3' = -Gamma u3 + w4
    w1' = 2 v1                 w2' = 2 v2 + 2 Gamma v1
                                w3' = 2 v3 + 2 Gamma v2
                                w4' = 2 Gamma v3

Section 7.3 of the paper (via Theorem 5.1 and Buchstaber, "The Mumford
dynamical system and hyperelliptic Kleinian functions", arXiv:2402.09218)
solves this system by genus-three hyperelliptic Kleinian functions on the
spectral curve

    y^2 = F(x) = 4 x^7 + lambda_4 x^5 + lambda_6 x^4
                 + lambda_8 x^3 + lambda_10 x^2 + lambda_12 x + lambda_14,

with the solution flowing along the z_1 direction of the Jacobian:

    u_i = -wp_{2i},    v_i = wp_{2i}'/2,    W_xi = (xi + 2 wp_2) p_I + p_III/2

where p_I = xi^3 - wp_2 xi^2 - wp_4 xi - wp_6, ' denotes d/dz_1 and

    w1 = wp_2,
    w2 = wp_2''/2 - wp_4 - 2 wp_2^2,
    w3 = wp_4''/2 - wp_6 - 2 wp_2 wp_4,
    w4 = wp_6''/2 - 2 wp_2 wp_6.

The corrections to the printed formulas established by the genus-two demo
in this directory carry over verbatim: the W_xi factor is p_I, not p_II;
the middle w' range runs to n; and w_{n+1}' = 2 Gamma v_n. The genus-three
identities checked below additionally exercise wp_6, wp_{3,3} and wp_{3,5},
which have no genus-two counterparts.

The paper's z_1, z_3, z_5 are Genera Abelian coordinates 2, 1, 0, so
wp_2 = -d^2 log sigma / dz_1^2 maps to kleinian_p indices (2, 2), wp_4 to
(1, 2) and wp_6 to (0, 2), with ' appending a 2 and '' appending two.

The curve below has seven real branch points summing to zero, so it is in
the canonical form with vanishing x^6 coefficient. The imaginary half-period
translate omega tau [1, 1, 1]^T places the trajectory on the compact real
component of the Jacobian, as in the genus-two demo; the first two Abelian
coordinates are constant phases, of which the demo varies the first.

The script also closes the state-to-divisor-to-Jacobian loop through
Curve.abel_map_kind_1: the roots of U(x) = x^3 + u1 x^2 + u2 x + u3 give
the divisor x-coordinates, the conserved relation U*W + V^2 = F/4 gives the
sheets y_i = 2*V(x_i), and the Abel image reproduces the selected phase
modulo the period lattice and all ten state coordinates.

Run from the Genera repository root with

    .venv/bin/python -m examples.neumann_moser.neumann_moser_genus_three_demo
"""

import argparse
from dataclasses import dataclass

from genera import algebraic_curve, kleinian_p
from mpmath import mp

from examples._rk4 import rk4_step, rk4_trajectory


def multiply_by_linear(coefficients, root):
    """Multiply an ascending coefficient vector by x - root."""
    result = [mp.zero] * (len(coefficients) + 1)
    for degree, coefficient in enumerate(coefficients):
        result[degree] -= root * coefficient
        result[degree + 1] += coefficient
    return result


def curve_coefficients():
    """Return coefficients of the real genus-three demonstration curve.

    The roots sum to zero, giving the canonical form 4x^7 + lambda_4 x^5
    + lambda_6 x^4 + ... + lambda_14 with no x^6 term.
    """
    roots = (mp.mpf(-5), mp.mpf(-3), mp.mpf("-1.2"), mp.mpf("-0.3"),
             mp.mpf("0.8"), mp.mpf("2.7"), mp.mpf(6))
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
    # component on which the ten coordinates remain bounded and real
    # throughout the comparison interval.
    lattice_vector = mp.matrix([1, 1, 1])
    u_offset = (
        mp.matrix([mp.zero, mp.mpf("0.27"), mp.mpf("0.41")])
        + omega * tau * lattice_vector
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
    """Return the ten Neumann-Moser coordinates from the Kleinian solution.

    The paper's z_1 is the last Genera Abelian coordinate (the one the
    flow advances; the genus-three analogue of the second coordinate in
    the genus-two demo), z_3 is the middle coordinate and z_5 the first,
    so wp_2, wp_4 and wp_6 map to the kleinian_p coordinate index pairs
    (2, 2), (1, 2) and (0, 2) respectively.
    """
    u = [data["u_offset"][0], data["u_offset"][1], data["u_offset"][2] + x]
    wp2, wp4, wp6, wp2p, wp4p, wp6p, wp2pp, wp4pp, wp6pp = kleinian_p(
        u, data["omega"], data["tau"], data["kappa"],
        ((2, 2), (1, 2), (0, 2), (2, 2, 2), (1, 2, 2), (0, 2, 2),
         (2, 2, 2, 2), (1, 2, 2, 2), (0, 2, 2, 2)),
        data["characteristic"],
    )
    coordinates = (
        -wp2,                               # u1
        -wp4,                               # u2
        -wp6,                               # u3
        wp2p / 2,                           # v1
        wp4p / 2,                           # v2
        wp6p / 2,                           # v3
        wp2,                                # w1
        wp2pp / 2 - wp4 - 2 * wp2 ** 2,     # w2
        wp4pp / 2 - wp6 - 2 * wp2 * wp4,    # w3
        wp6pp / 2 - 2 * wp2 * wp6,          # w4
    )
    return tuple(real_part(value, f"coordinate {index + 1}")
                 for index, value in enumerate(coordinates))


def divisor_from_state(state):
    """Recover the degree-three spectral divisor from a Mumford state.

    At a root x_i of U(x)=x^3+u1*x^2+u2*x+u3, the conserved relation
    U*W+V^2=F/4 gives the affine curve point y_i=2*V(x_i).
    """
    u1, u2, u3, v1, v2, v3 = state[:6]
    roots = mp.polyroots(
        [u3, u2, u1, mp.one], maxsteps=200, asc=True)
    return tuple((root, 2 * (v1 * root ** 2 + v2 * root + v3))
                 for root in roots)


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


def kleinian_residuals(data):
    """Check the genus-three Kleinian identities at the starting point.

    These are the paper's relation (55) specialised to i = 1, 2 and 3, in
    the sigma normalisation used by `kleinian_sigma`. The wp_8 term of the
    general relation is absent at genus three, and the curve constant
    enters only the i = 1 case, at half the printed weight, exactly as in
    the genus-two demo.
    """
    u = list(data["u_offset"])
    wp2, wp4, wp6, wp33, wp35, wp2pp, wp4pp, wp6pp = kleinian_p(
        u, data["omega"], data["tau"], data["kappa"],
        ((2, 2), (1, 2), (0, 2), (1, 1), (1, 0),
         (2, 2, 2, 2), (1, 2, 2, 2), (0, 2, 2, 2)),
        data["characteristic"],
    )
    lambda4 = data["coefficients"][5]
    residual_second = wp2pp - (6 * wp2 ** 2 + 4 * wp4 + lambda4 / 2)
    residual_cross = wp4pp - (6 * (wp2 * wp4 + wp6) - 2 * wp33)
    residual_high = wp6pp - (6 * wp2 * wp6 - 2 * wp35)
    return abs(residual_second), abs(residual_cross), abs(residual_high)


def spectral_polynomial(state, data):
    """Return the coefficients of H_xi = U_xi W_xi + V_xi^2 in xi.

    For the solution above this is identically the fixed spectral curve
    F(xi)/4, the conserved Hamiltonian family of the paper's section 4
    (the coefficients h_i of H_xi are the curve's lambda constants / 4,
    including the vanishing h_1 forced by the canonical curve form).
    """
    u1, u2, u3, v1, v2, v3, w1, w2, w3, w4 = state
    upper = [u3, u2, u1, mp.one]
    middle = [v3, v2, v1]
    lower = [w4, w3, w2, w1, mp.one]
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
    expected = [value / 4 for value in data["coefficients"]]
    actual = spectral_polynomial(state, data)
    return max(abs(a - b) for a, b in zip(actual, expected))


def ode_rhs(state, data):
    """Return the first-order form of the expanded Neumann-Moser system."""
    u1, u2, u3, v1, v2, v3, w1, w2, w3, w4 = state
    gamma = w1 - u1
    return (
        -2 * v1,
        -2 * v2,
        -2 * v3,
        -gamma * u1 - u2 + w2,
        -gamma * u2 - u3 + w3,
        -gamma * u3 + w4,
        2 * v1,
        2 * v2 + 2 * gamma * v1,
        2 * v3 + 2 * gamma * v2,
        2 * gamma * v3,
    )


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


@dataclass
class DocumentationComparison:
    """Genus-three trajectories and their enforced numerical checks."""

    rk4_times: list
    rk4_states: tuple
    analytic_times: list
    analytic_states: list
    errors: tuple
    numeric_spectral_error: object
    analytic_spectral_error: object
    kleinian_residuals: tuple
    abel_residuals: tuple
    rk4_error_ratio: object
    rk4_refinement: object


def documentation_example():
    """Check one 30-digit phase over [-0.5, 0.5] using sparse P-evaluations.

    The fine RK4 run uses 32,768 steps, saving 129 plotting samples.
    Seventeen analytic states suffice for the comparison markers. A coarse
    run uses half as many steps to verify fourth-order convergence.
    """
    with mp.workdps(30):
        data = data_for_phase(problem_data(), mp.mpf("0.13"))
        if data["curve"].genus != 3 or abs(data["coefficients"][6]) > mp.mpf("1e-24"):
            raise RuntimeError("unexpected canonical Neumann-Moser curve")
        identities = kleinian_residuals(data)
        abel_residuals = abel_map_residuals(data)
        if max((*identities, *abel_residuals)) > mp.mpf("1e-24"):
            raise RuntimeError("Neumann-Moser Kleinian or Abel-map check failed")

        rk4_times = [-mp.mpf("0.5") + mp.mpf(i) / 128 for i in range(129)]
        sample_indices = list(range(0, 129, 8))
        analytic_times = [rk4_times[i] for i in sample_indices]
        analytic = [analytic_state(t, data) for t in analytic_times]
        fine = rk4_trajectory(ode_rhs, analytic[0], rk4_times,
                              substeps=256, args=(data,))
        coarse = rk4_trajectory(ode_rhs, analytic[0], rk4_times,
                                substeps=128, args=(data,))
        errors = tuple(max(abs(exact[component] - fine[index][component])
                           for index, exact in zip(sample_indices, analytic))
                       for component in range(10))
        coarse_error = max(abs(a - b)
                           for index, exact in zip(sample_indices, analytic)
                           for a, b in zip(exact, coarse[index]))
        ratio = coarse_error / max(errors)
        refinement = max(abs(a - b) for left, right in zip(fine, coarse)
                         for a, b in zip(left, right))
        analytic_spectral = max(spectral_error(state, data) for state in analytic)
        numeric_spectral = max(spectral_error(state, data) for state in fine)
        if analytic_spectral > mp.mpf("1e-24"):
            raise RuntimeError("analytic Neumann-Moser spectral polynomial check failed")
        if max((*errors, numeric_spectral)) > mp.mpf("1e-14"):
            raise RuntimeError("Neumann-Moser solution failed its RK4 comparison")
        if refinement > mp.mpf("1e-13") or not 8 < ratio < 32:
            raise RuntimeError("Neumann-Moser RK4 step-refinement check failed")
        return DocumentationComparison(
            rk4_times, fine, analytic_times, analytic, errors,
            numeric_spectral, analytic_spectral, identities, abel_residuals,
            ratio, refinement,
        )


def make_figure(result):
    """Plot the U, V and W coefficients, importing Matplotlib on demand."""
    import matplotlib.pyplot as plt
    from examples._plotting import comparison_styles

    rk4_times = [float(t) for t in result.rk4_times]
    analytic_times = [float(t) for t in result.analytic_times]
    names = ("u1", "u2", "u3", "v1", "v2", "v3", "w1", "w2", "w3", "w4")
    figure, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    for ax, indices, label in ((axes[0], range(3), "U coefficients"),
                               (axes[1], range(3, 6), "V coefficients"),
                               (axes[2], range(6, 10), "W coefficients")):
        for style_index, index in enumerate(indices):
            line_style, marker_style = comparison_styles(style_index)
            ax.plot(rk4_times, [float(s[index]) for s in result.rk4_states],
                            **line_style, label=names[index] + " (RK4)")
            ax.plot(analytic_times, [float(s[index]) for s in result.analytic_states],
                    **marker_style,
                    label=names[index] + " (Kleinian)")
        ax.set_ylabel(label)
        ax.legend(fontsize=8, ncol=len(indices), loc="lower center",
                  bbox_to_anchor=(0.5, 1.02))
    axes[2].set_xlabel("Time t")
    figure.tight_layout()
    return figure


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
            max(row[3][2] for row in rows),
            max(row[3][9] for row in rows),
            max(row[4] for row in rows),
        ))
    reference_data = data_for_phase(common_data, mp.mpf("0.13"))
    residuals = kleinian_residuals(reference_data)
    abel_residuals = abel_map_residuals(reference_data)
    print(f"precision: {mp.dps} decimal places")
    print(f"RK4 steps: {arguments.steps}")
    print("lambda_4:", mp.nstr(common_data["coefficients"][5], 12))
    print("lambda_6:", mp.nstr(common_data["coefficients"][4], 12))
    print("phase     max u1 error   max u2 error   max u3 error"
          "   max w4 error   spectral drift")
    for phase, u1_error, u2_error, u3_error, w4_error, drift in summaries:
        print(
            f"{mp.nstr(phase, 5):>6}   {mp.nstr(u1_error, 8):>12}   "
            f"{mp.nstr(u2_error, 8):>12}   {mp.nstr(u3_error, 8):>12}   "
            f"{mp.nstr(w4_error, 8):>12}   {mp.nstr(drift, 8):>14}"
        )
    print("Kleinian wp2'' identity residual:",
          mp.nstr(residuals[0], 8))
    print("Kleinian wp4'' identity residual:",
          mp.nstr(residuals[1], 8))
    print("Kleinian wp6'' identity residual:",
          mp.nstr(residuals[2], 8))
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
