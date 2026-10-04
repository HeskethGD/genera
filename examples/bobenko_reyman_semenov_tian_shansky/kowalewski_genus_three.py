#!/usr/bin/env python3
"""Check the genus-three Kowalewski formula of Bobenko et al."""

import argparse
from dataclasses import dataclass

from generapy import rtheta, rtheta_jet
from mpmath import mp

from examples._rk4 import rk4_trajectory


mp.dps = 30

tau = mp.matrix([
    [mp.mpc("-1.4167975950708838", "0.7850043126602884"),
     mp.mpc("-0.3454058539987599", "0.0108719647228774"),
     mp.mpc("1.5832024049291242", "-0.7145033905769181")],
    [mp.mpc("-0.3454058539987745", "0.0108719647228998"),
     mp.mpc("0.0856784508858941", "0.0477513704373007"),
     mp.mpc("-0.3454058539987747", "0.0108719647228998")],
    [mp.mpc("1.5832024049291098", "-0.7145033905768952"),
     mp.mpc("-0.3454058539987707", "0.0108719647229234"),
     mp.mpc("-1.4167975950708986", "0.7850043126603116")],
])
# The independently integrated entries differ across the diagonal by about
# 2e-14.  Project the external numerical data onto its exact symmetric form
# before asking rtheta to validate it at 30-digit working precision.
tau = (tau + tau.T) / 2
velocity = mp.matrix([
    mp.mpc("-0.0496797438629430", "0.0048411751083748"),
    mp.mpc("0.0032224203740264", "-0.0195489988450696"),
    mp.mpc("-0.0496797438629420", "0.0048411751083599"),
])
r_from_abel = mp.matrix([
    mp.mpc(0, "0.5528062206626531"),
    0,
    mp.mpc(0, "-0.5528062206626164"),
])
r_vector = r_from_abel
zero_plus = mp.matrix([
    mp.mpc("-0.3616583527767097", "0.3508874460882212"),
    mp.mpc("-0.1780294655302932", "0.1264701281265881"),
    mp.mpc("-0.3616583527767367", "-0.2019187745750049"),
])
zero_minus = mp.matrix([
    mp.mpc("0.3616583527767067", "0.2019187745744802"),
    mp.mpc("0.1780294655303256", "-0.1264701281266026"),
    mp.mpc("0.3616583527767618", "-0.3508874460877124"),
])
third_scalar = mp.mpc(0, "0.19778648270564")
epsilon = ([0, 0, 0], [0, mp.mpf("0.5"), 0])
prym_tau = mp.matrix([
    [mp.mpc("0.3328096197164521", "0.1410018441667864"),
     mp.mpc("-0.6908117079975306", "0.0217439294458008")],
    [mp.mpc("-0.6908117079975490", "0.0217439294457996"),
     mp.mpc("0.0856784508858941", "0.0477513704373007")],
])
prym_tau = (prym_tau + prym_tau.T) / 2
# The difference coordinate has index-two lattice spacing in the split sum.
# Thus the elliptic theta in (7.64) maps to Generapy with twice the normalized
# b1 period of du1-du3 (the paper calls half of B0 the elliptic period).
b0 = mp.matrix([[2 * (tau[0, 0] - tau[2, 0])]])

def theta_data(argument, characteristic=None):
    """Return theta and its directional derivative along the flow."""
    jet = rtheta_jet(argument, tau, 1, characteristic)
    gradient = mp.matrix([jet[(1, 0, 0)], jet[(0, 1, 0)],
                          jet[(0, 0, 1)]])
    return jet[(0, 0, 0)], sum(
        velocity[index] * gradient[index] for index in range(3))


def solution(time, p1, p2):
    """Evaluate the six variables in Theorem 7.7, equation (7.42)."""
    phase = mp.matrix([p1 / 2, p2, p1 / 2]) + velocity * time
    theta, dtheta = theta_data(phase)
    theta_e, dtheta_e = theta_data(phase, epsilon)
    theta_r = rtheta(phase - r_vector, tau)
    theta_er = rtheta(phase - r_vector, tau, epsilon)

    a_value = rtheta(zero_minus + phase, tau)
    b_value = rtheta(zero_plus + phase, tau)
    # Both marked paths differ from the generated paths by an odd B2
    # coefficient.  For epsilon this supplies the relative multiplier -1.
    c_value = -rtheta(zero_minus + phase, tau, epsilon)
    d_value = -rtheta(zero_plus + phase, tau, epsilon)
    denominator = a_value * d_value + b_value * c_value

    l_minus = 2j * third_scalar * theta_r / theta
    l_plus = 2j * third_scalar * theta_er / theta_e
    ell1 = (l_plus + l_minus) / 2
    ell2 = (l_plus - l_minus) / (2j)
    ell3 = -1j * (dtheta_e / theta_e - dtheta / theta)

    g_minus = 2 * theta_e / theta * a_value * b_value / denominator
    g_plus = 2 * theta / theta_e * c_value * d_value / denominator
    g1 = (g_plus + g_minus) / 2
    g2 = (g_plus - g_minus) / (2j)
    g3 = -(a_value * d_value - b_value * c_value) / denominator
    return (ell1, ell2, ell3), (g1, g2, g3)


def invariants(ell, gravity):
    """Return |g|^2, H, I1 and I2 for the classical top."""
    ell1, ell2, ell3 = ell
    g1, g2, g3 = gravity
    norm = g1**2 + g2**2 + g3**2
    hamiltonian = (ell1**2 + ell2**2 + 2 * ell3**2) / 2 - g1
    i1 = (ell1 * g1 + ell2 * g2 + ell3 * g3)**2
    i2 = (ell1**2 - ell2**2 + 2 * g1)**2 \
        + 4 * (ell1 * ell2 + g2)**2
    return norm, hamiltonian, i1, i2


def euler_poisson_rhs(state):
    """Return the classical Kowalewski Euler--Poisson vector field."""
    ell1, ell2, ell3, gravity1, gravity2, gravity3 = state
    return (
        ell2 * ell3,
        -ell1 * ell3 - gravity3,
        gravity2,
        2 * gravity2 * ell3 - gravity3 * ell2,
        gravity3 * ell1 - 2 * gravity1 * ell3,
        gravity1 * ell2 - gravity2 * ell1,
    )


def rk4_comparison(stop, steps, samples):
    """Compare equation (7.42) with a direct Euler--Poisson integration."""
    if steps <= 0 or samples <= 1 or steps % (samples - 1):
        raise ValueError("steps must be positive and divisible by samples - 1")
    phase1 = mp.mpc("-0.86020578672750655767", "0.37129178949707584365")
    phase2 = mp.mpc("-0.82333617811367142382", "0.1887005817208524906")
    times = tuple(stop * index / (samples - 1) for index in range(samples))
    analytic = []
    for time in times:
        ell, gravity = solution(time, phase1, phase2)
        analytic.append(tuple(ell) + tuple(gravity))
    numeric = rk4_trajectory(
        euler_poisson_rhs, analytic[0], times,
        substeps=steps // (samples - 1),
    )
    component_errors = tuple(
        max(abs(exact[index] - approximation[index])
            for exact, approximation in zip(analytic, numeric))
        for index in range(6)
    )
    imaginary = max(abs(mp.im(value))
                    for state in analytic for value in state)
    invariant_drift = max(
        max(abs(value - reference)
            for value, reference in zip(
                invariants(state[:3], state[3:]),
                invariants(analytic[0][:3], analytic[0][3:])))
        for state in analytic
    )
    return component_errors, imaginary, invariant_drift


def decomposition_residuals(p1, p2):
    """Check the zero and epsilon cases of the decomposition (7.64)."""
    phase = mp.matrix([p1 / 2, p2, p1 / 2])
    w = [p1, p2]
    genus1_zero = rtheta([0], b0)
    genus1_half = rtheta(
        [0], b0, ([mp.mpf("0.5")], [0]))
    zero_rhs = (
        rtheta(w, prym_tau) * genus1_zero
        + rtheta(
            w, prym_tau,
            ([mp.mpf("0.5"), 0], [0, 0])) * genus1_half
    )
    epsilon_rhs = (
        rtheta(
            w, prym_tau,
            ([0, 0], [0, mp.mpf("0.5")])) * genus1_zero
        + rtheta(
            w, prym_tau,
            ([mp.mpf("0.5"), 0], [0, mp.mpf("0.5")])) * genus1_half
    )
    zero_lhs = rtheta(phase, tau)
    epsilon_lhs = rtheta(phase, tau, epsilon)
    return (
        abs(zero_lhs - zero_rhs) / max(1, abs(zero_lhs)),
        abs(epsilon_lhs - epsilon_rhs) / max(1, abs(epsilon_lhs)),
    )


def identity_760_ratio(time, p1, p2):
    """Return the phase-dependent ratio that (7.60) says is constant."""
    phase = mp.matrix([p1 / 2, p2, p1 / 2]) + velocity * time
    a_value = rtheta(zero_minus + phase, tau)
    b_value = rtheta(zero_plus + phase, tau)
    c_value = -rtheta(zero_minus + phase, tau, epsilon)
    d_value = -rtheta(zero_plus + phase, tau, epsilon)
    return ((a_value * d_value + b_value * c_value)
            / (rtheta(phase, tau)
               * rtheta(phase, tau, epsilon)))


@dataclass
class DocumentationComparison:
    """Checked sparse theta solution and dense independent RK4 trajectory."""

    rk4_times: tuple
    rk4_states: tuple
    analytic_times: tuple
    analytic_states: tuple
    errors: tuple
    imaginary: object
    analytic_invariant_drift: object
    numeric_invariant_drift: object
    decomposition_errors: tuple
    identity_error: object
    rk4_refinement: object
    rk4_refinement_ratio: object


def documentation_example():
    """Check the rounded genus-three fixture at 30-digit working precision.

    Save 257 RK4 states on [0, 2], evaluating theta at only 17 times.
    Three independent RK4 grids separate step convergence from fixture error.
    """
    with mp.workdps(30):
        p1 = mp.mpc("-0.86020578672750655767", "0.37129178949707584365")
        p2 = mp.mpc("-0.82333617811367142382", "0.1887005817208524906")
        times = tuple(mp.mpf(i) / 128 for i in range(257))
        indices = tuple(range(0, 257, 16))
        analytic_times = tuple(times[i] for i in indices)
        analytic = tuple(tuple(ell) + tuple(gravity)
                         for ell, gravity in
                         (solution(t, p1, p2) for t in analytic_times))
        fine = rk4_trajectory(euler_poisson_rhs, analytic[0], times, substeps=32)
        middle = rk4_trajectory(euler_poisson_rhs, analytic[0], times, substeps=16)
        coarse = rk4_trajectory(euler_poisson_rhs, analytic[0], times, substeps=8)

        def distance(left, right):
            return max(abs(a - b) for s, t in zip(left, right)
                       for a, b in zip(s, t))

        refinement = distance(fine, middle)
        ratio = distance(middle, coarse) / refinement
        errors = tuple(max(abs(state[j] - fine[i][j])
                           for i, state in zip(indices, analytic))
                       for j in range(6))
        imaginary = max(abs(mp.im(value)) for state in analytic for value in state)
        initial_integrals = invariants(analytic[0][:3], analytic[0][3:])

        def invariant_drift(states):
            return max(abs(value - reference) for state in states
                       for value, reference in zip(
                           invariants(state[:3], state[3:]), initial_integrals))

        analytic_drift = invariant_drift(analytic)
        numeric_drift = invariant_drift(fine)
        phases = ((mp.mpc("0.13", "0.07"), mp.mpc("-0.09", "0.04")),
                  (mp.mpc("0.21", "-0.03"), mp.mpc("0.08", "0.11")))
        decomposition = tuple(error for a, b in phases
                              for error in decomposition_residuals(a, b))
        ratios = tuple(identity_760_ratio(t, a, b) for a, b in phases
                       for t in (0, mp.mpf("0.17")))
        identity_error = max(abs(value - ratios[0]) / max(1, abs(ratios[0]))
                             for value in ratios)
        if max(decomposition) > mp.mpf("1e-12"):
            raise RuntimeError("Bobenko theta decomposition check failed")
        if max((*errors, imaginary, analytic_drift, identity_error,
                abs(initial_integrals[0] - 1))) > mp.mpf("1e-10"):
            raise RuntimeError("Bobenko rounded-data validation failed")
        if (numeric_drift > mp.mpf("1e-12")
                or refinement > mp.mpf("1e-12") or not 8 < ratio < 32):
            raise RuntimeError("Bobenko RK4 refinement check failed")
        return DocumentationComparison(
            times, fine, analytic_times, analytic, errors, imaginary,
            analytic_drift, numeric_drift, decomposition, identity_error,
            refinement, ratio,
        )


def make_figure(result):
    """Plot angular momentum and gravity; import Matplotlib only on demand."""
    import matplotlib.pyplot as plt
    from examples._plotting import comparison_styles

    figure, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    for ax, indices, label, symbol in (
            (axes[0], range(3), "Angular momentum", r"\ell"),
            (axes[1], range(3, 6), "Gravity direction", "g")):
        for style_index, j in enumerate(indices):
            line_style, marker_style = comparison_styles(style_index)
            name = rf"${symbol}_{j % 3 + 1}$"
            ax.plot([float(t) for t in result.rk4_times],
                            [float(mp.re(s[j])) for s in result.rk4_states],
                            **line_style, label=name + " (RK4)")
            ax.plot([float(t) for t in result.analytic_times],
                    [float(mp.re(s[j])) for s in result.analytic_states],
                    **marker_style,
                    label=name + " (theta)")
        ax.set_ylabel(label)
        ax.legend(fontsize=8, ncol=3, loc="lower center",
                  bbox_to_anchor=(0.5, 1.02))
    axes[1].set_xlabel("Time t")
    figure.suptitle("Bobenko–Reyman–Semenov-Tian-Shansky genus-three top")
    figure.tight_layout()
    return figure


def parse_arguments():
    """Parse numerical-comparison options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=30)
    parser.add_argument("--stop", default="2")
    parser.add_argument("--steps", type=int, default=400)
    parser.add_argument("--samples", type=int, default=41)
    parser.add_argument("--tol", default="1e-6")
    return parser.parse_args()


def main():
    """Print formula checks and the analytic-versus-RK4 comparison."""
    arguments = parse_arguments()
    mp.dps = arguments.dps
    phases = ((mp.mpc("0.13", "0.07"), mp.mpc("-0.09", "0.04")),
              (mp.mpc("0.21", "-0.03"), mp.mpc("0.08", "0.11")))
    for p1, p2 in phases:
        print("phase:", p1, p2)
        print(" decomposition residuals:", decomposition_residuals(p1, p2))
        print(" equation (7.60) ratios:",
              identity_760_ratio(0, p1, p2),
              identity_760_ratio(mp.mpf("0.17"), p1, p2))
        for time in (0, mp.mpf("0.17")):
            ell, gravity = solution(time, p1, p2)
            print(" t =", time)
            print(" ell =", tuple(ell))
            print(" g =", tuple(gravity))
            print(" invariants =", invariants(ell, gravity))
    errors, imaginary, invariant_drift = rk4_comparison(
        mp.mpf(arguments.stop), arguments.steps, arguments.samples)
    print("maximum component errors versus RK4:",
          tuple(mp.nstr(error, 8) for error in errors))
    print("maximum imaginary component:", mp.nstr(imaginary, 8))
    print("maximum analytic invariant drift:", mp.nstr(invariant_drift, 8))
    if max(errors + (imaginary, invariant_drift)) > mp.mpf(arguments.tol):
        raise SystemExit("genus-three formula failed its numerical tolerance")


if __name__ == "__main__":
    main()
