#!/usr/bin/env python3
"""Reproduce Bernatska's genus-three trigonal period and P-function data.

This is Example 3 and Examples 3a--3b of J. Bernatska, "Computation of
P-Functions on Plane Algebraic Curves" (arXiv:2407.05632v4).

Run from the Genera repository root with

    .venv/bin/python -m examples.bernatska.bernatska_trigonal_demo
"""

from genera import Curve, kleinian_p
from mpmath import mp


def trigonal_curve():
    """Return Bernatska's (3,4)-curve from equation (73)."""
    return {
        (0, 3): -1,
        (4, 0): 1,
        (3, 0): 3,
        (2, 0): 7,
        (1, 0): 16,
        (0, 0): 9,
        (2, 1): 4,
        (1, 1): 5,
        (0, 1): 11,
    }


def differentials():
    """Return Bernatska's three first- and three second-kind forms."""
    def denominator(x, y):
        return -3 * y**2 + 4 * x**2 + 5 * x + 11

    first_kind = (
        lambda x, y: y / denominator(x, y),
        lambda x, y: x / denominator(x, y),
        lambda x, y: 1 / denominator(x, y),
    )

    def r5(x, y):
        return (5 * x**2 * y + 9 * x * y + mp.mpf(32) / 3 * x**2
                + 7 * y + mp.mpf(40) / 3 * x)

    second_kind = (
        lambda x, y: x**2 / denominator(x, y),
        lambda x, y: 2 * x * y / denominator(x, y),
        lambda x, y: r5(x, y) / denominator(x, y),
    )
    return first_kind, second_kind


def paper_periods():
    """Return the six-digit period matrices printed in Example 3."""
    c = mp.mpc
    omega = mp.matrix([
        [-c(0, "0.646716"), c(0, "1.367235"), -c(0, "1.406214")],
        [c(0, "0.557691"), c(0, "0.662524"), c(0, "0.237700")],
        [-c(0, "0.425220"), -c(0, "0.085658"), c(0, "0.761823")],
    ])
    omega_prime = mp.matrix([
        [c("1.114221", "0.360259"), c("-0.838244", "0.360259"),
         c("0.830310", "-0.703107")],
        [c("-0.888801", "0.610108"), c("-0.725076", "0.610108"),
         c("-0.483530", "0.118850")],
        [c("0.212490", "-0.255439"), c("0.017209", "-0.255439"),
         c("-0.244951", "0.380911")],
    ])
    eta = mp.matrix([
        [-c(0, "0.541959"), -c(0, "0.385425"), -c(0, "0.722057")],
        [c(0, "1.52536"), -c(0, "0.88414"), c(0, "0.484784")],
        [c(0, "0.975636"), -c(0, "1.01088"), -c(0, "2.65892")],
    ])
    eta_prime = mp.matrix([
        [c("-1.357307", "-0.463692"), c("2.945439", "-0.463692"),
         c("-0.354124", "-0.361028")],
        [c("2.292766", "0.320609"), c("5.356432", "0.320609"),
         c("2.131611", "0.242392")],
        [c("-5.584050", "-0.017623"), c("4.038588", "-0.017623"),
         c("6.689080", "-1.329459")],
    ])
    return join_periods(omega, omega_prime), join_periods(eta, eta_prime)


def join_periods(left, right):
    result = mp.matrix(left.rows, left.cols + right.cols)
    for row in range(left.rows):
        for column in range(left.cols):
            result[row, column] = left[row, column]
            result[row, column + left.cols] = right[row, column]
    return result


def realify(periods):
    return mp.matrix([
        [mp.re(periods[row, column]) for column in range(6)]
        for row in range(3)
    ] + [
        [mp.im(periods[row, column]) for column in range(6)]
        for row in range(3)
    ])


def nearest_integer_matrix(matrix):
    return mp.matrix([
        [int(mp.nint(matrix[row, column]))
         for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ])


def maximum_abs(matrix):
    return max((abs(value) for value in matrix), default=mp.zero)


def standard_intersection(genus):
    result = mp.zeros(2 * genus)
    for index in range(genus):
        result[index, genus + index] = 1
        result[genus + index, index] = -1
    return result


def paper_example_3a_values(omega, tau, kappa):
    """Evaluate the eight P-functions printed in equation (86)."""
    c = mp.mpc
    u = mp.matrix([
        c("-0.270333", "-1.612257"),
        c("-1.116879", "0.562199"),
        c("0.258194", "0.268653"),
    ])
    characteristic = (
        (mp.mpf("0.5"), 0, mp.mpf("0.5")),
        (0, mp.mpf("0.5"), mp.mpf("0.5")),
    )
    indices = (
        (0, 0), (0, 1), (0, 2), (1, 1), (1, 2),
        (0, 0, 0), (0, 0, 1), (0, 0, 2),
    )
    values = kleinian_p(
        u, omega / 2, tau, -kappa, indices, characteristic)
    expected = (
        c("0.059654", "1.020925"), c("-0.793416", "0.889005"),
        c("0.885372", "-3.089764"), c("-0.269700", "1.472739"),
        c("-3.501466", "10.538856"), c("-2.156576", "3.543516"),
        c("-3.595029", "2.840859"), c("5.656516", "-0.559812"),
    )
    return values, expected


def paper_example_3b_values(omega, tau, kappa):
    """Evaluate the eight P-functions printed in equation (91)."""
    c = mp.mpc
    u = mp.matrix([
        c("-0.421105", "-2.303962"),
        c("-1.319230", "-1.997581"),
        c("-0.176345", "0.125109"),
    ])
    characteristic = (
        (mp.mpf("0.5"), 0, mp.mpf("0.5")),
        (0, mp.mpf("0.5"), mp.mpf("0.5")),
    )
    indices = (
        (0, 0), (0, 1), (0, 2), (1, 1), (1, 2),
        (0, 0, 0), (0, 0, 1), (0, 0, 2),
    )
    values = kleinian_p(
        u, omega / 2, tau, -kappa, indices, characteristic)
    expected = (
        c("-0.497171", "-1.306218"), c("0.485105", "2.618402"),
        c("2.083016", "-2.086324"), c("-2.356414", "10.869587"),
        c("15.590831", "2.902800"), c("1.678988", "8.731706"),
        c("-4.377331", "-0.119524"), c("2.198126", "13.211222"),
    )
    return values, expected


def main():
    mp.dps = 20

    curve = Curve(mp, trigonal_curve())
    first_kind, second_kind = differentials()
    first_data = curve.periods_kind_1(first_kind)
    second_data = curve.periods_kind_2(
        differentials=first_kind,
        second_differentials=second_kind,
    )
    points = curve.branch_locus.branch_values
    monodromy = curve.monodromy
    permutations = tuple(reversed(monodromy.permutations)) + (
        monodromy.infinity_permutation,)
    homology = curve.homology
    first = join_periods(
        2 * first_data.omega, 2 * first_data.omega_prime)
    second = join_periods(
        -2 * second_data.eta, -2 * second_data.eta_prime)
    paper_first, paper_second = paper_periods()

    numerical_change = realify(first)**-1 * realify(paper_first)
    cycle_change = nearest_integer_matrix(numerical_change)
    first_in_paper_basis = first * cycle_change
    second_in_paper_basis = second * cycle_change
    omega = first_in_paper_basis[:, :3]
    omega_prime = first_in_paper_basis[:, 3:]
    eta = second_in_paper_basis[:, :3]
    eta_prime = second_in_paper_basis[:, 3:]
    tau = omega**-1 * omega_prime
    kappa = eta * omega**-1
    tau_symmetry = mp.norm(tau - tau.T)
    kappa_symmetry = mp.norm(kappa - kappa.T)
    # Quadrature leaves antisymmetric roundoff well above mp.eps.  Projecting
    # onto the identities already checked below gives the exact matrix shape
    # expected by the theta evaluator without changing the reported test.
    tau = (tau + tau.T) / 2
    kappa = (kappa + kappa.T) / 2
    imaginary_tau = mp.matrix([
        [mp.im(tau[row, column]) for column in range(3)]
        for row in range(3)
    ])
    eigenvalues = mp.eigsy(imaginary_tau, eigvals_only=True)
    p_values, p_expected = paper_example_3a_values(omega, tau, kappa)
    p_values_3b, p_expected_3b = paper_example_3b_values(
        omega, tau, kappa)

    branch_reference = (
        mp.mpc("-4.58931"), mp.mpc("-1.17922", "-0.934455"),
        mp.mpc("-1.17922", "0.934455"),
        mp.mpc("-0.431732", "-2.20256"),
        mp.mpc("-0.431732", "2.20256"),
        mp.mpc("0.499118", "-1.57527"),
        mp.mpc("0.499118", "1.57527"), mp.mpc("0.812986"),
    )
    branch_residual = max(
        min(abs(point - reference) for point in points)
        for reference in branch_reference)
    integer_deviation = maximum_abs(numerical_change - cycle_change)
    intersection = standard_intersection(3)
    symplectic_residual = maximum_abs(
        cycle_change.T * intersection * cycle_change - intersection)
    first_residual = maximum_abs(first_in_paper_basis - paper_first)
    second_residual = maximum_abs(second_in_paper_basis - paper_second)
    p_residual = max(abs(value - expected)
                     for value, expected in zip(p_values, p_expected))
    p_residual_3b = max(
        abs(value - expected)
        for value, expected in zip(p_values_3b, p_expected_3b))

    print("finite branch permutations:", permutations[:-1])
    print("infinity permutation:", permutations[-1])
    print("genus / intersection rank:", homology.genus,
          homology.intersection_rank)
    print("cycle change to Bernatska's basis:")
    print(cycle_change)
    print("maximum distance to integer:", mp.nstr(integer_deviation, 8))
    print("symplectic cycle residual:", mp.nstr(symplectic_residual, 8))
    print("branch-point residual:", mp.nstr(branch_residual, 8))
    print("first-kind period residual:", mp.nstr(first_residual, 8))
    print("second-kind period residual:", mp.nstr(second_residual, 8))
    print("tau symmetry residual:", mp.nstr(tau_symmetry, 8))
    print("Im(tau) eigenvalues:",
          [mp.nstr(value, 8) for value in eigenvalues])
    print("kappa symmetry residual:", mp.nstr(kappa_symmetry, 8))
    print("Example 3a maximum P-function residual:",
          mp.nstr(p_residual, 8))
    print("Example 3b maximum P-function residual:",
          mp.nstr(p_residual_3b, 8))

    if (branch_residual > mp.mpf("8e-6")
            or integer_deviation > mp.mpf("2e-6")
            or abs(mp.det(cycle_change)) != 1
            or symplectic_residual != 0
            or first_residual > mp.mpf("2e-6")
            or second_residual > mp.mpf("2e-5")
            or tau_symmetry > mp.mpf("1e-12")
            or min(eigenvalues) <= 0
            or kappa_symmetry > mp.mpf("1e-11")
            or p_residual > mp.mpf("2e-3")
            or p_residual_3b > mp.mpf("2e-3")):
        raise RuntimeError("Bernatska trigonal validation failed")


if __name__ == "__main__":
    main()
