#!/usr/bin/env python3
"""Numerically verify Onishi's determinant formulas in arbitrary genus.

The example implements Theorems 7.2 and 8.3 of Y. Onishi, "Determinant
Expressions for Hyperelliptic Functions" (arXiv:math/0105189).  Theorem 7.2
equates a multipoint quotient of sigma-stratum derivatives to a determinant
of affine curve monomials.  Theorem 8.3 takes its confluent limit, equating a
sigma division function to a determinant of repeated curve derivatives.

The curve is written in the Genera convention

    Y**2 = 4 * product(x - e_i),

so Genera's holomorphic differentials x**j dx/Y are exactly Onishi's
x**j dx/(2*y), with y = Y/2.  All sigma values use Genera's canonical
``normalization="hyperelliptic"`` convention.  The constant multiplying
the sigma quotient is fixed from the leading Schur--Weierstrass polynomial
in that convention.  This is an exact local calculation, not a fitted
factor.  The script also prints the sign in Onishi's displayed statement
and the resulting convention conversion.

For n < g the numerator is sigma_natural^n; for n >= g it is sigma.
The special derivative is

    natural^n = {i : n + 1 <= i <= g and i == n + 1 (mod 2)},

with Onishi's one-based indices converted to zero-based sigma-jet keys.
The determinant columns are generated in increasing pole order at infinity:

    1, x, ..., x**g, y, x**(g+1), x*y, x**(g+2), ... .

For the Kiepert determinant, the initial constant column is removed and rows
contain derivative orders 1 through n-1.  Along the curve, the script applies

    d/du_j = Y/x**(j-1) * d/dx

exactly to expressions A(x) + B(x)*Y using Y**2=P(x).  Its sign is again
fixed independently from the local Schur--Weierstrass model.

Run from the Genera repository root, for example

    .venv/bin/python -m examples.onishi.onishi_determinant_demo
    .venv/bin/python -m examples.onishi.onishi_determinant_demo \
        --genus 4 --n 2 4 6 --kiepert-n 4 6 --coordinate 1 4 --dps 25

Theta-series cost grows exponentially with genus; the formulas and column
construction are genus-independent, but low genera are the practical
numerical examples.
"""

import argparse
from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations
from math import factorial

from genera import Curve, kleinian_sigma, kleinian_sigma_jet
from mpmath import mp


def multiply_by_linear(coefficients, root):
    """Multiply an ascending coefficient vector by x - root."""
    result = [mp.zero] * (len(coefficients) + 1)
    for degree, coefficient in enumerate(coefficients):
        result[degree] -= root * coefficient
        result[degree + 1] += coefficient
    return result


def curve_coefficients(genus):
    """Return 4*product(x-e) for the roots -g, ..., g."""
    coefficients = [mp.mpf(4)]
    for root in range(-genus, genus + 1):
        coefficients = multiply_by_linear(coefficients, mp.mpf(root))
    return coefficients


def polynomial_value(coefficients, x):
    """Evaluate an ascending coefficient vector."""
    return mp.fsum(coefficient * x ** degree
                   for degree, coefficient in enumerate(coefficients))


def natural_index(genus, n):
    """Return Onishi's natural^n as a zero-based derivative-count tuple."""
    return tuple(
        int(index >= n and (index - n) % 2 == 0)
        for index in range(genus)
    )


def onishi_monomials(x, y, genus, count):
    """Return the first ``count`` affine monomials in pole order."""
    result = []
    x_power = 0
    y_power = 0
    while len(result) < count:
        x_pole_order = 2 * x_power
        y_pole_order = 2 * genus + 1 + 2 * y_power
        if x_pole_order < y_pole_order:
            result.append(x ** x_power)
            x_power += 1
        else:
            result.append(y * x ** y_power)
            y_power += 1
    return result


def onishi_monomial_pairs(genus, count):
    """Represent the first monomials as A(x) + B(x)*Y, where Y=2y."""
    result = []
    x_power = 0
    y_power = 0
    while len(result) < count:
        x_pole_order = 2 * x_power
        y_pole_order = 2 * genus + 1 + 2 * y_power
        if x_pole_order < y_pole_order:
            result.append(({x_power: mp.one}, {}))
            x_power += 1
        else:
            result.append(({}, {y_power: mp.mpf("0.5")}))
            y_power += 1
    return result


def laurent_add(left, right):
    result = dict(left)
    for exponent, coefficient in right.items():
        result[exponent] = result.get(exponent, mp.zero) + coefficient
    return {exponent: coefficient for exponent, coefficient in result.items()
            if coefficient}


def laurent_scale(polynomial, scalar):
    return {exponent: scalar * coefficient
            for exponent, coefficient in polynomial.items()
            if scalar * coefficient}


def laurent_shift(polynomial, shift):
    return {exponent + shift: coefficient
            for exponent, coefficient in polynomial.items()}


def laurent_derivative(polynomial):
    return {exponent - 1: exponent * coefficient
            for exponent, coefficient in polynomial.items() if exponent}


def laurent_multiply(left, right):
    result = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = left_exponent + right_exponent
            result[exponent] = (
                result.get(exponent, mp.zero)
                + left_coefficient * right_coefficient
            )
    return {exponent: coefficient for exponent, coefficient in result.items()
            if coefficient}


def curve_derivative(pair, coefficients, coordinate):
    """Apply d/du_coordinate along Y**2=P(x) to A(x)+B(x)Y."""
    a_polynomial, b_polynomial = pair
    curve = {degree: coefficient
             for degree, coefficient in enumerate(coefficients)
             if coefficient}
    curve_prime = laurent_derivative(curve)
    new_a = laurent_add(
        laurent_multiply(laurent_derivative(b_polynomial), curve),
        laurent_scale(
            laurent_multiply(b_polynomial, curve_prime), mp.mpf("0.5")),
    )
    new_b = laurent_derivative(a_polynomial)
    shift = -(coordinate - 1)
    return laurent_shift(new_a, shift), laurent_shift(new_b, shift)


def evaluate_monomial_pair(pair, x, y):
    """Evaluate A(x)+B(x)Y at Y=2y."""
    a_polynomial, b_polynomial = pair
    a_value = mp.fsum(coefficient * x ** exponent
                      for exponent, coefficient in a_polynomial.items())
    b_value = mp.fsum(coefficient * x ** exponent
                      for exponent, coefficient in b_polynomial.items())
    return a_value + b_value * (2 * y)


# The following small exact-polynomial helpers construct the leading
# Schur--Weierstrass term of Genera's hyperelliptic sigma.  A polynomial is
# a mapping from exponent tuples to rational coefficients.


def polynomial_add(left, right):
    result = dict(left)
    for exponent, coefficient in right.items():
        result[exponent] = result.get(exponent, Fraction(0)) + coefficient
    return {exponent: coefficient for exponent, coefficient in result.items()
            if coefficient}


def polynomial_scale(polynomial, scalar):
    return {exponent: scalar * coefficient
            for exponent, coefficient in polynomial.items()
            if scalar * coefficient}


def polynomial_multiply(left, right):
    result = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(a + b for a, b in
                             zip(left_exponent, right_exponent))
            result[exponent] = (
                result.get(exponent, Fraction(0))
                + left_coefficient * right_coefficient
            )
    return {exponent: coefficient for exponent, coefficient in result.items()
            if coefficient}


def polynomial_derivative(polynomial, derivative):
    result = polynomial
    for coordinate, order in enumerate(derivative):
        for unused in range(order):
            differentiated = {}
            for exponent, coefficient in result.items():
                if exponent[coordinate]:
                    lowered = list(exponent)
                    lowered[coordinate] -= 1
                    differentiated[tuple(lowered)] = (
                        coefficient * exponent[coordinate]
                    )
            result = differentiated
    return result


def polynomial_evaluate(polynomial, argument):
    return sum(
        coefficient * rational_product(
            value ** exponent
            for value, exponent in zip(argument, exponents)
        )
        for exponents, coefficient in polynomial.items()
    )


def rational_product(values):
    result = Fraction(1)
    for value in values:
        result *= value
    return result


def permutation_sign(permutation):
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def schur_weierstrass_polynomial(genus):
    """Return the leading sigma polynomial in Genera normalization.

    Onishi Section 1 writes the Schur--Weierstrass polynomial as
    det(U_(g-2*i+j+1)).  The generating series in Abelian coordinates is
    exp(sum(u_j*z**(2*(g-j)+1))).  The final scale is chosen exactly as in
    ``kleinian_sigma(..., normalization="hyperelliptic")``: the coefficient
    of u_m**m is the reversing-permutation sign, m=floor((g+1)/2).
    """
    zero = (0,) * genus
    weights = [2 * (genus - coordinate) + 1
               for coordinate in range(1, genus + 1)]
    complete = [{zero: Fraction(1)}]
    for degree in range(1, 2 * genus):
        coefficient = {}
        for coordinate, weight in enumerate(weights):
            if weight > degree:
                continue
            multiplied = {}
            for exponent, value in complete[degree - weight].items():
                raised = list(exponent)
                raised[coordinate] += 1
                multiplied[tuple(raised)] = value
            coefficient = polynomial_add(
                coefficient,
                polynomial_scale(multiplied, Fraction(weight, degree)),
            )
        complete.append(coefficient)

    result = {}
    for permutation in permutations(range(genus)):
        term = {zero: Fraction(1)}
        for row, column in enumerate(permutation):
            degree = genus - 2 * (row + 1) + (column + 1) + 1
            if degree < 0:
                term = {}
                break
            term = polynomial_multiply(term, complete[degree])
        result = polynomial_add(
            result, polynomial_scale(term, permutation_sign(permutation)))

    normalization_degree = (genus + 1) // 2
    normalization_coordinate = normalization_degree - 1
    normalization_exponent = [0] * genus
    normalization_exponent[normalization_coordinate] = normalization_degree
    coefficient = result[tuple(normalization_exponent)]
    target = Fraction(
        -1 if normalization_degree * (normalization_degree - 1) // 2 % 2
        else 1
    )
    return polynomial_scale(result, target / coefficient)


def rational_determinant(rows):
    """Evaluate a small exact determinant by Gaussian elimination."""
    matrix = [list(row) for row in rows]
    determinant = Fraction(1)
    for column in range(len(matrix)):
        pivot = next(row for row in range(column, len(matrix))
                     if matrix[row][column])
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            determinant = -determinant
        pivot_value = matrix[column][column]
        determinant *= pivot_value
        for entry in range(column, len(matrix)):
            matrix[column][entry] /= pivot_value
        for row in range(column + 1, len(matrix)):
            factor = matrix[row][column]
            for entry in range(column, len(matrix)):
                matrix[row][entry] -= factor * matrix[column][entry]
    return determinant


def local_abel_point(genus, parameter):
    """Return the leading one-point Abel expansion with u_g=parameter."""
    return tuple(
        parameter ** (2 * (genus - coordinate) + 1)
        / Fraction(2 * (genus - coordinate) + 1)
        for coordinate in range(1, genus + 1)
    )


def local_monomials(genus, parameter, count):
    """Return leading affine monomials for x=t^-2, y=-t^-(2g+1)."""
    x = parameter ** -2
    y = -parameter ** (-(2 * genus + 1))
    return onishi_monomials(x, y, genus, count)


def genera_prefactor(genus, n, schur):
    """Derive the formula sign from Genera's leading sigma polynomial."""
    parameters = [Fraction(index + 1) for index in range(n)]
    arguments = [local_abel_point(genus, parameter)
                 for parameter in parameters]
    total = tuple(sum(argument[index] for argument in arguments)
                  for index in range(genus))
    sharp = polynomial_derivative(schur, natural_index(genus, 1))
    flat = polynomial_derivative(schur, natural_index(genus, 2))
    numerator_index = (natural_index(genus, n) if n < genus
                       else (0,) * genus)
    numerator = polynomial_evaluate(
        polynomial_derivative(schur, numerator_index), total)
    for left in range(n):
        for right in range(left + 1, n):
            difference = tuple(
                arguments[left][index] - arguments[right][index]
                for index in range(genus)
            )
            numerator *= polynomial_evaluate(flat, difference)
    denominator = rational_product(
        polynomial_evaluate(sharp, argument) ** n
        for argument in arguments
    )
    quotient = numerator / denominator
    determinant = rational_determinant([
        local_monomials(genus, parameter, n)
        for parameter in parameters
    ])
    conversion = determinant / quotient
    if conversion not in (Fraction(-1), Fraction(1)):
        raise ArithmeticError(
            f"unexpected local normalization factor {conversion}")
    return int(conversion)


def local_monomial_terms(genus, count):
    """Return (coefficient, t-exponent) for local affine monomials."""
    result = []
    x_power = 0
    y_power = 0
    while len(result) < count:
        x_pole_order = 2 * x_power
        y_pole_order = 2 * genus + 1 + 2 * y_power
        if x_pole_order < y_pole_order:
            result.append((Fraction(1), -2 * x_power))
            x_power += 1
        else:
            result.append((Fraction(-1), -y_pole_order))
            y_power += 1
    return result


def local_curve_derivative(term, genus, coordinate):
    """Apply d/du_coordinate to coefficient*t**exponent locally."""
    coefficient, exponent = term
    weight = 2 * (genus - coordinate) + 1
    return coefficient * exponent, exponent - weight


def genera_kiepert_prefactor(genus, n, coordinate, schur):
    """Derive the Kiepert sign from the exact local sigma model."""
    argument = local_abel_point(genus, Fraction(1))
    scaled_argument = tuple(n * value for value in argument)
    sharp = polynomial_derivative(schur, natural_index(genus, 1))
    numerator_index = (natural_index(genus, n) if n < genus
                       else (0,) * genus)
    numerator = polynomial_evaluate(
        polynomial_derivative(schur, numerator_index), scaled_argument)
    denominator = polynomial_evaluate(sharp, argument) ** (n ** 2)
    psi = numerator / denominator

    columns = local_monomial_terms(genus, n)[1:]
    rows = []
    differentiated = list(columns)
    for unused_order in range(1, n):
        differentiated = [
            local_curve_derivative(term, genus, coordinate)
            for term in differentiated
        ]
        # The local parameter is t=1, so only the coefficient remains.
        rows.append([coefficient for coefficient, unused_exponent
                     in differentiated])
    determinant = rational_determinant(rows)
    factorial_product = rational_product(
        Fraction(factorial(order)) for order in range(1, n)
    )
    conversion = determinant / (factorial_product * psi)
    if conversion not in (Fraction(-1), Fraction(1)):
        raise ArithmeticError(
            f"unexpected local Kiepert normalization factor {conversion}")
    return int(conversion)


def printed_onishi_prefactor(genus, n):
    """Return c_n as displayed in Theorem 7.2."""
    if n < genus:
        exponent = (
            genus + 1 + (n - 1) * (n - 2) * (n - 3) // 2
        )
        return -1 if exponent % 2 else 1
    table = {
        1: {1: 1, 2: 1, 3: -1, 0: -1},
        2: {1: -1, 2: -1, 3: -1, 0: -1},
        3: {1: -1, 2: 1, 3: 1, 0: -1},
        0: {1: -1, 2: 1, 3: -1, 0: 1},
    }
    return table[genus % 4][n % 4]


def printed_kiepert_prefactor(genus, n):
    """Return c'_n as displayed in Theorem 8.3."""
    genus_class = genus % 4
    if genus_class == 1:
        exponent = n - 1
        return -1 if exponent % 2 else 1
    if genus_class == 2:
        exponent = n * (n - 1) // 2
        return 1 if exponent % 2 else -1
    if genus_class == 3:
        return -1
    exponent = n * (n + 1) // 2
    return -1 if exponent % 2 else 1


def sigma_derivative(argument, derivative, data):
    """Evaluate one canonically normalized sigma derivative."""
    omega, tau, kappa, characteristic = data
    jet = kleinian_sigma_jet(
        argument, omega, tau, kappa, sum(derivative), characteristic,
        normalization="hyperelliptic",
    )
    return jet[derivative]


def evaluate_formula(points, images, data, genus, prefactor):
    """Return the two sides of Theorem 7.2 in Genera conventions."""
    n = len(points)
    total = [mp.fsum(image[index] for image in images)
             for index in range(genus)]
    if n < genus:
        numerator = sigma_derivative(
            total, natural_index(genus, n), data)
    else:
        omega, tau, kappa, characteristic = data
        numerator = kleinian_sigma(
            total, omega, tau, kappa, characteristic,
            normalization="hyperelliptic",
        )

    flat_index = natural_index(genus, 2)
    for left in range(n):
        for right in range(left + 1, n):
            difference = [
                images[left][index] - images[right][index]
                for index in range(genus)
            ]
            numerator *= sigma_derivative(difference, flat_index, data)

    sharp_index = natural_index(genus, 1)
    denominator = mp.fprod(
        sigma_derivative(image, sharp_index, data) ** n
        for image in images
    )
    left_side = prefactor * numerator / denominator
    right_side = mp.det(mp.matrix([
        onishi_monomials(x, y, genus, n) for x, y in points
    ]))
    return left_side, right_side


def evaluate_kiepert_formula(point, image, coefficients, data, genus, n,
                             coordinate, prefactor):
    """Return the two sides of Onishi's Kiepert determinant formula."""
    x, y = point
    scaled_image = [n * value for value in image]
    omega, tau, kappa, characteristic = data
    if n < genus:
        numerator = sigma_derivative(
            scaled_image, natural_index(genus, n), data)
    else:
        numerator = kleinian_sigma(
            scaled_image, omega, tau, kappa, characteristic,
            normalization="hyperelliptic",
        )
    sharp = sigma_derivative(image, natural_index(genus, 1), data)
    psi = numerator / sharp ** (n ** 2)
    factorial_product = mp.fprod(factorial(order)
                                 for order in range(1, n))
    left_side = prefactor * factorial_product * psi

    columns = onishi_monomial_pairs(genus, n)[1:]
    rows = []
    differentiated = list(columns)
    for unused_order in range(1, n):
        differentiated = [
            curve_derivative(pair, coefficients, coordinate)
            for pair in differentiated
        ]
        rows.append([
            evaluate_monomial_pair(pair, x, y)
            for pair in differentiated
        ])
    determinant = mp.det(mp.matrix(rows))
    power = (coordinate - 1) * n * (n - 1) // 2
    right_side = x ** power * determinant
    return left_side, right_side


def default_sizes(genus):
    """Select nontrivial cases on both sides of the n=g boundary."""
    sizes = [2]
    if genus > 2:
        sizes.append(genus)
    sizes.append(genus + 2)
    return tuple(dict.fromkeys(sizes))


@dataclass
class DocumentationComparison:
    """Four-point determinant comparisons and confluent identity checks."""

    sample_x: tuple
    sigma_values: tuple
    determinant_values: tuple
    dense_x: tuple
    dense_determinants: tuple
    scaled_errors: tuple
    kiepert_errors: tuple
    prefactor: int
    kiepert_prefactors: tuple
    stratum_residual: object


def documentation_example():
    """Check genus-two four-point identities at 30 decimal digits.

    The first three curve points have x=4,5,6; the fourth varies over [7,9].
    The Kiepert limit is checked at x=4 in both Abelian coordinates.
    Signs are fixed by exact local calculations before theta evaluation.
    """
    with mp.workdps(30):
        genus = 2
        coefficients = curve_coefficients(genus)
        curve = Curve({(0, 2): 1, **{(i, 0): -c for i, c in enumerate(coefficients)}})
        if curve.genus != genus:
            raise RuntimeError("unexpected Onishi curve genus")
        first = curve.periods_kind_1()
        second = curve.periods_kind_2()
        data = (first.omega, first.tau, second.kappa,
                curve.riemann_constant().characteristic)
        schur = schur_weierstrass_polynomial(genus)
        prefactor = genera_prefactor(genus, 4, schur)
        kiepert_prefactors = tuple(genera_kiepert_prefactor(genus, 4, j, schur)
                                   for j in (1, 2))

        def point_at(x):
            return x, mp.sqrt(polynomial_value(coefficients, x)) / 2

        def image_of(point):
            x, y = point
            return curve.abel_map_kind_1((x, 2 * y)).value

        points = [point_at(mp.mpf(x)) for x in (4, 5, 6)]
        images = [image_of(point) for point in points]
        sample_x = tuple(mp.mpf(7) + mp.mpf(i) / 8 for i in range(17))
        sigma_values, determinant_values, errors = [], [], []
        all_images = list(images)
        for x in sample_x:
            point = point_at(x)
            image = image_of(point)
            all_images.append(image)
            left, right = evaluate_formula(points + [point], images + [image],
                                            data, genus, prefactor)
            sigma_values.append(left)
            determinant_values.append(right)
            errors.append(abs(left - right) / max(1, abs(left), abs(right)))

        kiepert_errors = []
        for j, factor in zip((1, 2), kiepert_prefactors):
            left, right = evaluate_kiepert_formula(
                points[0], images[0], coefficients, data, genus, 4, j, factor)
            kiepert_errors.append(abs(left - right) / max(1, abs(left), abs(right)))
        # A curve point lies on sigma=0. The special derivative sigma_2,
        # rather than sigma itself, supplies the denominators of the theorem.
        stratum_residual = max(abs(sigma_derivative(image, (0, 0), data))
                               for image in all_images)
        if max((*errors, *kiepert_errors, stratum_residual)) > mp.mpf("1e-24"):
            raise RuntimeError("Onishi sigma or determinant identity check failed")

        dense_x = tuple(mp.mpf(7) + mp.mpf(i) / 64 for i in range(129))
        dense_determinants = tuple(mp.det(mp.matrix([
            onishi_monomials(px, py, genus, 4)
            for px, py in points + [point_at(x)]
        ])) for x in dense_x)
        return DocumentationComparison(
            sample_x, tuple(sigma_values), tuple(determinant_values),
            dense_x, dense_determinants, tuple(errors), tuple(kiepert_errors),
            prefactor, kiepert_prefactors, stratum_residual,
        )


def make_figure(result):
    """Plot algebraic determinants and sigma evaluations; import on demand."""
    import matplotlib.pyplot as plt
    from examples._plotting import TEAL

    figure, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot([float(x) for x in result.dense_x],
            [float(mp.re(y)) for y in result.dense_determinants], "--",
            color=TEAL, label="Algebraic determinant")
    ax.plot([float(x) for x in result.sample_x],
            [float(mp.re(y)) for y in result.sigma_values], "o",
            color=TEAL, fillstyle="none", label="Sigma quotient")
    ax.set_xlabel("Fourth point's x-coordinate")
    ax.set_ylabel("Four-point determinant")
    ax.legend()
    figure.tight_layout()
    return figure


def default_kiepert_sizes(genus):
    """Select a boundary case and a determinant containing later columns."""
    return tuple(dict.fromkeys((max(2, genus), genus + 2)))


def default_coordinates(genus):
    """Exercise the first and last holomorphic coordinates."""
    return tuple(dict.fromkeys((1, genus)))


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genus", type=int, default=3,
                        help="curve genus (default: 3)")
    parser.add_argument("--n", type=int, nargs="+",
                        help="numbers of points in Theorem 7.2")
    parser.add_argument("--kiepert-n", type=int, nargs="+",
                        help="multipliers n in Theorem 8.3")
    parser.add_argument("--coordinate", type=int, nargs="+",
                        help="one-based derivative coordinates for 8.3")
    parser.add_argument("--dps", type=int, default=30,
                        help="working decimal precision (default: 30)")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    if arguments.genus < 1:
        raise ValueError("genus must be positive")
    sizes = tuple(arguments.n or default_sizes(arguments.genus))
    if any(size < 1 for size in sizes):
        raise ValueError("every point count must be positive")
    kiepert_sizes = tuple(
        arguments.kiepert_n or default_kiepert_sizes(arguments.genus))
    if any(size < max(2, arguments.genus) for size in kiepert_sizes):
        raise ValueError(
            "Kiepert determinant checks require n >= max(2, genus)")
    coordinates = tuple(
        arguments.coordinate or default_coordinates(arguments.genus))
    if any(coordinate < 1 or coordinate > arguments.genus
           for coordinate in coordinates):
        raise ValueError("Kiepert coordinates must lie between 1 and genus")

    mp.dps = arguments.dps
    genus = arguments.genus
    coefficients = curve_coefficients(genus)
    curve = Curve({(0, 2): 1, **{(i, 0): -c for i, c in enumerate(coefficients)}})
    first = curve.periods_kind_1()
    second = curve.periods_kind_2()
    data = (first.omega, first.tau, second.kappa,
            curve.riemann_constant().characteristic)
    schur = schur_weierstrass_polynomial(genus)

    maximum_size = max(sizes)
    points = []
    images = []
    for offset in range(maximum_size):
        x = mp.mpf(genus + 2 + offset)
        y = mp.sqrt(polynomial_value(coefficients, x)) / 2
        points.append((x, y))
        images.append(curve.abel_map_kind_1((x, 2 * y)).value)

    print("Onishi Theorem 7.2: arbitrary-genus determinant check")
    print(f"genus: {genus}")
    print(f"curve: Y^2 = 4*product(x-e), e = -{genus}, ..., {genus}")
    print(f"working precision: {mp.dps} digits")
    print(f"sigma_sharp jet key: {natural_index(genus, 1)}")
    print(f"sigma_flat jet key:  {natural_index(genus, 2)}")
    print()

    for n in sizes:
        prefactor = genera_prefactor(genus, n, schur)
        printed = printed_onishi_prefactor(genus, n)
        left, right = evaluate_formula(
            points[:n], images[:n], data, genus, prefactor)
        absolute_error = abs(left - right)
        comparison_scale = max(abs(left), abs(right))
        relative_error = (absolute_error / comparison_scale
                          if comparison_scale else absolute_error)
        numerator_name = (
            f"sigma_natural^{n}" if n < genus else "sigma"
        )
        print(f"n = {n} ({numerator_name} numerator)")
        print(f"  printed c_n:            {printed:+d}")
        print(f"  Genera prefactor:       {prefactor:+d}")
        print(f"  convention conversion: {prefactor * printed:+d}")
        print(f"  sigma side:             {mp.nstr(left, 16)}")
        print(f"  determinant side:       {mp.nstr(right, 16)}")
        print(f"  absolute residual:      {mp.nstr(absolute_error, 6)}")
        print(f"  relative residual:      {mp.nstr(relative_error, 6)}")
        print()

    print("Onishi Theorem 8.3: Kiepert-type confluent determinant")
    print()
    for n in kiepert_sizes:
        for coordinate in coordinates:
            prefactor = genera_kiepert_prefactor(
                genus, n, coordinate, schur)
            printed = printed_kiepert_prefactor(genus, n)
            left, right = evaluate_kiepert_formula(
                points[0], images[0], coefficients, data, genus, n,
                coordinate, prefactor)
            absolute_error = abs(left - right)
            comparison_scale = max(abs(left), abs(right))
            relative_error = (absolute_error / comparison_scale
                              if comparison_scale else absolute_error)
            print(f"n = {n}, derivative coordinate j = {coordinate}")
            print(f"  printed c'_n:           {printed:+d}")
            print(f"  Genera prefactor:       {prefactor:+d}")
            print(f"  convention conversion: {prefactor * printed:+d}")
            print(f"  sigma side:             {mp.nstr(left, 16)}")
            print(f"  determinant side:       {mp.nstr(right, 16)}")
            print(f"  absolute residual:      {mp.nstr(absolute_error, 6)}")
            print(f"  relative residual:      {mp.nstr(relative_error, 6)}")
            print()


if __name__ == "__main__":
    main()
