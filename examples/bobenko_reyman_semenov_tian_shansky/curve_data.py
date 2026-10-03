#!/usr/bin/env python3
"""Compute marked Abel images and flow residues for the Kowalewski curve.

This companion to ``kowalewski_genus_three`` constructs the geometric data
for H=3/2, I1=1/5, I2=27/5, following BRS (1989), equations (5.2)-(5.4)
and (7.9). All periods and integrals use Genera's geometric marking.
No external period matrix or marked-point values enter the calculation.

Run from the repository root:

    python -m examples.bobenko_reyman_semenov_tian_shansky.curve_data

The existing theta solution uses a different, involution-adapted marking.
Conversion to that marking, the third-kind normalization scalar, and
recovery of physical initial conditions are outside this companion.
"""

import argparse
from dataclasses import dataclass

from genera import Curve
from mpmath import mp


def curve_polynomial():
    """Return F(x,y) for the fixed genus-three Kowalewski model."""
    return {
        (2, 4): 1, (1, 2): -2, (2, 2): 6, (3, 2): -4,
        (0, 0): 1, (1, 0): -mp.mpf(26) / 5, (2, 0): mp.mpf(27) / 5,
    }


def holomorphic_differentials():
    """Return the three forms x/F_y, xy/F_y, (xy^2-1)/F_y times dx."""
    def fy(x, y):
        return 4 * x**2 * y**3 - 4 * (x - 3*x**2 + 2*x**3) * y

    return (lambda x, y: x / fy(x, y),
            lambda x, y: x*y / fy(x, y),
            lambda x, y: (x*y**2 - 1) / fy(x, y))


def make_curve():
    """Bind the holomorphic basis to the ambient polynomial."""
    return Curve(polynomial=curve_polynomial(),
                 differentials_kind_1=holomorphic_differentials())


def zero_chart(curve):
    """Resolve the two places above x=0 using x=t^2, y=1/(t(1+tw))."""
    fifth = mp.mpf(1) / 5
    polynomial = {
        (0, 0): 4*fifth, (2, 0): 7*fifth,
        (1, 1): -44*fifth, (3, 1): 68*fifth,
        (0, 2): 4, (2, 2): -126*fifth, (4, 2): 142*fifth,
        (1, 3): 4, (3, 3): -104*fifth, (5, 3): 108*fifth,
        (2, 4): 1, (4, 4): -26*fifth, (6, 4): 27*fifth,
    }
    return curve.chart(
        polynomial, lambda t, w: (t**2, 1/(t*(1+t*w)), 2*t))


def velocity_residues(curve, chart, place, radius, *, steps=24):
    """Compute Res(y*du/2) by integrating a closed loop in a local chart."""
    circle = tuple(radius * mp.exp(2*mp.pi*mp.j*index/steps)
                   for index in range(steps))
    circle += (circle[0],)
    forms = tuple(lambda x, y, du=du: y*du(x, y)/2
                  for du in curve.differentials_kind_1)
    # In this chart x=t^-2, y=t^-1*w, so w=t*y at the junction.
    integral = curve.chart_integral(chart, forms, circle, radius*place.y)
    return mp.matrix([value/(2*mp.pi*mp.j) for value in integral.values])


def residue_example():
    """Return the raw flow residues without the global period calculation."""
    curve = make_curve()
    chart = curve.monomial_chart(-2, -1)
    radius = mp.mpf("0.05")
    place = curve.chart_place(chart, 2, radius)
    return velocity_residues(curve, chart, place, radius)


@dataclass
class CurveData:
    """Computed data in a single geometric marking and differential basis."""

    tau: object
    marked_shifts: dict
    velocity: object
    raw_residues: object
    residue_error: object
    radius_error: object
    divisor_error: object


def compute_curve_data(*, cutoff="0.05"):
    """Compute periods, four marked places, their Abel shifts and velocity."""
    radius = mp.mpf(cutoff)
    if not 0 < radius <= mp.mpf("0.1"):
        raise ValueError("cutoff must be positive and at most 0.1")
    curve = make_curve()
    periods = curve.periods_kind_1()
    if periods.genus != 3 or not curve.validate(periods).passed:
        raise RuntimeError("Kowalewski first-kind period validation failed")
    local_zero = zero_chart(curve)
    zero_places = tuple(curve.chart_place(local_zero, root, radius)
                        for root in curve.chart_fibre(local_zero, 0))
    growing_chart = curve.monomial_chart(-2, -1)
    growing = curve.chart_place(growing_chart, 2, radius)
    vanishing_chart = curve.monomial_chart(-2, 1)
    vanishing = curve.chart_place(
        vanishing_chart, mp.sqrt(mp.mpf(27)/20), radius)
    places = {
        "zero_0": zero_places[0], "zero_1": zero_places[1],
        "infinity_growing": growing, "infinity_vanishing": vanishing,
    }
    # The common origin is the growing infinity place. The same global
    # cycle marking and form order apply to every integral below.
    abel = {name: curve.abel_map_kind_1(place).value
            for name, place in places.items()}
    raw_shifts = {name: value - abel["infinity_growing"]
                  for name, value in abel.items()}
    inverse_a = (2*periods.omega)**-1
    shifts = {name: inverse_a*value for name, value in raw_shifts.items()}
    residues = velocity_residues(curve, growing_chart, growing, radius)
    smaller = radius/2
    inner_place = curve.chart_place(growing_chart, 2, smaller)
    inner_residues = velocity_residues(
        curve, growing_chart, inner_place, smaller)
    # div(x) = 2*zero_0 + 2*zero_1 - 2*infinity_growing
    #          - 2*infinity_vanishing. Abel's theorem gives a period vector.
    relation = 2*(raw_shifts["zero_0"] + raw_shifts["zero_1"]
                  - raw_shifts["infinity_vanishing"])
    divisor_error = mp.norm(curve.lattice_reduce(relation, periods).value)
    residue_error = mp.norm(residues - mp.matrix([0, 0, -mp.mpf(1)/2]))
    radius_error = mp.norm(residues-inner_residues)
    if max(residue_error, radius_error, divisor_error) > mp.mpf("1e-10"):
        raise RuntimeError("marked-divisor or flow-residue check failed")
    return CurveData(periods.tau, shifts, inverse_a*residues, residues,
                     residue_error, radius_error, divisor_error)


def main():
    """Print the computed curve data and mathematical validation residuals."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=20)
    parser.add_argument("--cutoff", default="0.05")
    args = parser.parse_args()
    mp.dps = args.dps
    data = compute_curve_data(cutoff=args.cutoff)
    print("normalized period matrix (geometric marking):")
    print(mp.nstr(data.tau, 12))
    for name, shift in data.marked_shifts.items():
        print(name, "normalized Abel shift:", mp.nstr(shift, 12))
    print("normalized velocity:", mp.nstr(data.velocity, 12))
    print("raw residues:", mp.nstr(data.raw_residues, 12))
    print("residue residual:", mp.nstr(data.residue_error, 6))
    print("radius-change residual:", mp.nstr(data.radius_error, 6))
    print("principal-divisor residual:", mp.nstr(data.divisor_error, 6))


if __name__ == "__main__":
    main()
