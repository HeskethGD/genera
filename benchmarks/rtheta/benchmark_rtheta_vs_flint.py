#!/usr/bin/env python3
"""Compare generapy rtheta with Python-FLINT's Riemann theta bindings.

Run this script with the dedicated comparison environment::

    PYTHONPATH="$PWD" venv/flint/bin/python \
        benchmarks/rtheta/benchmark_rtheta_vs_flint.py

Inputs and result construction are excluded from timed regions. The repeated
workloads retain one fixed period matrix while evaluating a sequence of
arguments, as in plotting and integrable-system applications.
"""

from __future__ import annotations

import argparse
import datetime
from fractions import Fraction
from math import factorial
from pathlib import Path
import platform
import statistics
import sys
import timeit


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_REPORT = SCRIPT_DIR / "rtheta_vs_flint.md"
sys.path.insert(0, str(REPO_ROOT/'src'))

try:
    import flint
    from flint import acb, acb_mat, fmpq
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "Run with venv/flint/bin/python; python-flint is not available"
    ) from exc

import generapy
from mpmath import mp
generapy_version = generapy.__version__

from benchmark_rtheta_vs_wolfram import (  # noqa: E402
    BASE_POINTS,
    CASES as PREVIOUS_CASES,
    case_points,
    characteristic,
    number,
)


def tridiagonal_case(name, genus, scale):
    """Return a deterministic coupled comparison case."""
    rows = []
    for i in range(genus):
        row = []
        for j in range(genus):
            if i == j:
                row.append((Fraction(i % 3, 20),
                            Fraction(scale) + Fraction(i % 4, 10)))
            elif abs(i - j) == 1:
                row.append(number("0.02", "0.15"))
            else:
                row.append(number("0"))
        rows.append(tuple(row))
    return {
        "name": name,
        "genus": genus,
        "kind": "coupled higher-genus control",
        "tau": tuple(rows),
    }


EXTRA_CASES = (
    tridiagonal_case("g4-coupled", 4, "1.4"),
    tridiagonal_case("g5-damped", 5, "3.0"),
)
CASES = PREVIOUS_CASES + EXTRA_CASES


def as_fraction(value):
    """Return an exact Fraction for a benchmark decimal."""
    return value if isinstance(value, Fraction) else Fraction(value)


def fraction_pair(value):
    """Convert a complex benchmark description to exact fractions."""
    real, imag = value
    return as_fraction(real), as_fraction(imag)


def mp_number(value):
    """Convert an exact benchmark pair to an mpmath number."""
    real, imag = fraction_pair(value)
    return mp.mpc(mp.mpf(real.numerator) / real.denominator,
                  mp.mpf(imag.numerator) / imag.denominator)


def flint_rational(value):
    """Convert a benchmark rational to an exact FLINT rational."""
    value = as_fraction(value)
    return fmpq(value.numerator, value.denominator)


def flint_number(value):
    """Convert an exact benchmark pair to an acb number."""
    real, imag = fraction_pair(value)
    return acb(flint_rational(real), flint_rational(imag))


def mp_inputs(case, points):
    """Construct mpmath inputs outside the timed region."""
    tau = [[mp_number(value) for value in row] for row in case["tau"]]
    zs = [[mp_number(value) for value in point] for point in points]
    return tau, zs


def flint_inputs(case, points):
    """Construct FLINT ball inputs outside the timed region."""
    tau = acb_mat([[flint_number(value) for value in row]
                   for row in case["tau"]])
    zs = [acb_mat([[flint_number(value)] for value in point])
          for point in points]
    return tau, zs


def characteristic_data(genus, half):
    """Return one selected mpmath/FLINT characteristic pair."""
    data = characteristic(genus) if half else None
    if data is None:
        return None, 0
    a, b = data
    mp_data = (
        [mp_number(value).real for value in a],
        [mp_number(value).real for value in b],
    )
    bits = [int(2 * fraction_pair(value)[0]) for value in a + b]
    index = 0
    for bit in bits:
        index = 2 * index + bit
    return mp_data, index


def jet_tuples(genus, order):
    """Return FLINT jet multi-indices in total-degree reverse lexicographic order."""
    result = []

    def compositions(prefix, remaining, position):
        if position == genus - 1:
            result.append(tuple(prefix + [remaining]))
            return
        for value in range(remaining, -1, -1):
            compositions(prefix + [value], remaining - value, position + 1)

    for degree in range(order + 1):
        compositions([], degree, 0)
    return result


def jet_column(genus, derivative):
    """Return a FLINT jet column and its Taylor-coefficient scale."""
    order = sum(derivative)
    column = jet_tuples(genus, order).index(tuple(derivative))
    scale = 1
    for value in derivative:
        scale *= factorial(value)
    return column, scale


def selected_point(case):
    """Choose a nonzero argument from an established benchmark case."""
    points = case_points(case)
    if case.get("points"):
        return points[0][:case["genus"]]
    point = list(points[1][:case["genus"]])
    extra = (
        number("-0.09", "0.015"),
        number("0.04", "-0.02"),
    )
    while len(point) < case["genus"]:
        point.append(extra[(len(point) - 3) % len(extra)])
    return tuple(point)


def line_points(case, count):
    """Return exact points along a line for a fixed-tau workload."""
    genus = case["genus"]
    base = selected_point(case)
    directions = (
        number("0.31", "0.07"),
        number("-0.23", "0.05"),
        number("0.17", "-0.04"),
        number("-0.11", "0.03"),
        number("0.07", "-0.02"),
    )
    result = []
    for index in range(count):
        parameter = (Fraction(0) if count == 1 else
                     Fraction(2 * index, count - 1) - 1)
        point = []
        for j in range(genus):
            base_real, base_imag = fraction_pair(base[j])
            direction_real, direction_imag = fraction_pair(
                directions[j % len(directions)])
            point.append((base_real + parameter * direction_real,
                          base_imag + parameter * direction_imag))
        result.append(tuple(point))
    return tuple(result)


def workload_data(label, genus):
    """Return whether a half-characteristic and which derivative are used."""
    zero = (0,) * genus
    if label == "value-zero":
        return False, zero
    if label == "value-half":
        return True, zero
    if label == "dz0":
        return False, (1,) + (0,) * (genus - 1)
    if label == "half-dz0":
        return True, (1,) + (0,) * (genus - 1)
    if label == "mixed2" and genus >= 2:
        return False, (1, 1) + (0,) * (genus - 2)
    if label == "third" and genus >= 2:
        return False, (2, 1) + (0,) * (genus - 2)
    raise ValueError(f"workload {label!r} is not available in genus {genus}")


def workloads(case):
    """Return bounded pointwise workloads for one case."""
    labels = ["value-zero", "value-half"]
    if case["genus"] <= 3 and (
            "coupled" in case["name"] or "unreduced" in case["name"]):
        labels.extend(("dz0", "half-dz0", "mixed2", "third"))
    return tuple(labels)


def flint_value(tau, z, characteristic_index, derivative):
    """Evaluate one selected FLINT value or derivative ball."""
    genus = tau.nrows()
    order = sum(derivative)
    if not order:
        return tau.theta(z)[0, characteristic_index]
    column, scale = jet_column(genus, derivative)
    return tau.theta_jets(z, order)[characteristic_index, column] * scale


def mpmath_value(tau, z, characteristic_value, derivative):
    """Evaluate one selected mpmath value or derivative."""
    return mp.rtheta(z, tau, characteristic_value, derivative)


def compare_value(actual, enclosure):
    """Return midpoint error, ball radius and containment information."""
    midpoint = mp.mpc(enclosure.mid())
    scale = max(mp.one, abs(midpoint))
    error = abs(actual - midpoint) / scale
    radius = mp.mpf(enclosure.rad()) / scale
    return error, radius, enclosure.contains(acb(actual))


def best_time(function, repeat, number):
    """Return the minimum warmed time per invocation."""
    function()
    timings = timeit.repeat(function, repeat=repeat, number=number)
    return min(timings) / number


def pointwise_row(case, label, dps, repeat, number):
    """Compare one pointwise workload."""
    point = selected_point(case)
    mp_tau, mp_zs = mp_inputs(case, (point,))
    flint_tau, flint_zs = flint_inputs(case, (point,))
    half, derivative = workload_data(label, case["genus"])
    mp_char, flint_char = characteristic_data(case["genus"], half)

    actual = mpmath_value(mp_tau, mp_zs[0], mp_char, derivative)
    enclosure = flint_value(
        flint_tau, flint_zs[0], flint_char, derivative)
    error, radius, contains = compare_value(actual, enclosure)
    if not contains:
        raise ArithmeticError("mpmath value is outside the FLINT enclosure")
    mp_time = best_time(
        lambda: mpmath_value(mp_tau, mp_zs[0], mp_char, derivative),
        repeat, number)
    flint_time = best_time(
        lambda: flint_value(
            flint_tau, flint_zs[0], flint_char, derivative),
        repeat, number)
    return {
        "dps": dps,
        "case": case["name"],
        "genus": case["genus"],
        "workload": label,
        "points": 1,
        "mp_time": mp_time,
        "flint_time": flint_time,
        "error": error,
        "radius": radius,
        "contains": contains,
    }


def repeated_row(case, label, dps, count, repeat):
    """Compare repeated arguments with a fixed period matrix."""
    points = line_points(case, count)
    mp_tau, mp_zs = mp_inputs(case, points)
    flint_tau, flint_zs = flint_inputs(case, points)
    half, derivative = workload_data(label, case["genus"])
    mp_char, flint_char = characteristic_data(case["genus"], half)

    actual = [mpmath_value(mp_tau, z, mp_char, derivative) for z in mp_zs]
    enclosures = [flint_value(flint_tau, z, flint_char, derivative)
                  for z in flint_zs]
    comparisons = [compare_value(value, enclosure)
                   for value, enclosure in zip(actual, enclosures)]
    if not all(item[2] for item in comparisons):
        raise ArithmeticError(
            "an mpmath value is outside its FLINT enclosure")

    def eval_mpmath():
        return [mpmath_value(mp_tau, z, mp_char, derivative) for z in mp_zs]

    def eval_flint():
        return [flint_value(flint_tau, z, flint_char, derivative)
                for z in flint_zs]

    mp_time = best_time(eval_mpmath, repeat, 1) / count
    flint_time = best_time(eval_flint, repeat, 1) / count
    return {
        "dps": dps,
        "case": case["name"],
        "genus": case["genus"],
        "workload": f"repeated-{label}",
        "points": count,
        "mp_time": mp_time,
        "flint_time": flint_time,
        "error": max(item[0] for item in comparisons),
        "radius": max(item[1] for item in comparisons),
        "contains": all(item[2] for item in comparisons),
    }


def format_number(value, digits=4):
    """Format an mpmath value compactly for Markdown."""
    return mp.nstr(value, digits, min_fixed=0, max_fixed=0)


def render_table(rows):
    """Render benchmark rows as a Markdown table."""
    lines = [
        "| dps | case | g | workload | z count | mpmath ms | FLINT ms | "
        "mpmath/FLINT | midpoint error | ball radius | enclosed |",
        "|---:|---|---:|---|---:|---:|---:|---:|---:|---:|:---:|",
    ]
    for row in rows:
        ratio = row["mp_time"] / row["flint_time"]
        lines.append(
            f"| {row['dps']} | {row['case']} | {row['genus']} | "
            f"{row['workload']} | {row['points']} | "
            f"{1000 * row['mp_time']:.4f} | "
            f"{1000 * row['flint_time']:.4f} | {ratio:.3f} | "
            f"{format_number(row['error'])} | "
            f"{format_number(row['radius'])} | "
            f"{'yes' if row['contains'] else 'no'} |"
        )
    return "\n".join(lines)


def render_report(args, pointwise, repeated, failures):
    """Return the complete Markdown report."""
    revision = "unknown"
    try:
        import subprocess
        revision = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=REPO_ROOT, text=True).strip()
        dirty = subprocess.run(
            ["git", "diff", "--quiet"], cwd=REPO_ROOT,
            check=False).returncode != 0
        if dirty:
            revision += " (dirty)"
    except (OSError, subprocess.CalledProcessError):
        pass
    lines = [
        "# Riemann theta: mpmath versus Python-FLINT",
        "",
        "The comparison uses exact rational inputs and matched decimal "
        "precision. Input and result construction are outside timed regions.",
        "",
        f"- Date: {datetime.date.today().isoformat()}",
        f"- Python: {platform.python_version()}",
        f"- mpmath: {generapy_version}",
        f"- Python-FLINT: {flint.__version__}",
        f"- FLINT: {flint.__FLINT_VERSION__}",
        f"- Revision: `{revision}`",
        f"- Command: `{' '.join(sys.argv)}`",
        f"- Repeats: {args.repeat}; pointwise calls per repeat: {args.number}",
        "- Threads: 1",
        "",
        "Python-FLINT's public `acb_mat.theta` and `theta_jets` methods "
        "return all half-characteristics. The tables extract one matching "
        "characteristic, so FLINT timings include work not requested from "
        "`rtheta`. FLINT jet Taylor coefficients are rescaled to ordinary "
        "partial derivatives before comparison.",
        "A ratio above one means that mpmath took longer than FLINT.",
        "",
        "`midpoint error` is scaled by `max(1, abs(FLINT midpoint))`. "
        "`enclosed` records whether the FLINT complex ball contains the "
        "mpmath value rounded at the same working precision.",
        "",
        "## Pointwise calls",
        "",
        render_table(pointwise) if pointwise else "No completed rows.",
        "",
        "## Fixed tau, varying z",
        "",
        "Each row evaluates a preconstructed line of arguments using the "
        "same period-matrix object. Times are reported per argument.",
        "",
        render_table(repeated) if repeated else "No completed rows.",
    ]
    if failures:
        lines.extend(("", "## Incomplete workloads", ""))
        lines.extend(f"- {failure}" for failure in failures)
    lines.append("")
    return "\n".join(lines)


def parse_args():
    """Parse command-line controls."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", nargs="+", type=int, default=[30, 50])
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--number", type=int, default=1)
    parser.add_argument("--plot-points", type=int, default=25)
    parser.add_argument("--skip-repeated", action="store_true")
    parser.add_argument(
        "--case", action="append", choices=[case["name"] for case in CASES])
    parser.add_argument(
        "--workload", action="append",
        choices=("value-zero", "value-half", "dz0", "half-dz0",
                 "mixed2", "third"))
    parser.add_argument("--output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def main():
    """Run comparisons and write the report incrementally."""
    args = parse_args()
    if args.plot_points < 1:
        raise SystemExit("--plot-points must be positive")
    selected = [case for case in CASES
                if not args.case or case["name"] in args.case]
    pointwise = []
    repeated = []
    failures = []
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        args.output.write_text(
            render_report(args, pointwise, repeated, failures),
            encoding="utf-8")

    for dps in args.dps:
        mp.dps = dps
        flint.ctx.dps = dps
        flint.ctx.threads = 1
        for case in selected:
            for label in workloads(case):
                if args.workload and label not in args.workload:
                    continue
                print(f"starting {dps} dps, {case['name']}, {label}",
                      flush=True)
                try:
                    pointwise.append(pointwise_row(
                        case, label, dps, args.repeat, args.number))
                except (ArithmeticError, ValueError) as exc:
                    failures.append(
                        f"{dps} dps, {case['name']}, {label}: {exc}")
                save()

        if not args.skip_repeated:
            for name in ("g2-coupled", "g3-coupled"):
                case = next((item for item in selected
                             if item["name"] == name), None)
                if case is None:
                    continue
                for label in ("value-zero", "dz0"):
                    if args.workload and label not in args.workload:
                        continue
                    print(f"starting repeated {dps} dps, {name}, {label}",
                          flush=True)
                    try:
                        repeated.append(repeated_row(
                            case, label, dps, args.plot_points, args.repeat))
                    except (ArithmeticError, ValueError) as exc:
                        failures.append(
                            f"repeated {dps} dps, {name}, {label}: {exc}")
                    save()

    save()
    print(f"wrote {args.output} "
          f"({len(pointwise) + len(repeated)} completed, "
          f"{len(failures)} incomplete)")
    if failures and args.strict:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
