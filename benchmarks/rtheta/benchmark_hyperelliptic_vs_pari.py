#!/usr/bin/env python3
"""Compare genera hyperelliptic periods with PARI/GP.

Run from the repository root with the normal mpmath development environment::

    .venv/bin/python \
        benchmarks/rtheta/benchmark_hyperelliptic_vs_pari.py

A PARI build exposing ``hyperellperiods`` is required (tested with the
2.18.1 alpha). Set ``PARI_GP`` when ``gp`` is not on PATH. Interpreter startup
and input construction are excluded from both timings.
"""

from __future__ import annotations

import argparse
import datetime
from fractions import Fraction
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import sys
import timeit


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_REPORT = SCRIPT_DIR / "hyperelliptic_vs_pari.md"
DEFAULT_PARI = Path.home() / ".local/pari-2.18.1-alpha/bin/gp"
sys.path.insert(0, str(REPO_ROOT/'src'))

import genera
from mpmath import mp
genera_version = genera.__version__
from genera.curves._stages import (  # noqa: E402
    _stage_hyperelliptic_periods,
)


CASES = (
    ("g1-symmetric", (Fraction(-1), Fraction(0), Fraction(1))),
    ("g1-even", tuple(map(Fraction, (-3, -1, 2, 4)))),
    ("g2-symmetric", tuple(map(Fraction, (-2, -1, 0, 1, 2)))),
    ("g2-irregular", (
        Fraction(-3), Fraction(-1), Fraction(1, 2), Fraction(2), Fraction(5))),
    ("g2-even", tuple(map(Fraction, (-3, -2, -1, 1, 2, 4)))),
    ("g3-symmetric", tuple(map(Fraction, (-3, -2, -1, 0, 1, 2, 3)))),
    ("g3-irregular", (
        Fraction(-4), Fraction(-5, 2), Fraction(-7, 10), Fraction(1, 5),
        Fraction(11, 10), Fraction(3), Fraction(6))),
)


def polynomial_coefficients(roots):
    """Return exact ascending coefficients of a monic polynomial."""
    coefficients = [Fraction(1)]
    for root in roots:
        result = [Fraction(0)] * (len(coefficients) + 1)
        for index, value in enumerate(coefficients):
            result[index] -= root * value
            result[index + 1] += value
        coefficients = result
    return tuple(coefficients)


def gp_rational(value):
    """Format a Fraction as exact GP input."""
    if value.denominator == 1:
        return str(value.numerator)
    return f"({value.numerator}/{value.denominator})"


def gp_polynomial(coefficients):
    """Format ascending coefficients as a GP polynomial."""
    terms = []
    for power, value in enumerate(coefficients):
        if not value:
            continue
        coefficient = gp_rational(value)
        if power == 0:
            terms.append(coefficient)
        elif power == 1:
            terms.append(f"({coefficient})*x")
        else:
            terms.append(f"({coefficient})*x^{power}")
    return "+".join(terms)


def find_gp(argument):
    """Resolve the requested or locally installed PARI/GP executable."""
    candidates = (argument, os.environ.get("PARI_GP"),
                  str(DEFAULT_PARI) if DEFAULT_PARI.exists() else None,
                  shutil.which("gp"))
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return str(Path(candidate).resolve())
    raise SystemExit(
        "PARI/GP 2.18.1 or newer was not found; pass --pari-gp or set PARI_GP")


def run_pari(gp, coefficients, digits, repetitions):
    """Compute and internally time PARI hyperellperiods in one GP process."""
    polynomial = gp_polynomial(coefficients)
    script = f"""\\p {digits}
P={polynomial};
M=hyperellperiods(P,1);
t=getwalltime();
for(k=1,{repetitions},M=hyperellperiods(P,1));
elapsed=getwalltime()-t;
print(\"VERSION|\",version());
print(\"TIME|\",elapsed,\"|{repetitions}\");
for(i=1,matsize(M)[1],for(j=1,matsize(M)[2],print(\"ENTRY|\",i,\"|\",j,\"|\",real(M[i,j]),\"|\",imag(M[i,j]))));
quit
"""
    result = subprocess.run(
        [gp, "-q", "-f"], input=script, text=True,
        capture_output=True, check=True)
    version = None
    elapsed = None
    actual_repetitions = repetitions
    entries = {}
    for line in result.stdout.splitlines():
        if line.startswith("VERSION|"):
            version = line.split("|", 1)[1]
        elif line.startswith("TIME|"):
            _, milliseconds, count = line.split("|")
            elapsed = float(milliseconds) / int(count)
            actual_repetitions = int(count)
        elif line.startswith("ENTRY|"):
            _, row, column, real, imag = line.split("|")
            entries[int(row) - 1, int(column) - 1] = mp.mpc(real, imag)
    if version is None or elapsed is None or not entries:
        raise RuntimeError(
            "could not parse PARI output:\n" + result.stdout + result.stderr)
    genus = (len(coefficients) - 2) // 2
    matrix = mp.matrix(genus, 2 * genus)
    # PARI flag 1 interleaves a- and b-period columns.
    for row in range(genus):
        for column in range(genus):
            matrix[row, column] = entries[row, 2 * column]
            matrix[row, genus + column] = entries[row, 2 * column + 1]
    return matrix, elapsed, actual_repetitions, version


def realification(matrix):
    """Stack real and imaginary parts of a complex period matrix."""
    rows, columns = matrix.rows, matrix.cols
    result = mp.matrix(2 * rows, columns)
    for row in range(rows):
        for column in range(columns):
            result[row, column] = mp.re(matrix[row, column])
            result[rows + row, column] = mp.im(matrix[row, column])
    return result


def lattice_comparison(mpmath_periods, pari_periods):
    """Recover and assess the integral change of homology basis."""
    source = realification(mpmath_periods)
    target = realification(pari_periods)
    approximate = mp.inverse(source) * target
    transform = approximate.apply(mp.nint)
    residual = mp.norm(mpmath_periods * transform - pari_periods) / max(
        mp.one, mp.norm(pari_periods))
    integer_error = mp.norm(approximate - transform)
    genus = mpmath_periods.rows
    symplectic_form = mp.zeros(2 * genus)
    for index in range(genus):
        symplectic_form[index, genus + index] = 1
        symplectic_form[genus + index, index] = -1
    preserved = transform.T * symplectic_form * transform
    symplectic_error = mp.norm(preserved - symplectic_form)
    antisymplectic_error = mp.norm(preserved + symplectic_form)
    if antisymplectic_error < symplectic_error:
        orientation = "opposite"
        form_error = antisymplectic_error
    else:
        orientation = "same"
        form_error = symplectic_error
    return residual, integer_error, orientation, form_error


def geometric_mean(values):
    """Return the geometric mean of positive finite values."""
    values = [value for value in values if value > 0 and math.isfinite(value)]
    return math.exp(statistics.fmean(math.log(value) for value in values))


def scientific(value, digits=3):
    """Format an mpmath or float value compactly."""
    return mp.nstr(value, digits + 1, min_fixed=0, max_fixed=0)


def benchmark(args):
    """Run all selected cases and return result rows plus metadata."""
    gp = find_gp(args.pari_gp)
    rows = []
    pari_version = None
    for digits in args.digits:
        with mp.workdps(digits):
            for name, roots in CASES:
                coefficients_exact = polynomial_coefficients(roots)
                coefficients = [
                    mp.mpf(value.numerator) / value.denominator
                    for value in coefficients_exact]
                def mpmath_period_data():
                    from genera import Curve
                    _stage_hyperelliptic_periods.cache_clear()
                    return Curve(
                        {(0, 2): 1, **{(i, 0): -c for i, c in enumerate(coefficients)}}).periods_kind_1()

                # Warm up root finding and quadrature allocation paths.
                mpmath_period_data()
                samples = timeit.repeat(
                    mpmath_period_data,
                    repeat=args.trials, number=args.repetitions)
                mpmath_ms = 1000 * statistics.median(samples) / args.repetitions
                pari_periods, pari_ms, pari_repetitions, pari_version = run_pari(
                    gp, coefficients_exact, digits, args.pari_repetitions)
                data = mpmath_period_data()
                omega, omega_prime = data.omega, data.omega_prime
                periods = mp.matrix(omega.rows, 2 * omega.cols)
                for row in range(omega.rows):
                    for column in range(omega.cols):
                        periods[row, column] = omega[row, column]
                        periods[row, omega.cols + column] = omega_prime[row, column]
                residual, integer_error, orientation, form_error = (
                    lattice_comparison(periods, pari_periods / 2))
                rows.append({
                    "name": name,
                    "genus": omega.rows,
                    "digits": digits,
                    "mpmath_ms": mpmath_ms,
                    "pari_ms": pari_ms,
                    "ratio": mpmath_ms / pari_ms if pari_ms else math.inf,
                    "residual": residual,
                    "integer_error": integer_error,
                    "orientation": orientation,
                    "form_error": form_error,
                    "pari_repetitions": pari_repetitions,
                })
    return rows, gp, pari_version


def markdown(rows, args, gp, pari_version):
    """Render a reproducible Markdown benchmark report."""
    command = (
        ".venv/bin/python "
        "benchmarks/rtheta/benchmark_hyperelliptic_vs_pari.py")
    lines = [
        "# Hyperelliptic period comparison with PARI/GP",
        "",
        "This compares mpmath's pure-Python period calculation with PARI's "
        "`hyperellperiods` at matched decimal precision. PARI chooses a "
        "different homology basis, so accuracy is measured after recovering "
        "the nearest integral change of basis for the complete period lattice.",
        "",
        "## Reproduction",
        "",
        f"```console\n{command}\n```",
        "",
        f"- Date: {datetime.datetime.now(datetime.timezone.utc).isoformat()}",
        f"- Platform: {platform.platform()}",
        f"- Python: {platform.python_version()}",
        f"- mpmath: {genera_version}",
        f"- PARI/GP: {pari_version} (`{gp}`)",
        f"- mpmath timing: median of {args.trials} trials, "
        f"{args.repetitions} call(s) per trial",
        f"- PARI timing: internal `getwalltime`, "
        f"{args.pari_repetitions} call(s) per case (interpreter startup excluded)",
        "",
        "## Accuracy",
        "",
        "| Case | Genus | dps | Relative lattice residual | Integer-map error | Orientation | Form error |",
        "|---|---:|---:|---:|---:|---|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['name']} | {row['genus']} | {row['digits']} | "
            f"{scientific(row['residual'])} | "
            f"{scientific(row['integer_error'])} | {row['orientation']} | "
            f"{scientific(row['form_error'])} |")
    lines.extend([
        "",
        "`opposite` means that the integral basis change reverses the sign of "
        "the standard symplectic form. This reflects the two implementations' "
        "opposite cycle-orientation conventions, not an accuracy discrepancy.",
        "",
        "## Performance",
        "",
        "| Case | Genus | dps | mpmath (ms) | PARI (ms) | mpmath / PARI |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    for row in rows:
        lines.append(
            f"| {row['name']} | {row['genus']} | {row['digits']} | "
            f"{row['mpmath_ms']:.3f} | {row['pari_ms']:.3f} | "
            f"{row['ratio']:.2f}x |")
    lines.extend([
        "",
        "### Geometric means",
        "",
        "| dps | mpmath (ms) | PARI (ms) | mpmath / PARI |",
        "|---:|---:|---:|---:|",
    ])
    for digits in args.digits:
        selected = [row for row in rows if row["digits"] == digits]
        lines.append(
            f"| {digits} | "
            f"{geometric_mean([row['mpmath_ms'] for row in selected]):.3f} | "
            f"{geometric_mean([row['pari_ms'] for row in selected]):.3f} | "
            f"{geometric_mean([row['ratio'] for row in selected]):.2f}x |")
    lines.extend([
        "",
        "PARI/GP is GPL-2.0-or-later. This script invokes an independently "
        "installed executable and records numerical output; it does not copy "
        "or distribute PARI source code. For provenance and reproducibility, "
        "reports using these results should cite the tested PARI/GP version.",
        "",
    ])
    return "\n".join(lines)


def parse_args():
    """Parse bounded benchmark options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--digits", nargs="+", type=int, default=(30, 50, 100))
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--repetitions", type=int, default=1)
    parser.add_argument("--pari-repetitions", type=int, default=20)
    parser.add_argument("--pari-gp")
    parser.add_argument("--output", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def main():
    """Run the benchmark and write its Markdown report."""
    args = parse_args()
    rows, gp, pari_version = benchmark(args)
    report = markdown(rows, args, gp, pari_version)
    args.output.write_text(report, encoding="utf-8")
    print(report)
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
