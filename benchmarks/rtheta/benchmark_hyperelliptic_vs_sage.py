#!/usr/bin/env python3
"""Compare genera hyperelliptic periods with two Sage implementations.

Run from the repository root with the normal mpmath environment::

    .venv/bin/python \
        benchmarks/rtheta/benchmark_hyperelliptic_vs_sage.py

The comparison covers Sage's arbitrary-precision Riemann-surface machinery
(also used by the nbruin/RiemannTheta examples) and abelfunctions' independent
machine-precision period implementation. Sage startup is excluded from all
reported timings.
"""

from __future__ import annotations

import argparse
import datetime
from fractions import Fraction
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import timeit


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_REPORT = SCRIPT_DIR / "hyperelliptic_vs_sage.md"
DEFAULT_SAGE = "/usr/local/bin/sage"
sys.path.insert(0, str(REPO_ROOT/'src'))

import genera
from mpmath import mp
genera_version = genera.__version__
from genera.curves._stages import (  # noqa: E402
    _stage_hyperelliptic_periods,
)

from benchmark_hyperelliptic_vs_pari import (  # noqa: E402
    CASES as REAL_ROOT_CASES,
    geometric_mean,
    lattice_comparison,
    polynomial_coefficients,
    scientific,
)


CASES = tuple(
    (name, polynomial_coefficients(roots))
    for name, roots in REAL_ROOT_CASES
) + (
    # (x^2+1)*((x-2)^2+1)*(x-4)
    ("g2-complex", tuple(map(Fraction, (-20, 21, -28, 22, -8, 1)))),
    # (x^2+1)*((x-2)^2+1)*((x-4)^2+1)*(x-6)
    ("g3-complex", tuple(map(
        Fraction, (-510, 733, -942, 859, -450, 127, -18, 1)))),
)


SAGE_WORKER = r'''
import json
import math
import statistics
import sys
import time

from sage.all import Curve, PolynomialRing, QQ


def rational(field, pair):
    return field(pair[0]) / field(pair[1])


def serialise(matrix):
    rows = matrix.rows() if callable(getattr(matrix, "rows", None)) else matrix
    result = []
    for row in rows:
        values = []
        for value in row:
            real = value.real() if callable(value.real) else value.real
            imag = value.imag() if callable(value.imag) else value.imag
            values.append([str(real), str(imag)])
        result.append(values)
    return result


request = json.load(sys.stdin)
implementation = request["implementation"]
if implementation == "abelfunctions":
    from abelfunctions import RiemannSurface as AbelRiemannSurface

ring = PolynomialRing(QQ, names=("x", "y"))
x, y = ring.gens()

# Import and one-time Sage initialization are intentionally outside timings.
warmup_polynomial = x**3 - x
if implementation == "sage":
    Curve(y**2 - warmup_polynomial).riemann_surface(
        prec=53).period_matrix()
else:
    AbelRiemannSurface(y**2 - warmup_polynomial).period_matrix()

rows = []
for case in request["cases"]:
    dps = case["dps"]
    bits = max(53, math.ceil(dps * math.log2(10)))
    coefficients = [rational(QQ, value) for value in case["coefficients"]]
    polynomial = sum(value * x**power
                     for power, value in enumerate(coefficients))

    def evaluate():
        if implementation == "sage":
            surface = Curve(y**2 - polynomial).riemann_surface(
                prec=bits, integration_method="rigorous")
        else:
            surface = AbelRiemannSurface(y**2 - polynomial)
        return surface.period_matrix()

    samples = []
    matrix = None
    for unused in range(request["trials"]):
        start = time.perf_counter()
        matrix = evaluate()
        samples.append(time.perf_counter() - start)
    rows.append({
        "implementation": implementation,
        "case": case["name"],
        "genus": (len(coefficients) - 2) // 2,
        "dps": dps,
        "seconds": statistics.median(samples),
        "matrix": serialise(matrix),
    })

print("HYPERELL_SAGE_JSON=" + json.dumps(rows))
'''


def fraction_pair(value):
    """Encode a Fraction for lossless transfer to Sage."""
    return value.numerator, value.denominator


def sage_request(implementation, digits, trials, selected_cases):
    """Build bounded exact input for one Sage-side implementation."""
    selected_digits = (15,) if implementation == "abelfunctions" else digits
    cases = []
    for dps in selected_digits:
        for name, coefficients in CASES:
            if name not in selected_cases:
                continue
            cases.append({
                "name": name,
                "dps": dps,
                "coefficients": [fraction_pair(value)
                                 for value in coefficients],
            })
    return {"implementation": implementation, "trials": trials, "cases": cases}


def run_sage(sage, implementation, digits, trials, selected_cases, timeout):
    """Run one Sage worker and parse its machine-readable output."""
    request = sage_request(
        implementation, digits, trials, selected_cases)
    environment = os.environ.copy()
    dot_sage = SCRIPT_DIR / ".sage-period-benchmark"
    dot_sage.mkdir(exist_ok=True)
    environment["DOT_SAGE"] = str(dot_sage)
    result = subprocess.run(
        [sage, "-python", "-c", SAGE_WORKER],
        input=json.dumps(request), text=True, capture_output=True,
        cwd="/tmp", env=environment, timeout=timeout)
    if result.returncode:
        raise RuntimeError(
            f"{implementation} worker failed:\n" + result.stdout + result.stderr)
    marker = "HYPERELL_SAGE_JSON="
    for line in reversed(result.stdout.splitlines()):
        if line.startswith(marker):
            return json.loads(line[len(marker):])
    raise RuntimeError(
        f"could not parse {implementation} output:\n"
        + result.stdout + result.stderr)


def decoded_matrix(data):
    """Convert a serialised Sage matrix to an mpmath matrix."""
    return mp.matrix([
        [mp.mpc(real, imag) for real, imag in row]
        for row in data
    ])


def mpmath_periods(coefficients):
    """Return the genera first-kind half-period matrix."""
    from genera import algebraic_curve
    _stage_hyperelliptic_periods.cache_clear()
    data = algebraic_curve(coefficients).periods_kind_1()
    omega, omega_prime = data.omega, data.omega_prime
    periods = mp.matrix(omega.rows, 2 * omega.cols)
    periods[:, :omega.cols] = omega
    periods[:, omega.cols:] = omega_prime
    return periods


def benchmark(args):
    """Run selected Sage implementations and compare with local mpmath."""
    sage = str(Path(args.sage).resolve())
    if not Path(sage).is_file():
        raise SystemExit(f"Sage executable not found: {args.sage}")
    sage_rows = []
    for implementation in args.implementations:
        sage_rows.extend(run_sage(
            sage, implementation, tuple(args.digits), args.sage_trials,
            set(args.cases), args.sage_timeout))

    rows = []
    for sage_row in sage_rows:
        name = sage_row["case"]
        exact = next(coefficients for case_name, coefficients in CASES
                     if case_name == name)
        dps = sage_row["dps"]
        with mp.workdps(dps):
            coefficients = [mp.mpf(value.numerator) / value.denominator
                            for value in exact]
            mpmath_periods(coefficients)
            samples = timeit.repeat(
                lambda: mpmath_periods(coefficients),
                repeat=args.mpmath_trials, number=1)
            ours = mpmath_periods(coefficients)
            # Sage and abelfunctions use x^k dx/(df/dy) = x^k dx/(2y),
            # whose complete cycle integrals equal our half-periods.
            theirs = decoded_matrix(sage_row["matrix"])
            residual, integer_error, orientation, form_error = (
                lattice_comparison(ours, theirs))
        rows.append({
            **sage_row,
            "mpmath_seconds": statistics.median(samples),
            "ratio": statistics.median(samples) / sage_row["seconds"],
            "residual": residual,
            "integer_error": integer_error,
            "orientation": orientation,
            "form_error": form_error,
        })
    return rows, sage


def markdown(rows, args, sage):
    """Render the Sage comparisons as a Markdown report."""
    command = (
        ".venv/bin/python "
        "benchmarks/rtheta/benchmark_hyperelliptic_vs_sage.py "
        "--cases " + " ".join(args.cases))
    lines = [
        "# Hyperelliptic period comparison with Sage",
        "",
        "This compares mpmath with Sage's arbitrary-precision Riemann-surface "
        "period calculation and abelfunctions' machine-precision calculation. "
        "The former is the period engine used by the nbruin/RiemannTheta "
        "README; RiemannTheta itself begins with an already computed period "
        "matrix.",
        "",
        "Sage and abelfunctions represent the standard algebraic differential "
        "as `x^k dx/(df/dy) = x^k dx/(2y)`. Their matrices are therefore "
        "multiplied by the known factor two before recovering the nearest "
        "integral change of cycle basis. No numerical scale is fitted.",
        "",
        "## Reproduction",
        "",
        "```console",
        command,
        "```",
        "",
        f"- Date: {datetime.datetime.now(datetime.timezone.utc).isoformat()}",
        f"- Platform: {platform.platform()}",
        f"- Python: {platform.python_version()}",
        f"- mpmath: {genera_version}",
        f"- Sage: `{sage}`",
        f"- mpmath timing: median of {args.mpmath_trials} trials",
        f"- Sage timing: median of {args.sage_trials} trials inside one Sage "
        "process; interpreter startup and the initial warmup are excluded",
        "",
        "## Accuracy",
        "",
        "| Implementation | Case | Genus | dps | Relative lattice residual | Integer-map error | Orientation | Form error |",
        "|---|---|---:|---:|---:|---:|---|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['implementation']} | {row['case']} | {row['genus']} | "
            f"{row['dps']} | {scientific(row['residual'])} | "
            f"{scientific(row['integer_error'])} | {row['orientation']} | "
            f"{scientific(row['form_error'])} |")
    lines.extend([
        "",
        "abelfunctions stores periods in NumPy complex128 arrays, so its dps "
        "column denotes the approximately 15-digit comparison target rather "
        "than user-selectable arbitrary precision.",
        "",
        "## Performance",
        "",
        "The timed workload is the complete curve-to-period calculation, "
        "including topology/root work and integration but excluding imports "
        "and construction of the polynomial input.",
        "",
        "| Implementation | Case | Genus | dps | mpmath (ms) | Other (ms) | mpmath / other |",
        "|---|---|---:|---:|---:|---:|---:|",
    ])
    for row in rows:
        lines.append(
            f"| {row['implementation']} | {row['case']} | {row['genus']} | "
            f"{row['dps']} | {1000 * row['mpmath_seconds']:.3f} | "
            f"{1000 * row['seconds']:.3f} | {row['ratio']:.3g}x |")
    lines.extend([
        "",
        "### Geometric means",
        "",
        "| Implementation | dps | mpmath (ms) | Other (ms) | mpmath / other |",
        "|---|---:|---:|---:|---:|",
    ])
    groups = sorted({(row["implementation"], row["dps"]) for row in rows})
    for implementation, dps in groups:
        selected = [row for row in rows
                    if row["implementation"] == implementation
                    and row["dps"] == dps]
        lines.append(
            f"| {implementation} | {dps} | "
            f"{1000 * geometric_mean([row['mpmath_seconds'] for row in selected]):.3f} | "
            f"{1000 * geometric_mean([row['seconds'] for row in selected]):.3f} | "
            f"{geometric_mean([row['ratio'] for row in selected]):.3g}x |")
    lines.append("")
    return "\n".join(lines)


def parse_args():
    """Parse deliberately bounded benchmark options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sage", default=DEFAULT_SAGE)
    parser.add_argument("--implementations", nargs="+",
                        choices=("sage", "abelfunctions"),
                        default=("sage", "abelfunctions"))
    parser.add_argument("--digits", nargs="+", type=int,
                        default=(30, 50, 100))
    parser.add_argument(
        "--cases", nargs="+", choices=tuple(name for name, unused in CASES),
        default=("g1-symmetric", "g2-irregular"))
    parser.add_argument("--sage-trials", type=int, default=1)
    parser.add_argument("--sage-timeout", type=float, default=90)
    parser.add_argument("--mpmath-trials", type=int, default=3)
    parser.add_argument("--output", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def main():
    """Run the comparison and write its Markdown report."""
    args = parse_args()
    rows, sage = benchmark(args)
    report = markdown(rows, args, sage)
    args.output.write_text(report, encoding="utf-8")
    print(report)
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
