#!/usr/bin/env python3
"""Compare genera rtheta with Wolfram Engine's SiegelTheta.

The script evaluates genus-one, genus-two and genus-three cases, checks the
returned values, times warmed fixed-tau calls, and writes a Markdown report
next to this file.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import io
import json
import math
import os
import platform
import re
import statistics
import subprocess
import sys
import timeit
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_REPORT = SCRIPT_DIR / "rtheta_vs_wolfram.md"
WOLFRAMSCRIPT = "/usr/local/bin/wolframscript"
WOLFRAM_KERNEL = (
    "/Applications/Wolfram Engine.app/Contents/Resources/"
    "Wolfram Player.app/Contents/MacOS/WolframKernel"
)
sys.path.insert(0, str(REPO_ROOT/'src'))

import genera
from mpmath import mp  # noqa: E402


def number(real, imag="0"):
    """Describe a complex number using exact decimal strings."""
    return real, imag


CASES = (
    {
        "name": "g1",
        "genus": 1,
        "kind": "genus-one control",
        "tau": ((number("0.20", "1.20"),),),
    },
    {
        "name": "g2-diagonal",
        "genus": 2,
        "kind": "factorisation control",
        "tau": (
            (number("0.20", "1.00"), number("0")),
            (number("0"), number("-0.15", "1.20")),
        ),
    },
    {
        "name": "g2-coupled",
        "genus": 2,
        "kind": "coupled matrix",
        "tau": (
            (number("0.10", "0.90"), number("-0.20", "0.12")),
            (number("-0.20", "0.12"), number("0.05", "1.10")),
        ),
        "workloads": ("zero", "half", "dz0", "half-dz0"),
    },
    {
        "name": "g3-diagonal",
        "genus": 3,
        "kind": "factorisation control",
        "tau": (
            (number("0.10", "0.90"), number("0"), number("0")),
            (number("0"), number("-0.20", "1.10"), number("0")),
            (number("0"), number("0"), number("0.05", "1.25")),
        ),
    },
    {
        "name": "g3-coupled",
        "genus": 3,
        "kind": "coupled matrix",
        "tau": (
            (number("0.10", "1.00"), number("-0.12", "0.08"),
             number("0.05", "0.04")),
            (number("-0.12", "0.08"), number("0.20", "1.15"),
             number("0.08", "0.06")),
            (number("0.05", "0.04"), number("0.08", "0.06"),
             number("-0.15", "0.90")),
        ),
        "workloads": ("zero", "half", "dz0", "half-dz0"),
        "derivative_max_dps": 50,
    },
    {
        "name": "g2-near-boundary",
        "genus": 2,
        "kind": "near-boundary direct-path control",
        "tau": (
            (number("0.20", "0.92"), number("-0.12", "0.06")),
            (number("-0.12", "0.06"), number("-0.18", "1.08")),
        ),
    },
    {
        "name": "g2-unreduced",
        "genus": 2,
        "kind": "strongly unreduced matrix",
        "tau": (
            (number("0.20", "0.15"), number("0.12", "0.03")),
            (number("0.12", "0.03"), number("-0.10", "0.20")),
        ),
        "workloads": ("zero", "half", "real"),
        "points": (
            (number("0.13", "0.07"), number("-0.21", "0.04")),
            (number("2.31", "-1.14"), number("-1.72", "0.83")),
        ),
    },
    {
        "name": "g3-unreduced",
        "genus": 3,
        "kind": "strongly unreduced matrix",
        "tau": (
            (number("0.30", "0.12"), number("0.17", "0.03"),
             number("0", "0.02")),
            (number("0.17", "0.03"), number("-0.20", "0.16"),
             number("0", "0.02")),
            (number("0", "0.02"), number("0", "0.02"),
             number("0.10", "0.20")),
        ),
        "workloads": ("zero", "half", "real"),
    },
    {
        "name": "g3-block",
        "genus": 3,
        "kind": "genus-one plus genus-two block control",
        "tau": (
            (number("0.15", "0.85"), number("0"), number("0")),
            (number("0"), number("0.10", "0.95"),
             number("-0.07", "0.04")),
            (number("0"), number("-0.07", "0.04"),
             number("-0.12", "1.15")),
        ),
    },
)


BASE_POINTS = (
    (number("0.00"), number("0.00"), number("0.00")),
    (number("0.11", "0.03"), number("-0.07", "0.02"),
     number("0.05", "-0.01")),
    (number("-0.18", "0.06"), number("0.13", "-0.04"),
     number("-0.09", "0.02")),
    (number("0.27", "-0.08"), number("-0.16", "0.05"),
     number("0.12", "0.04")),
)


def characteristic(genus):
    """Return a deterministic nonzero half-characteristic."""
    a = [number("0.5") if index % 2 == 0 else number("0")
         for index in range(genus)]
    b = [number("0.5") if index == genus - 1 else number("0")
         for index in range(genus)]
    return a, b


def real_characteristic(genus):
    """Return a deterministic non-half-integer characteristic."""
    a_values = ("0.3", "-0.2", "0.1")
    b_values = ("0.1", "0.4", "-0.15")
    return ([number(a_values[index]) for index in range(genus)],
            [number(b_values[index]) for index in range(genus)])


def case_points(case):
    """Return the varying-z workload for one case."""
    return case.get("points", BASE_POINTS)


def workloads(case):
    """Return the requested value and derivative workloads for a case."""
    return case.get("workloads", ("zero", "half"))


def workload_points(case, label):
    """Return points for a workload, omitting trivial derivative zeros."""
    points = case_points(case)
    return points[1:] if "dz" in label else points


def workload_data(label, genus):
    """Return characteristic and derivative data for a workload label."""
    derivative = (0,) * genus
    if label == "zero":
        return None, derivative
    if label == "half":
        return characteristic(genus), derivative
    if label == "real":
        return real_characteristic(genus), derivative
    if label == "dz0":
        return None, (1,) + (0,) * (genus - 1)
    if label == "half-dz0":
        return characteristic(genus), (1,) + (0,) * (genus - 1)
    raise ValueError("unknown workload %s" % label)


def git_revision():
    """Return the checked-out revision and dirty-tree state."""
    try:
        revision = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=REPO_ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        dirty = subprocess.run(
            ["git", "diff", "--quiet"], cwd=REPO_ROOT, check=False
        ).returncode != 0
        untracked = subprocess.check_output(
            ["git", "ls-files", "--others", "--exclude-standard"],
            cwd=REPO_ROOT,
            text=True,
        ).strip()
        return revision + (" (dirty)" if dirty or untracked else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def mp_number(value):
    """Convert a decimal pair to an mpmath scalar."""
    real, imag = value
    return mp.mpc(mp.mpf(real), mp.mpf(imag))


def mp_case(case):
    """Construct one case at the active mpmath precision."""
    genus = case["genus"]
    tau = [[mp_number(value) for value in row] for row in case["tau"]]
    points = [[mp_number(value) for value in point[:genus]]
              for point in case_points(case)]
    return tau, points


def mp_characteristic(data):
    """Convert characteristic decimal pairs to mpmath scalars."""
    if data is None:
        return None
    a_data, b_data = data
    return ([mp_number(value).real for value in a_data],
            [mp_number(value).real for value in b_data])


def wl_real(text, digits):
    """Represent a decimal as an exact rational, then set precision later."""
    # The precision mark prevents the decimal from first becoming a binary
    # machine number. Rationalize then preserves the intended decimal value.
    return "Rationalize[%s`%d, 0]" % (text, digits)


def wl_number(value, digits):
    """Represent a complex decimal pair in Wolfram Language syntax."""
    real, imag = value
    if mp.mpf(imag):
        return "(%s + I (%s))" % (wl_real(real, digits),
                                   wl_real(imag, digits))
    return wl_real(real, digits)


def wl_vector(values, digits):
    """Represent a vector in Wolfram Language syntax."""
    return "{%s}" % ", ".join(wl_number(value, digits) for value in values)


def wl_matrix(rows, digits):
    """Represent a matrix in Wolfram Language syntax."""
    return "{%s}" % ", ".join(wl_vector(row, digits) for row in rows)


def parse_wl_real(text):
    """Parse a real InputForm number emitted by Wolfram Language."""
    text = text.strip().replace("*^", "e")
    text = re.sub(r"`[0-9.]*", "", text)
    return mp.mpf(text)


def wl_function(characteristic_data, derivative, genus, digits):
    """Return Wolfram code for one value or partial-derivative function."""
    if characteristic_data is None:
        theta = "SiegelTheta[tau, %s]"
    else:
        a, b = characteristic_data
        char = "{%s, %s}" % (wl_vector(a, digits), wl_vector(b, digits))
        theta = "SiegelTheta[%s, tau, %%s]" % char
    if not any(derivative):
        return "Function[z, %s]" % (theta % "z")
    if sum(derivative) != 1:
        raise ValueError("Wolfram oracle currently supports first derivatives")
    index = derivative.index(1)
    variable = "x%d" % (index + 1)
    shifted = "ReplacePart[z, %d -> %s]" % (index + 1, variable)
    return (
        "Function[z, Quiet[NumericalCalculus`ND[%s, %s, z[[%d]], "
        "Method -> NIntegrate, WorkingPrecision -> %d]]]"
        % (theta % shifted, variable, index + 1, digits + 30)
    )


def wolfram_row(case, label, digits, batches, target_seconds, timeout,
                 correctness_only=False):
    """Evaluate one workload in an isolated, time-limited Wolfram process."""
    genus = case["genus"]
    tau = wl_matrix(case["tau"], digits)
    points = "{%s}" % ", ".join(
        wl_vector(point[:genus], digits)
        for point in workload_points(case, label)
    )
    characteristic_data, derivative = workload_data(label, genus)
    function = wl_function(characteristic_data, derivative, genus, digits)
    working_digits = digits + 30 if any(derivative) else digits
    timing = "False" if correctness_only else "True"
    code = r"""
result = TimeConstrained[
  Block[{$MaxExtraPrecision = 10000},
    Module[{tau, zs, f, values, elapsed},
      Needs["NumericalCalculus`"];
      tau = N[%(tau)s, %(working_digits)d];
      zs = N[%(points)s, %(working_digits)d];
      f = %(function)s;
      values = N[f /@ zs, %(digits)d];
      f /@ zs;
      elapsed = If[%(timing)s,
        First@RepeatedTiming[
          Do[Scan[f, zs], {%(batches)d}], %(target)s
        ]/(%(batches)d Length[zs]),
        Missing["NotTimed"]
      ];
      Join[
        {%(name)s, %(label)s,
         If[MissingQ[elapsed], "NA",
            ToString[N[elapsed, 17], InputForm]]},
        (ToString[#, InputForm] & /@
          Flatten[({Re[#], Im[#]} & /@ values)])
      ]
    ]
  ], %(timeout)s, $Failed];
If[result === $Failed, "RT_BENCHMARK_TIMEOUT", ExportString[{result}, "CSV"]]
""" % {
        "tau": tau,
        "digits": digits,
        "working_digits": working_digits,
        "points": points,
        "function": function,
        "timing": timing,
        "batches": batches,
        "target": target_seconds,
        "name": json.dumps(case["name"]),
        "label": json.dumps(label),
        "timeout": timeout,
    }
    environment = os.environ.copy()
    environment["WolframKernel"] = WOLFRAM_KERNEL
    try:
        process = subprocess.run(
            [WOLFRAMSCRIPT, "-code", code],
            cwd=REPO_ROOT,
            env=environment,
            text=True,
            capture_output=True,
            check=True,
            timeout=timeout + 15,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError("Wolfram process timed out") from exc
    if "RT_BENCHMARK_TIMEOUT" in process.stdout:
        raise RuntimeError("Wolfram evaluation timed out")
    rows = [row for row in csv.reader(io.StringIO(process.stdout)) if row]
    if len(rows) != 1 or rows[0][:2] != [case["name"], label]:
        raise RuntimeError(
            "unexpected Wolfram output:\nstdout=%s\nstderr=%s"
            % (process.stdout, process.stderr)
        )
    row = rows[0]
    elapsed = None if row[2] == "NA" else float(row[2])
    values = [
        mp.mpc(parse_wl_real(row[index]), parse_wl_real(row[index + 1]))
        for index in range(3, len(row), 2)
    ]
    return elapsed, values


def median_time(function, number, repeat, calls_per_batch):
    """Return median wall time per mathematical function call."""
    samples = timeit.repeat(function, number=number, repeat=repeat)
    return statistics.median(samples) / (number * calls_per_batch)


def scaled_error(actual, expected):
    """Return absolute error scaled by max(1, abs(expected))."""
    return abs(actual - expected) / max(mp.one, abs(expected))


def diagonal_product(tau, z, char, derivative):
    """Evaluate the product of genus-one factors for diagonal tau."""
    factors = []
    for index in range(len(z)):
        component_char = ([char[0][index]], [char[1][index]])
        factors.append(mp.rtheta(
            [z[index]], [[tau[index][index]]], component_char,
            derivative[index]
        ))
    return mp.fprod(factors)


def block_product(tau, z, char, derivative):
    """Evaluate the genus-one times genus-two block factorisation."""
    first_char = ([char[0][0]], [char[1][0]])
    second_char = (char[0][1:], char[1][1:])
    return mp.rtheta(
        z[:1], [[tau[0][0]]], first_char, derivative[0]
    ) * mp.rtheta(
        z[1:], [row[1:] for row in tau[1:]], second_char, derivative[1:]
    )


def mpmath_row(case, label, digits, number, repeat,
                correctness_only=False):
    """Evaluate and optionally time one mpmath workload."""
    with mp.workdps(digits):
        tau, unused_points = mp_case(case)
        points = [[mp_number(value) for value in point[:case["genus"]]]
                  for point in workload_points(case, label)]
        characteristic_data, derivative = workload_data(label, case["genus"])
        selected_char = mp_characteristic(characteristic_data)
        values = [mp.rtheta(z, tau, selected_char, derivative)
                  for z in points]

        def evaluate_batch():
            for z in points:
                mp.rtheta(z, tau, selected_char, derivative)

        elapsed = None
        if not correctness_only:
            evaluate_batch()
            elapsed = median_time(
                evaluate_batch, number, repeat, len(points)
            )
        factor_error = None
        actual_char = selected_char or (
            [mp.zero] * case["genus"], [mp.zero] * case["genus"]
        )
        if case["kind"] == "factorisation control":
            factor_error = max(
                scaled_error(value,
                             diagonal_product(tau, z, actual_char,
                                              derivative))
                for z, value in zip(points, values)
            )
        elif case["name"] == "g3-block":
            factor_error = max(
                scaled_error(value,
                             block_product(tau, z, actual_char, derivative))
                for z, value in zip(points, values)
            )
    return elapsed, values, factor_error


def short_complex(value, digits=10):
    """Format a compact complex value for Markdown."""
    sign = "+" if value.imag >= 0 else "-"
    return "%s %s %sj" % (mp.nstr(value.real, digits), sign,
                           mp.nstr(abs(value.imag), digits))


def geometric_mean(values):
    """Return the geometric mean of positive values."""
    return math.exp(statistics.mean(math.log(value) for value in values))


def render_report(args, results, failures):
    """Render all benchmark results as Markdown."""
    timed = [row for row in results
             if row["mp_time"] is not None and row["wl_time"] is not None]
    all_ratios = [row["mp_time"] / row["wl_time"] for row in timed]
    all_errors = [row["error"] for row in results]
    lines = [
        "# Riemann theta baseline: mpmath `rtheta` vs Wolfram `SiegelTheta`",
        "",
        "This compares the current generic mpmath direct sum with Wolfram "
        "Engine's multidimensional `SiegelTheta`. It includes genus 1, 2 "
        "and 3, zero and nonzero half-characteristics, diagonal "
        "factorisation controls, and genuinely coupled Riemann matrices.",
        "",
        "Timed rows are warmed fixed-tau evaluations over varying z vectors. "
        "Wolfram process startup and result formatting are excluded.",
        "Derivative rows are correctness-only because Wolfram 14.3 computes "
        "them here through `NumericalCalculus`ND`, not a directly comparable "
        "native `SiegelTheta` derivative call.",
        "",
        "- Generated: `%s`" % datetime.datetime.now().astimezone().isoformat(),
        "- mpmath revision: `%s`" % git_revision(),
        "- Python: `%s`" % platform.python_version(),
        "- Platform: `%s`" % platform.platform(),
        "- Decimal precisions: `%s`" % ", ".join(map(str, args.dps)),
        "- mpmath timing: median of `%d` repeats, `%d` batch per repeat"
        % (args.repeat, args.number),
        "- Wolfram timing: `RepeatedTiming` target `%g` seconds, `%d` batch"
        % (args.wolfram_target, args.wolfram_batches),
        "- Arguments per batch: up to `%d`" % len(BASE_POINTS),
        "",
        "## Summary",
        "",
    ]
    if all_ratios:
        lines.extend([
            "Across all timed rows, the geometric-mean `mpmath / Wolfram` "
            "time ratio was **%.2fx**. The range was **%.2fx to %.2fx**."
            % (geometric_mean(all_ratios), min(all_ratios), max(all_ratios)),
            "",
        ])
    if all_errors:
        lines.extend([
            "The largest scaled numerical difference was `%s`."
            % mp.nstr(max(all_errors), 6),
            "",
        ])
    for dps in args.dps:
        selected = [row for row in results if row["dps"] == dps]
        ratios = [row["mp_time"] / row["wl_time"] for row in selected
                  if row["mp_time"] is not None
                  and row["wl_time"] is not None]
        errors = [row["error"] for row in selected]
        if ratios and errors:
            lines.append(
                "- At %d dps: geometric-mean ratio **%.2fx**; largest "
                "scaled difference `%s`."
                % (dps, geometric_mean(ratios), mp.nstr(max(errors), 6))
            )
        elif errors:
            lines.append("- At %d dps: largest scaled difference `%s`."
                         % (dps, mp.nstr(max(errors), 6)))

    lines.extend([
        "",
        "## Timing",
        "",
        "| dps | Case | Genus | Characteristic | mpmath (ms/call) | Wolfram (ms/call) | mpmath / Wolfram |",
        "|---:|---|---:|---|---:|---:|---:|",
    ])
    for row in results:
        if row["mp_time"] is None or row["wl_time"] is None:
            mp_time = wl_time = ratio = "—"
        else:
            mp_time = "%.3f" % (row["mp_time"] * 1e3)
            wl_time = "%.3f" % (row["wl_time"] * 1e3)
            ratio = "%.2fx" % (row["mp_time"] / row["wl_time"])
        lines.append(
            "| %(dps)d | %(case)s | %(genus)d | %(label)s | %(mp)s | "
            "%(wl)s | %(ratio)s |"
            % {
                **row,
                "mp": mp_time,
                "wl": wl_time,
                "ratio": ratio,
            }
        )

    lines.extend([
        "",
        "A ratio above 1 means Wolfram was faster; below 1 means mpmath was "
        "faster. These are steady-state call costs, not command-line latency.",
        "",
        "## Numerical agreement",
        "",
        "Error is `abs(mpmath - Wolfram) / max(1, abs(Wolfram))`, maximised "
        "over the varying z vectors.",
        "",
        "| dps | Case | Characteristic | Maximum scaled difference | Diagonal factorisation error |",
        "|---:|---|---|---:|---:|",
    ])
    for row in results:
        factor_error = ("—" if row["factor_error"] is None else
                        "`%s`" % mp.nstr(row["factor_error"], 7))
        lines.append(
            "| %(dps)d | %(case)s | %(label)s | `%(error)s` | %(factor)s |"
            % {**row, "error": mp.nstr(row["error"], 7),
               "factor": factor_error}
        )

    lines.extend([
        "",
        "## Representative coupled-case values",
        "",
        "Values below are included to make convention or parsing errors "
        "visible rather than hiding them behind aggregate error figures.",
        "",
        "| dps | Case | Characteristic | mpmath first value | Wolfram first value |",
        "|---:|---|---|---:|---:|",
    ])
    for row in results:
        if "coupled" not in row["case"]:
            continue
        lines.append(
            "| %(dps)d | %(case)s | %(label)s | `%(mp_value)s` | "
            "`%(wl_value)s` |"
            % {**row,
               "mp_value": short_complex(row["mp_values"][0]),
               "wl_value": short_complex(row["wl_values"][0])}
        )

    if args.include_values:
        lines.extend([
            "",
            "## Full Wolfram oracle values",
            "",
            "These values are emitted for promoting selected independently "
            "generated cases into the self-contained test suite.",
            "",
            "| dps | Case | Workload | z index | Real part | Imaginary part |",
            "|---:|---|---|---:|---:|---:|",
        ])
        for row in results:
            for index, value in enumerate(row["wl_values"]):
                lines.append(
                    "| %(dps)d | %(case)s | %(label)s | %(index)d | "
                    "`%(real)s` | `%(imag)s` |"
                    % {**row, "index": index,
                       "real": mp.nstr(value.real, row["dps"]),
                       "imag": mp.nstr(value.imag, row["dps"])}
                )

    lines.extend([
        "",
        "## Cases",
        "",
        "| Case | Genus | Type |",
        "|---|---:|---|",
    ])
    for case in CASES:
        lines.append("| %s | %d | %s |" % (
            case["name"], case["genus"], case["kind"]
        ))
    lines.extend([
        "",
        "The diagonal cases independently verify that a genus-g theta value "
        "equals the product of its genus-one components, for both tested "
        "characteristics. The coupled cases have nonzero off-diagonal real "
        "and imaginary parts and therefore do not factor.",
        "",
        "Wolfram calls use `SiegelTheta[tau, z]` and "
        "`SiegelTheta[{a, b}, tau, z]`; no `EllipticTheta` calls are made.",
        "",
    ])
    if failures:
        lines.extend([
            "## Incomplete workloads",
            "",
            "A failure here does not discard successful rows from the same "
            "run.",
            "",
            "| dps | Case | Workload | Reason |",
            "|---:|---|---|---|",
        ])
        for failure in failures:
            lines.append("| %(dps)d | %(case)s | %(label)s | %(reason)s |"
                         % failure)
        lines.append("")
    return "\n".join(lines)


def parse_args():
    """Parse benchmark controls."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", nargs="+", type=int, default=[30, 100])
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--number", type=int, default=1)
    parser.add_argument("--wolfram-batches", type=int, default=1)
    parser.add_argument("--wolfram-target", type=float, default=0.2)
    parser.add_argument(
        "--wolfram-timeout", type=float, default=180,
        help="maximum seconds for each isolated Wolfram workload",
    )
    parser.add_argument(
        "--case", action="append", choices=[case["name"] for case in CASES],
        help="run only this case; repeat the option to select several",
    )
    parser.add_argument(
        "--workload", action="append",
        choices=("zero", "half", "real", "dz0", "half-dz0"),
        help="run only this workload where available; repeat to select more",
    )
    parser.add_argument(
        "--correctness-only", action="store_true",
        help="compare values without collecting steady-state timings",
    )
    parser.add_argument(
        "--strict", action="store_true",
        help="return a failing status if any workload times out or fails",
    )
    parser.add_argument(
        "--include-values", action="store_true",
        help="include full-precision Wolfram values in the Markdown report",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def main():
    """Run both implementations, validate them and write the report."""
    args = parse_args()
    results = []
    failures = []
    selected_cases = [case for case in CASES
                      if not args.case or case["name"] in args.case]
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save_progress():
        args.output.write_text(
            render_report(args, results, failures), encoding="utf-8")

    for dps in args.dps:
        for case in selected_cases:
            for label in workloads(case):
                if args.workload and label not in args.workload:
                    continue
                unused_char, derivative = workload_data(label, case["genus"])
                derivative_limit = case.get("derivative_max_dps")
                if (any(derivative) and derivative_limit is not None
                        and dps > derivative_limit):
                    failures.append({
                        "dps": dps,
                        "case": case["name"],
                        "label": label,
                        "reason": ("skipped: Wolfram contour differentiation "
                                   "is impractical above %d dps"
                                   % derivative_limit),
                    })
                    save_progress()
                    continue
                print("starting %d dps, %s, %s" % (
                    dps, case["name"], label), flush=True)
                try:
                    validation_only = args.correctness_only or any(derivative)
                    # Parse decimal strings inside the matching context so
                    # they are not rounded at the default 15 dps.
                    with mp.workdps(dps):
                        wl_time, wl_values = wolfram_row(
                            case, label, dps, args.wolfram_batches,
                            args.wolfram_target, args.wolfram_timeout,
                            validation_only,
                        )
                    mp_time, mp_values, factor_error = mpmath_row(
                        case, label, dps, args.number, args.repeat,
                        validation_only,
                    )
                    with mp.workdps(dps):
                        errors = [
                            scaled_error(actual, expected)
                            for actual, expected in zip(mp_values, wl_values)
                        ]
                        tolerance = mp.power(10, -max(5, dps - 5))
                        if max(errors) > tolerance:
                            worst = max(
                                range(len(errors)), key=errors.__getitem__)
                            raise RuntimeError(
                                "%s %s at %d dps failed validation: %s; "
                                "mpmath=%s; Wolfram=%s"
                                % (case["name"], label, dps,
                                   mp.nstr(errors[worst], 8),
                                   mp.nstr(mp_values[worst], dps),
                                   mp.nstr(wl_values[worst], dps))
                            )
                        results.append({
                            "dps": dps,
                            "case": case["name"],
                            "genus": case["genus"],
                            "label": label,
                            "mp_time": mp_time,
                            "wl_time": wl_time,
                            "mp_values": mp_values,
                            "wl_values": wl_values,
                            "error": max(errors),
                            "factor_error": factor_error,
                        })
                except (OSError, subprocess.CalledProcessError,
                        RuntimeError, ValueError) as exc:
                    reason = str(exc).replace("|", "\\|").replace("\n", " ")
                    failures.append({
                        "dps": dps,
                        "case": case["name"],
                        "label": label,
                        "reason": reason,
                    })
                    print("failed %d dps, %s, %s: %s" % (
                        dps, case["name"], label, reason), flush=True)
                else:
                    print("finished %d dps, %s, %s" % (
                        dps, case["name"], label), flush=True)
                save_progress()
    print("wrote %s (%d successful, %d incomplete)" % (
        args.output, len(results), len(failures)))
    if failures and args.strict:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
