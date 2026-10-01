#!/usr/bin/env python3
"""Benchmark genus-one ``rtheta`` against the specialised ``jtheta``.

Run this script from anywhere in the checkout. It imports the local mpmath
source tree and writes a Markdown report next to this file by default.
"""

from __future__ import annotations

import argparse
import datetime
import math
import platform
import statistics
import subprocess
import sys
import timeit
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_REPORT = SCRIPT_DIR / "rtheta_vs_jtheta.md"
sys.path.insert(0, str(REPO_ROOT/'src'))

import genera
from mpmath import mp  # noqa: E402


TAU_CASES = (
    ("comfortable", "0.20", "1.20", "rapid direct convergence"),
    ("generic", "0.35", "0.75", "moderate direct convergence"),
    ("near-boundary", "0.20", "0.30",
     "transformation overhead outweighs the small predicted saving"),
    ("strongly-unreduced", "0.17", "0.003",
     "benefits substantially from modular reduction"),
)

W_VALUES = (
    ("0.00", "0.00"),
    ("0.11", "0.03"),
    ("-0.23", "0.07"),
    ("0.37", "-0.09"),
    ("-0.41", "-0.12"),
    ("0.52", "0.16"),
)

WORKLOADS = (
    ("theta1", 1, 0, ("0.5", "0.5"), -1),
    ("theta2", 2, 0, ("0.5", "0.0"), 1),
    ("theta3", 3, 0, ("0.0", "0.0"), 1),
    ("theta4", 4, 0, ("0.0", "0.5"), 1),
    ("theta3 d1", 3, 1, ("0.0", "0.0"), 1),
    ("theta3 d2", 3, 2, ("0.0", "0.0"), 1),
)

SWEEP_POINT_COUNT = 50


def git_revision():
    """Return the checked-out revision and whether the tree is dirty."""
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


def median_time(function, number, repeat, calls_per_batch):
    """Return the median wall time per mathematical function call."""
    samples = timeit.repeat(function, number=number, repeat=repeat)
    return statistics.median(samples) / (number * calls_per_batch)


def scaled_error(actual, expected):
    """Return an absolute error scaled by max(1, abs(expected))."""
    return abs(actual - expected) / max(mp.one, abs(expected))


def make_inputs(tau_case, workload):
    """Construct exact decimal benchmark inputs at the active precision."""
    unused_name, tau_re, tau_im, unused_note = tau_case
    unused_label, theta_index, derivative, characteristic, sign = workload
    tau = mp.mpc(tau_re, tau_im)
    q = mp.expjpi(tau)
    points = tuple(mp.mpc(real, imag) for real, imag in W_VALUES)
    a, b = (mp.mpf(value) for value in characteristic)
    return tau, q, points, ([a], [b]), theta_index, derivative, sign


def make_sweep_points(count):
    """Construct a deterministic plot-like sequence of angular arguments."""
    golden_ratio_conjugate = (mp.sqrt(5) - 1) / 2
    return tuple(
        mp.mpf("0.1") + mp.mpf("0.8") * (index + mp.mpf("0.5")) / count
        + mp.j * mp.mpf("0.1") * mp.frac(
            (index + mp.mpf("0.5")) * golden_ratio_conjugate
        )
        for index in range(count)
    )


def reference_value(w, q, theta_index, derivative, sign):
    """Evaluate the jtheta equivalent in the rtheta z convention."""
    return (sign * mp.pi ** derivative
            * mp.jtheta(theta_index, w, q, derivative=derivative))


def collect_row(dps, tau_case, workload, number, repeat):
    """Validate and time one precision, tau and workload combination."""
    with mp.workdps(dps):
        tau, q, points, characteristic, theta_index, derivative, sign = (
            make_inputs(tau_case, workload)
        )

        rtheta_values = [
            mp.rtheta([w / mp.pi], [[tau]], characteristic,
                      derivative=derivative)
            for w in points
        ]
        # Construct q above the measured precision for the accuracy oracle.
        # When |q| is close to one, rounding q at the target precision can
        # noticeably perturb the tau represented by log(q)/(pi*i).
        with mp.extradps(30):
            (unused_tau, reference_q, reference_points,
             unused_characteristic, reference_index,
             reference_derivative, reference_sign) = make_inputs(
                tau_case, workload)
            jtheta_values = [
                reference_value(w, reference_q, reference_index,
                                reference_derivative, reference_sign)
                for w in reference_points
            ]
        errors = [scaled_error(actual, expected)
                  for actual, expected in zip(rtheta_values, jtheta_values)]
        tolerance = mp.power(10, -max(5, dps - 5))
        if max(errors) > tolerance:
            raise RuntimeError(
                "%s at %d dps failed accuracy validation: %s"
                % (workload[0], dps, mp.nstr(max(errors), 8))
            )

        def rtheta_batch():
            for w in points:
                mp.rtheta([w / mp.pi], [[tau]], characteristic,
                          derivative=derivative)

        def cold_rtheta_batch():
            for w in points:
                mp._rtheta_tau_data.cache_clear()
                mp._rtheta_reduction_data.cache_clear()
                mp._rtheta_derivative_radius.cache_clear()
                mp.rtheta([w / mp.pi], [[tau]], characteristic,
                          derivative=derivative)

        def jtheta_batch():
            for w in points:
                reference_value(w, q, theta_index, derivative, sign)

        # Ensure the fixed-tau path starts with populated preprocessing data.
        mp._rtheta_tau_data.cache_clear()
        mp._rtheta_reduction_data.cache_clear()
        mp._rtheta_derivative_radius.cache_clear()
        rtheta_batch()
        warm_seconds = median_time(rtheta_batch, number, repeat, len(points))
        cold_seconds = median_time(
            cold_rtheta_batch, number, repeat, len(points)
        )
        jtheta_seconds = median_time(jtheta_batch, number, repeat, len(points))

        return {
            "dps": dps,
            "case": tau_case[0],
            "tau": tau,
            "note": tau_case[3],
            "workload": workload[0],
            "warm": warm_seconds,
            "cold": cold_seconds,
            "jtheta": jtheta_seconds,
            "error": max(errors),
        }


def collect_sweep_row(dps, tau_case, workload, count, number, repeat):
    """Validate and time a fixed-tau sweep over preconstructed arguments."""
    with mp.workdps(dps):
        (tau, q, unused_points, characteristic, theta_index, derivative,
         sign) = make_inputs(tau_case, workload)
        points = make_sweep_points(count)
        tau_matrix = [[tau]]
        rtheta_points = tuple([w / mp.pi] for w in points)

        rtheta_values = [
            mp.rtheta(z, tau_matrix, characteristic, derivative=derivative)
            for z in rtheta_points
        ]
        with mp.extradps(30):
            (unused_tau, reference_q, unused_reference_points,
             unused_characteristic, reference_index, reference_derivative,
             reference_sign) = make_inputs(tau_case, workload)
            reference_points = make_sweep_points(count)
            jtheta_values = [
                reference_value(w, reference_q, reference_index,
                                reference_derivative, reference_sign)
                for w in reference_points
            ]
        errors = [scaled_error(actual, expected)
                  for actual, expected in zip(rtheta_values, jtheta_values)]
        tolerance = mp.power(10, -max(5, dps - 5))
        if max(errors) > tolerance:
            raise RuntimeError(
                "%s sweep at %d dps failed accuracy validation: %s"
                % (workload[0], dps, mp.nstr(max(errors), 8))
            )

        def rtheta_batch():
            for z in rtheta_points:
                mp.rtheta(z, tau_matrix, characteristic,
                          derivative=derivative)

        def jtheta_batch():
            for w in points:
                reference_value(w, q, theta_index, derivative, sign)

        mp._rtheta_tau_data.cache_clear()
        mp._rtheta_reduction_data.cache_clear()
        mp._rtheta_derivative_radius.cache_clear()
        rtheta_batch()
        rtheta_seconds = median_time(
            rtheta_batch, number, repeat, len(points)
        )
        jtheta_seconds = median_time(
            jtheta_batch, number, repeat, len(points)
        )

        return {
            "dps": dps,
            "case": tau_case[0],
            "tau": tau,
            "workload": workload[0],
            "derivative": derivative,
            "rtheta": rtheta_seconds,
            "jtheta": jtheta_seconds,
            "error": max(errors),
        }


def format_time(seconds):
    """Format a per-call duration in microseconds."""
    return f"{seconds * 1e6:.2f}"


def geometric_mean(values):
    """Return a geometric mean without depending on Python version details."""
    return math.exp(statistics.mean(math.log(value) for value in values))


def make_report(args, rows, sweep_rows):
    """Render benchmark metadata and results as Markdown."""
    lines = [
        "# Riemann theta genus-one baseline: `rtheta` vs `jtheta`",
        "",
        "This compares the generic Riemann theta evaluator, including its "
        "cost-selected Siegel reduction, with mpmath's specialised genus-one "
        "`jtheta` implementation. It is not expected to outperform `jtheta` "
        "on ordinary genus-one inputs.",
        "",
        "Each timing evaluates %d deterministic complex arguments. Times are "
        "median wall-clock microseconds per individual function call; lower "
        "is better." % len(W_VALUES),
        "",
        "- Generated: `%s`" % datetime.datetime.now().astimezone().isoformat(),
        "- mpmath revision: `%s`" % git_revision(),
        "- Python: `%s`" % platform.python_version(),
        "- Platform: `%s`" % platform.platform(),
        "- Decimal precisions: `%s`" % ", ".join(map(str, args.dps)),
        "- Timing repeats: `%d` (median reported)" % args.repeat,
        "- Batches per repeat: `%d`" % args.number,
        "- Arguments per batch: `%d`" % len(W_VALUES),
        "- Fixed-tau sweep arguments per batch: `%d`" % args.sweep_points,
        "- Fixed-tau sweep repeats: `%d` (median reported)"
        % args.sweep_repeat,
        "",
        "## Summary",
        "",
    ]

    for dps in args.dps:
        selected = [row for row in rows if row["dps"] == dps]
        warm_ratios = [row["warm"] / row["jtheta"] for row in selected]
        cold_ratios = [row["cold"] / row["jtheta"] for row in selected]
        largest_error = max(row["error"] for row in selected)
        lines.append(
            "- At %d dps, the geometric-mean `rtheta / jtheta` ratio was "
            "**%.2fx warm** and **%.2fx cold**. The largest scaled error "
            "was `%s`."
            % (dps, geometric_mean(warm_ratios), geometric_mean(cold_ratios),
               mp.nstr(largest_error, 5))
        )

    lines.extend([
        "",
        "## Timing results",
        "",
        "| dps | Case | tau | Workload | rtheta warm (us) | rtheta cold (us) | jtheta (us) | Warm / jtheta | Cold / jtheta |",
        "|---:|---|---:|---|---:|---:|---:|---:|---:|",
    ])
    for row in rows:
        lines.append(
            "| %(dps)d | %(case)s | `%(tau)s` | %(workload)s | %(warm)s | "
            "%(cold)s | %(jtheta)s | %(warm_ratio).2fx | "
            "%(cold_ratio).2fx |"
            % {
                **row,
                "tau": mp.nstr(row["tau"], 8),
                "warm": format_time(row["warm"]),
                "cold": format_time(row["cold"]),
                "jtheta": format_time(row["jtheta"]),
                "warm_ratio": row["warm"] / row["jtheta"],
                "cold_ratio": row["cold"] / row["jtheta"],
            }
        )
    for dps in args.dps:
        selected = [row for row in rows if row["dps"] == dps]
        warm = geometric_mean([row["warm"] for row in selected])
        cold = geometric_mean([row["cold"] for row in selected])
        jtheta = geometric_mean([row["jtheta"] for row in selected])
        lines.append(
            "| %d | **geometric mean** | — | all workloads | %s | %s | "
            "%s | **%.2fx** | **%.2fx** |"
            % (dps, format_time(warm), format_time(cold),
               format_time(jtheta), warm / jtheta, cold / jtheta)
        )

    lines.extend([
        "",
        "The geometric-mean rows aggregate every case and workload at each "
        "precision.",
        "",
        "A ratio above 1 means `jtheta` was faster. The cold measurement "
        "includes clearing and rebuilding `rtheta`'s precision-aware "
        "tau-preprocessing and derivative-radius caches before every call. "
        "The warm measurement uses one fixed tau while varying z.",
        "",
        "## Numerical agreement",
        "",
        "The reported error is `abs(rtheta - converted_jtheta) / "
        "max(1, abs(converted_jtheta))`. It is measured before timing.",
        "",
        "| dps | Case | Workload | Maximum scaled error |",
        "|---:|---|---|---:|",
    ])
    for row in rows:
        lines.append(
            "| %(dps)d | %(case)s | %(workload)s | `%(error)s` |"
            % {**row, "error": mp.nstr(row["error"], 8)}
        )

    lines.extend([
        "",
        "## Cases and conventions",
        "",
        "| Case | tau | Purpose |",
        "|---|---:|---|",
    ])
    for name, real, imag, note in TAU_CASES:
        lines.append("| %s | `%s + %s i` | %s |" % (name, real, imag, note))

    lines.extend([
        "",
        "The six angular `jtheta` arguments are:",
        "",
        "```text",
        ", ".join("%s%s%si" % (
            real, "+" if not imag.startswith("-") else "", imag
        ) for real, imag in W_VALUES),
        "```",
        "",
        "They are converted using `z = w/pi` and `q = exp(pi*i*tau)`. "
        "The four standard half-characteristics use the agreed sign for "
        "theta1. A derivative of order d is compared using "
        "`d^d rtheta/dz^d = pi^d d^d jtheta/dw^d`.",
        "",
        "The near-boundary case checks that the cost gate retains direct "
        "summation when a modular transform would save too little work. The "
        "strongly-unreduced case checks the selected transformation path. "
        "`jtheta` applies its own specialised genus-one transformations.",
        "",
        "## Fixed-tau, varying-argument sweeps",
        "",
        "This plot-like workload evaluates a deterministic sequence of %d "
        "arguments while keeping tau and the characteristic fixed. Both "
        "implementations receive preconstructed arguments, and `rtheta`'s "
        "precision-aware tau and derivative-radius caches are populated "
        "before timing."
        % args.sweep_points,
        "",
    ])
    if not sweep_rows:
        lines.extend(["Sweeps were skipped for this run.", ""])
        return "\n".join(lines)

    for dps in args.dps:
        selected = [row for row in sweep_rows if row["dps"] == dps]
        value_ratios = [
            row["rtheta"] / row["jtheta"]
            for row in selected if row["derivative"] == 0
        ]
        derivative_ratios = [
            row["rtheta"] / row["jtheta"]
            for row in selected if row["derivative"] != 0
        ]
        lines.append(
            "- At %d dps, the geometric-mean `rtheta / jtheta` ratio was "
            "**%.2fx for values** and **%.2fx for derivatives**."
            % (dps, geometric_mean(value_ratios),
               geometric_mean(derivative_ratios))
        )

    lines.extend([
        "",
        "| dps | Case | Workload | rtheta (us) | jtheta (us) | rtheta / jtheta | Maximum scaled error |",
        "|---:|---|---|---:|---:|---:|---:|",
    ])
    for row in sweep_rows:
        lines.append(
            "| %(dps)d | %(case)s | %(workload)s | %(rtheta)s | "
            "%(jtheta)s | %(ratio).2fx | `%(error)s` |"
            % {
                **row,
                "rtheta": format_time(row["rtheta"]),
                "jtheta": format_time(row["jtheta"]),
                "ratio": row["rtheta"] / row["jtheta"],
                "error": mp.nstr(row["error"], 8),
            }
        )
    for dps in args.dps:
        selected = [row for row in sweep_rows if row["dps"] == dps]
        for label, derivative in (("values", False), ("derivatives", True)):
            subset = [
                row for row in selected
                if bool(row["derivative"]) == derivative
            ]
            rtheta = geometric_mean([row["rtheta"] for row in subset])
            jtheta = geometric_mean([row["jtheta"] for row in subset])
            largest_error = max(row["error"] for row in subset)
            lines.append(
                "| %d | **geometric mean** | %s | %s | %s | **%.2fx** | "
                "`%s` |"
                % (dps, label, format_time(rtheta), format_time(jtheta),
                   rtheta / jtheta, mp.nstr(largest_error, 8))
            )

    lines.extend([
        "",
        "Sweep geometric means are separated into values and derivatives "
        "at each precision.",
        "",
        "The sweep arguments follow the same deterministic distribution as "
        "the existing `benchmark_jtheta_boundaries_pyperf.py` benchmark: "
        "their real parts cover 0.1 to 0.9 and their small imaginary parts "
        "are distributed using the golden-ratio conjugate. This represents "
        "repeated evaluation along a plot or trajectory without adding "
        "vector-valued semantics to either function.",
        "",
    ])
    return "\n".join(lines)


def parse_args():
    """Parse command-line benchmark controls."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", nargs="+", type=int, default=[15, 50, 100])
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--number", type=int, default=1)
    parser.add_argument("--sweep-points", type=int, default=SWEEP_POINT_COUNT)
    parser.add_argument("--sweep-repeat", type=int, default=3)
    parser.add_argument("--skip-sweeps", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def main():
    """Run all cases and write the Markdown report."""
    args = parse_args()
    rows = []
    for dps in args.dps:
        for tau_case in TAU_CASES:
            for workload in WORKLOADS:
                rows.append(collect_row(
                    dps, tau_case, workload, args.number, args.repeat
                ))
                print("finished %d dps, %s, %s"
                      % (dps, tau_case[0], workload[0]), flush=True)
    sweep_rows = []
    if not args.skip_sweeps:
        for dps in args.dps:
            for tau_case in TAU_CASES:
                for workload in WORKLOADS:
                    sweep_rows.append(collect_sweep_row(
                        dps, tau_case, workload, args.sweep_points,
                        args.number, args.sweep_repeat
                    ))
                    print("finished sweep %d dps, %s, %s"
                          % (dps, tau_case[0], workload[0]), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        make_report(args, rows, sweep_rows), encoding="utf-8"
    )
    print("wrote %s" % args.output)


if __name__ == "__main__":
    main()
