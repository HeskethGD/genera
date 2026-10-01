#!/usr/bin/env python3
"""Quickly compare four Riemann theta implementations."""

from __future__ import annotations

import argparse
import datetime
import json
import math
import os
import platform
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_REPORT = SCRIPT_DIR / "rtheta_implementation_derivative_comparison.md"
SAGE = "/usr/local/bin/sage"
WOLFRAMSCRIPT = "/usr/local/bin/wolframscript"
WOLFRAM_KERNEL = (
    "/Applications/Wolfram Engine.app/Contents/Resources/"
    "Wolfram Player.app/Contents/MacOS/WolframKernel"
)


CASES = (
    {
        "name": "g1",
        "tau": [[("0.20", "1.20")]],
        "zs": [
            [("0.11", "0.03")],
            [("-0.18", "0.06")],
        ],
    },
    {
        "name": "g2-coupled",
        "tau": [
            [("0.10", "0.90"), ("-0.20", "0.12")],
            [("-0.20", "0.12"), ("0.05", "1.10")],
        ],
        "zs": [
            [("0.11", "0.03"), ("-0.07", "0.02")],
            [("-0.18", "0.06"), ("0.13", "-0.04")],
        ],
    },
    {
        "name": "g3-coupled",
        "tau": [
            [("0.10", "1.00"), ("-0.12", "0.08"),
             ("0.05", "0.04")],
            [("-0.12", "0.08"), ("0.20", "1.15"),
             ("0.08", "0.06")],
            [("0.05", "0.04"), ("0.08", "0.06"),
             ("-0.15", "0.90")],
        ],
        "zs": [
            [("0.11", "0.03"), ("-0.07", "0.02"),
             ("0.05", "-0.01")],
            [("-0.18", "0.06"), ("0.13", "-0.04"),
             ("-0.09", "0.02")],
        ],
    },
    {
        "name": "g2-unreduced",
        "tau": [
            [("0.20", "0.15"), ("0.12", "0.03")],
            [("0.12", "0.03"), ("-0.10", "0.20")],
        ],
        "zs": [
            [("0.13", "0.07"), ("-0.21", "0.04")],
            [("2.31", "-1.14"), ("-1.72", "0.83")],
        ],
    },
    {
        "name": "g3-unreduced",
        "tau": [
            [("0.30", "0.12"), ("0.17", "0.03"),
             ("0", "0.02")],
            [("0.17", "0.03"), ("-0.20", "0.16"),
             ("0", "0.02")],
            [("0", "0.02"), ("0", "0.02"),
             ("0.10", "0.20")],
        ],
        "zs": [
            [("0.11", "0.03"), ("-0.07", "0.02"),
             ("0.05", "-0.01")],
            [("-0.18", "0.06"), ("0.13", "-0.04"),
             ("-0.09", "0.02")],
        ],
    },
)

WORKLOAD_ORDER = ("value", "d0", "d01", "jet2", "jet3")


def multiindices(genus, degree):
    """Return graded multi-indices through the requested total degree."""
    result = []
    index = [0] * genus

    def generate(position, remaining):
        if position == genus - 1:
            index[position] = remaining
            result.append(tuple(index))
            return
        for value in range(remaining + 1):
            index[position] = value
            generate(position + 1, remaining - value)

    for total in range(degree + 1):
        generate(0, total)
    return tuple(result)


def workloads(case, selected):
    """Return named derivative sets for one benchmark case."""
    genus = len(case["tau"])
    zero = (0,) * genus
    available = {
        "value": (zero,),
        "d0": ((1,) + (0,) * (genus - 1),),
        "d01": (((1, 1) + (0,) * (genus - 2),)
                 if genus > 1 else ((2,),)),
        "jet2": multiindices(genus, 2),
        "jet3": multiindices(genus, 3),
    }
    # Genus-one derivatives already have a dedicated jtheta benchmark. Keep
    # this comparison focused on the higher-genus target use case.
    return tuple((name, available[name]) for name in selected
                 if name == "value" or genus > 1)


def coordinate_indices(derivative):
    """Convert a multi-index to nbruin's repeated-coordinate convention."""
    return [i for i, order in enumerate(derivative) for unused in range(order)]


def coordinate_directions(derivative):
    """Convert a multi-index to abelfunctions direction vectors."""
    genus = len(derivative)
    return [[int(i == j) for j in range(genus)]
            for i in coordinate_indices(derivative)]


def sage_number(field, value):
    """Construct a Sage complex number from decimal strings."""
    real, imag = value
    return field(real) + field.gen() * field(imag)


def sage_worker(implementation, dps_values, repeat, selected):
    """Run inside Sage and emit machine-readable results."""
    from sage.all import ComplexField, matrix, vector

    if implementation == "nbruin":
        from riemann_theta.riemann_theta import RiemannTheta
    else:
        from abelfunctions.riemann_theta import RiemannTheta

    rows = []
    for dps in dps_values:
        if implementation == "abelfunctions" and dps > 15:
            continue
        bits = math.ceil(dps * math.log2(10))
        field = ComplexField(bits)
        for case in CASES:
            tau = matrix(field, [
                [sage_number(field, value) for value in row]
                for row in case["tau"]
            ])
            zs = [vector(field, [sage_number(field, value) for value in z])
                  for z in case["zs"]]
            theta = RiemannTheta(tau) if implementation == "nbruin" else None
            epsilon = 10.0 ** -dps
            for workload, derivatives in workloads(case, selected):
                nbruin_derivatives = [coordinate_indices(derivative)
                                      for derivative in derivatives]
                abel_derivatives = [coordinate_directions(derivative)
                                    for derivative in derivatives]

                def evaluate():
                    values = []
                    for z in zs:
                        if implementation == "nbruin":
                            if workload == "value":
                                values.append(theta(z=z))
                            else:
                                if len(derivatives) == 1:
                                    values.append(theta(
                                        z=z, derivs=nbruin_derivatives[0]))
                                else:
                                    result = theta(
                                        z=z, derivs=nbruin_derivatives)
                                    values.extend(result)
                        else:
                            for directions in abel_derivatives:
                                values.append(RiemannTheta(
                                    z, tau, epsilon=epsilon,
                                    derivs=directions))
                    return values

                values = evaluate()
                samples = []
                for unused in range(repeat):
                    start = time.perf_counter()
                    evaluate()
                    samples.append((time.perf_counter() - start) / len(zs))
                serialised = []
                for value in values:
                    if hasattr(value, "real") and callable(value.real):
                        serialised.append([str(value.real()), str(value.imag())])
                    else:
                        serialised.append([repr(value.real), repr(value.imag)])
                rows.append({
                    "implementation": implementation,
                    "dps": dps,
                    "case": case["name"],
                    "workload": workload,
                    "outputs": len(derivatives),
                    "seconds": statistics.median(samples),
                    "values": serialised,
                })
    print("RT_JSON=" + json.dumps(rows))


def mp_number(mp, value):
    """Construct an mpmath complex number from decimal strings."""
    return mp.mpc(mp.mpf(value[0]), mp.mpf(value[1]))


def run_mpmath(dps_values, repeat, selected):
    """Evaluate and time the local mpmath implementation."""
    sys.path.insert(0, str(REPO_ROOT/'src'))
    import genera
    from mpmath import mp
    from genera.riemann_theta import _rtheta_derivatives

    rows = []
    for dps in dps_values:
        with mp.workdps(dps):
            for case in CASES:
                tau = [[mp_number(mp, value) for value in row]
                       for row in case["tau"]]
                zs = [[mp_number(mp, value) for value in z]
                      for z in case["zs"]]

                for workload, derivatives in workloads(case, selected):

                    def evaluate():
                        values = []
                        for z in zs:
                            if workload in ("jet2", "jet3"):
                                values.extend(_rtheta_derivatives(
                                    mp, z, tau, None, derivatives))
                            else:
                                values.append(mp.rtheta(
                                    z, tau, derivative=derivatives[0]))
                        return values

                    values = evaluate()
                    samples = []
                    for unused in range(repeat):
                        start = time.perf_counter()
                        evaluate()
                        samples.append(
                            (time.perf_counter() - start) / len(zs))
                    seconds = statistics.median(samples) if samples else 0
                    rows.append({
                        "implementation": "mpmath",
                        "dps": dps,
                        "case": case["name"],
                        "workload": workload,
                        "outputs": len(derivatives),
                        "seconds": seconds,
                        "values": [[mp.nstr(value.real, dps + 5),
                                    mp.nstr(value.imag, dps + 5)]
                                   for value in values],
                    })
    return rows


def mpmath_oracles(dps, selected):
    """Compute higher-precision references without timing them."""
    rows = run_mpmath([dps], 0, selected)
    return {(row["case"], row["workload"]): row["values"]
            for row in rows}


def wl_real(value, digits):
    """Represent an exact decimal input in Wolfram Language."""
    return "Rationalize[%s`%d, 0]" % (value, digits)


def wl_number(value, digits):
    """Represent one complex input in Wolfram Language."""
    real, imag = value
    return "(%s + I %s)" % (
        wl_real(real, digits), wl_real(imag, digits))


def wl_vector(values, digits):
    """Represent a vector in Wolfram Language."""
    return "{%s}" % ",".join(wl_number(value, digits) for value in values)


def wl_matrix(rows, digits):
    """Represent a matrix in Wolfram Language."""
    return "{%s}" % ",".join(wl_vector(row, digits) for row in rows)


def wolfram_code(dps_values, target, row_timeout, oracle_dps):
    """Build one bounded Wolfram process containing every workload."""
    expressions = []
    for dps in dps_values:
        for case in CASES:
            tau = wl_matrix(case["tau"], dps)
            zs = "{%s}" % ",".join(wl_vector(z, dps) for z in case["zs"])
            oracle_tau = wl_matrix(case["tau"], oracle_dps)
            oracle_zs = "{%s}" % ",".join(
                wl_vector(z, oracle_dps) for z in case["zs"])
            expressions.append(r"""
TimeConstrained[
 Module[{tau=N[%(tau)s,%(dps)d],zs=N[%(zs)s,%(dps)d],f,values,
         oracleValues,elapsed},
  f=Function[z,SiegelTheta[tau,z]];
  values=N[f/@zs,%(dps)d]; f/@zs;
  elapsed=First@RepeatedTiming[Scan[f,zs],%(target)s]/Length[zs];
  oracleValues=N[
   SiegelTheta[N[%(oracle_tau)s,%(oracle_dps)d],#]&/@
    N[%(oracle_zs)s,%(oracle_dps)d],%(oracle_dps)d];
  {"ok","wolfram",%(dps)d,%(name)s,
   ToString[N[elapsed,17],InputForm],
   ({ToString[Re[#],InputForm],ToString[Im[#],InputForm]}&/@values),
   ({ToString[Re[#],InputForm],ToString[Im[#],InputForm]}&/@oracleValues)}],
 %(timeout)s,{"timeout","wolfram",%(dps)d,%(name)s}]
""" % {
                "tau": tau,
                "zs": zs,
                "dps": dps,
                "target": target,
                "oracle_tau": oracle_tau,
                "oracle_zs": oracle_zs,
                "oracle_dps": oracle_dps,
                "timeout": row_timeout,
                "name": json.dumps(case["name"]),
            })
    return "Print[ExportString[{%s},\"RawJSON\"]]" % ",".join(expressions)


def run_wolfram(dps_values, target, row_timeout, process_timeout,
                oracle_dps):
    """Run all Wolfram rows once under a hard wall-clock limit."""
    environment = os.environ.copy()
    environment["WolframKernel"] = WOLFRAM_KERNEL
    process = subprocess.run(
        [WOLFRAMSCRIPT, "-code",
         wolfram_code(dps_values, target, row_timeout, oracle_dps)],
        cwd="/tmp", env=environment, text=True, capture_output=True,
        check=True, timeout=process_timeout,
    )
    try:
        start = process.stdout.index("[")
        raw_rows, unused_end = json.JSONDecoder().raw_decode(
            process.stdout[start:])
    except (ValueError, json.JSONDecodeError) as exc:
        raise ValueError(
            "unexpected Wolfram output: stdout=%r stderr=%r"
            % (process.stdout, process.stderr)) from exc
    rows = []
    oracles = {}
    failures = []
    for raw in raw_rows:
        if raw[0] != "ok":
            failures.append("%s dps %s timed out" % (raw[2], raw[3]))
            continue
        unused, implementation, dps, case, seconds, values, oracle = raw
        rows.append({
            "implementation": implementation,
            "dps": dps,
            "case": case,
            "workload": "value",
            "outputs": 1,
            "seconds": float(seconds.replace("*^", "e").split("`")[0]),
            "values": values,
        })
        oracles[(case, "value")] = oracle
    return rows, oracles, failures


def parse_json_output(output):
    """Extract the worker payload while ignoring Sage startup messages."""
    for line in reversed(output.splitlines()):
        if line.startswith("RT_JSON="):
            return json.loads(line[len("RT_JSON="):])
    raise RuntimeError("Sage worker did not produce an RT_JSON result")


def run_sage(implementation, dps_values, repeat, timeout, dot_sage,
             selected):
    """Run one Sage implementation outside the repository working tree."""
    environment = os.environ.copy()
    environment["DOT_SAGE"] = str(dot_sage)
    command = [SAGE, "-python", str(Path(__file__).resolve()),
               "--sage-worker", implementation,
               "--dps", *map(str, dps_values), "--repeat", str(repeat),
               "--workloads", *selected]
    process = subprocess.run(
        command, cwd="/tmp", env=environment, text=True,
        capture_output=True, check=True, timeout=timeout,
    )
    return parse_json_output(process.stdout)


def parse_decimal(mp, value):
    """Parse Sage or Wolfram decimal output with precision marks removed."""
    value = re.sub(r"`[0-9.]*", "", value).replace("*^", "e")
    return mp.mpf(value)


def achieved_digits(row, oracle):
    """Return the worst-case decimal accuracy for a result row."""
    import genera
    from mpmath import mp

    with mp.workdps(70):
        errors = []
        for actual, expected in zip(row["values"], oracle):
            actual = mp.mpc(parse_decimal(mp, actual[0]),
                            parse_decimal(mp, actual[1]))
            expected = mp.mpc(parse_decimal(mp, expected[0]),
                              parse_decimal(mp, expected[1]))
            errors.append(abs(actual - expected) / max(1, abs(expected)))
        error = max(errors)
        return math.inf if not error else float(-mp.log10(error))


def render_report(args, rows, oracles, failures):
    """Render the combined speed and accuracy table."""
    lines = [
        "# Riemann theta derivative and jet comparison",
        "",
        "This bounded benchmark compares selected zero-characteristic values,",
        "partial derivatives, and complete jets at two varying arguments for",
        "each higher-genus period matrix. Derivative workloads omit genus one;",
        "it appears only when the value baseline is explicitly selected.",
        "Reported timings are warmed steady-state costs and exclude "
        "interpreter",
        "startup. Lower time and higher achieved digits are better.",
        "",
        "- Generated: `%s`" % datetime.datetime.now().astimezone().isoformat(),
        "- Python: `%s`" % platform.python_version(),
        "- Precisions: `%s` decimal digits" % ", ".join(map(str, args.dps)),
        "- Workloads: `%s`" % ", ".join(args.workloads),
        "- Timing samples: median of `%d` batches" % args.repeat,
        "- Reference precision: `%d` decimal digits" % args.oracle_dps,
        "",
        "`d0` is the first derivative in coordinate zero; `d01` is the mixed",
        "second derivative in coordinates zero and one. `jet2` and `jet3`",
        "contain the value and every distinct partial derivative through total",
        "order two or three. Time is per argument for the complete workload.",
        "",
        "| dps | Case | Workload | Outputs | Implementation | Time "
        "(ms/argument) | Achieved digits |",
        "|---:|---|---|---:|---|---:|---:|",
    ]
    order = {"mpmath": 0, "nbruin": 1, "abelfunctions": 2, "wolfram": 3}
    workload_order = {name: i for i, name in enumerate(WORKLOAD_ORDER)}
    for row in sorted(rows, key=lambda item: (
            item["dps"], item["case"], workload_order[item["workload"]],
            order[item["implementation"]])):
        digits = achieved_digits(
            row, oracles[(row["case"], row["workload"])])
        shown_digits = ">= %d" % row["dps"] if digits >= row["dps"] \
            else "%.1f" % max(0, digits)
        lines.append("| %d | %s | %s | %d | %s | %.3f | %s |" % (
            row["dps"], row["case"], row["workload"], row["outputs"],
            row["implementation"], row["seconds"] * 1000, shown_digits))

    timings = {(row["implementation"], row["dps"], row["case"],
                row["workload"]): row["seconds"] for row in rows}
    if any(row["implementation"] == "nbruin" for row in rows):
        lines.extend([
            "",
            "## Geometric-mean speed comparison",
            "",
            "The ratio is nbruin RiemannTheta time divided by mpmath time;",
            "values above one mean mpmath is faster.",
            "",
            "| dps | Cases | Matched workloads | nbruin / mpmath |",
            "|---:|---|---:|---:|",
        ])
        for dps in args.dps:
            for label, case_filter in (
                    ("all", lambda name: True),
                    ("coupled", lambda name: "coupled" in name),
                    ("unreduced", lambda name: "unreduced" in name)):
                ratios = []
                for row in rows:
                    if (row["implementation"] != "mpmath"
                            or row["dps"] != dps
                            or not case_filter(row["case"])):
                        continue
                    key = ("nbruin", dps, row["case"], row["workload"])
                    if key in timings:
                        ratios.append(timings[key] / row["seconds"])
                if ratios:
                    ratio = math.exp(statistics.fmean(
                        math.log(value) for value in ratios))
                    lines.append("| %d | %s | %d | %.3fx |" % (
                        dps, label, len(ratios), ratio))

    if args.skip_wolfram:
        value_reference = (
            "All references are this branch's mpmath implementation at the "
            "stated higher precision.")
    else:
        value_reference = (
            "Value references are Wolfram `SiegelTheta` at the stated oracle "
            "precision. Derivative references are this branch's mpmath "
            "implementation at that higher precision; independent Wolfram "
            "derivative checks remain in the separate correctness benchmark.")
    lines.extend([
        "",
        "Accuracy is `-log10(abs(value-reference)/max(1,abs(reference)))`,",
        "minimised over every output at both arguments and capped at the",
        "requested precision. " + value_reference,
        "",
        "`mpmath` uses its warmed precision-aware period-matrix cache. The",
        "nbruin timing excludes `RiemannTheta(tau)` construction and uses its",
        "native multi-derivative call for jets.",
        "Abelfunctions has no corresponding reusable period-matrix object, so",
        "its public function's preprocessing remains part of each call. Its",
        "API accepts one derivative at a time, so a jet is assembled from",
        "separate public calls. It is only run at 15 dps because its numerical",
        "core returns machine doubles.",
        "",
        "The default run includes favourable and strongly unreduced matrices",
        "and does not launch Wolfram. Wolfram is",
        "invoked only if `value` is selected and is never asked for numerical",
        "derivatives. No characteristics or 100/120-digit Wolfram workloads",
        "are included.",
    ])
    if failures:
        lines.extend(["", "## Incomplete workloads", ""])
        lines.extend("- " + failure for failure in failures)
    lines.append("")
    return "\n".join(lines)


def parse_args():
    """Parse benchmark controls."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", nargs="+", type=int, default=[15, 30])
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--workloads", nargs="+", choices=WORKLOAD_ORDER,
                        default=list(WORKLOAD_ORDER[1:]))
    parser.add_argument("--wolfram-target", type=float, default=0.05)
    parser.add_argument("--wolfram-row-timeout", type=float, default=5)
    parser.add_argument("--wolfram-timeout", type=float, default=45)
    parser.add_argument("--sage-timeout", type=float, default=45)
    parser.add_argument("--mpmath-only", action="store_true",
                        help="skip Sage and Wolfram subprocesses")
    parser.add_argument("--sage-implementations", nargs="+",
                        choices=("nbruin", "abelfunctions"),
                        default=("nbruin", "abelfunctions"),
                        help="select Sage implementations to run")
    parser.add_argument("--skip-wolfram", action="store_true",
                        help="use mpmath value oracles without running Wolfram")
    parser.add_argument("--oracle-dps", type=int, default=50)
    parser.add_argument("--dot-sage", type=Path,
                        default=SCRIPT_DIR / ".sage-cache")
    parser.add_argument("--output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--sage-worker", choices=("nbruin", "abelfunctions"))
    return parser.parse_args()


def main():
    """Run the bounded comparison or an internal Sage worker."""
    args = parse_args()
    if args.sage_worker:
        sage_worker(args.sage_worker, args.dps, args.repeat, args.workloads)
        return
    args.dot_sage.mkdir(parents=True, exist_ok=True)
    failures = []
    rows = run_mpmath(args.dps, args.repeat, args.workloads)
    derivative_oracles = mpmath_oracles(args.oracle_dps, args.workloads)
    if not args.mpmath_only:
        for implementation in args.sage_implementations:
            print("running %s" % implementation, flush=True)
            try:
                rows.extend(run_sage(
                    implementation, args.dps, args.repeat,
                    args.sage_timeout, args.dot_sage, args.workloads))
            except (OSError, subprocess.CalledProcessError,
                    subprocess.TimeoutExpired, RuntimeError) as exc:
                failures.append("%s: %s" % (implementation, exc))
    oracles = derivative_oracles
    if ("value" in args.workloads and not args.mpmath_only
            and not args.skip_wolfram):
        print("running Wolfram", flush=True)
        try:
            wolfram_rows, value_oracles, wolfram_failures = run_wolfram(
                args.dps, args.wolfram_target, args.wolfram_row_timeout,
                args.wolfram_timeout, args.oracle_dps)
            rows.extend(wolfram_rows)
            oracles.update(value_oracles)
            failures.extend(wolfram_failures)
        except (OSError, subprocess.CalledProcessError,
                subprocess.TimeoutExpired, ValueError) as exc:
            raise SystemExit(
                "Wolfram failed before producing the common oracle: %s" % exc)
    args.output.write_text(render_report(args, rows, oracles, failures),
                           encoding="utf-8")
    print("wrote %s" % args.output)


if __name__ == "__main__":
    main()
