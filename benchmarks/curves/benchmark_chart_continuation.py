"""Time chart continuation and integration against independent references.

Run unchanged and candidate sources in separate processes with --source.
Construction, reference quadrature, and warm-up are outside timed regions.
Example: python benchmark_chart_continuation.py --source src --trials 5
"""

import argparse
import json
import platform
from pathlib import Path
from statistics import median
import sys
from time import perf_counter


def cases(ctx, Curve):
    curve = Curve({(0, 2): 1, (5, 0): -4, (3, 0): 20, (1, 0): -16}, ctx=ctx)
    chart = curve.chart(
        {(0, 2): 1, (0, 0): -96, (2, 0): -200, (4, 0): -140, (6, 0): -40, (8, 0): -4},
        lambda t, w: (2+t*t, t*w, 2*t))
    end, seed = ctx.sqrt(6), ctx.sqrt(20160)
    forms = (lambda x, y: 1/y, lambda x, y: x/y)
    reference = tuple(2*ctx.quad(
        lambda t, k=k: (2+t*t)**k / ctx.sqrt((2+t*t)*((2+t*t)**2-1)*(4+t*t)),
        [0, end]) for k in (0, 1))
    for segments in (1, 32):
        path = tuple(-end+2*end*j/segments for j in range(segments+1))
        yield f"turning-point-{segments}-segments", curve, chart, forms, path, seed, reference

    for degree in (1, 3):
        terms = {(0, degree): 1, (0, 0): -1, (2, 0): -1}
        curve = Curve(terms, ctx=ctx)
        chart = curve.chart(terms, lambda t, w: (t, w, 1))
        left, right = ctx.mpf('-.8'), ctx.mpf('1.2')
        seed = ctx.root(1+left*left, degree)
        reference = (ctx.quad(lambda t: 1/ctx.root(1+t*t, degree), [left, 0, right]),)
        yield f"degree-{degree}", curve, chart, (lambda x, y: 1/y,), (left, right), seed, reference

    terms = {(0, 2): 1, (1, 0): -1}
    curve = Curve(terms, ctx=ctx)
    chart = curve.chart(terms, lambda t, w: (t, w, 1))
    path = (ctx.one, ctx.j, -ctx.one)
    reference = (2*(ctx.j**3-1)/3,)
    yield 'complex-square-root', curve, chart, (lambda x, y: y,), path, ctx.one, reference

    curve = Curve({(0, 2): 1, (1, 0): 1, (3, 0): -1}, ctx=ctx)
    chart = curve.monomial_chart(-2, -3)
    cutoff = ctx.mpf('.1')
    reference = (-2*ctx.quad(lambda t: 1/ctx.sqrt(1-t**4), [0, cutoff]),)
    yield 'infinity-chart-tail', curve, chart, (lambda x, y: 1/y,), (ctx.zero, cutoff), ctx.one, reference

    terms = {(2, 2): 1, (0, 1): -1, (0, 0): -1}
    curve = Curve(terms, ctx=ctx)
    chart = curve.chart(terms, lambda t, w: (t, w, 1))
    cutoff = ctx.mpf('.2')
    reference = (ctx.quad(lambda t: -2/(1+ctx.sqrt(1+4*t*t)), [0, cutoff]),)
    yield 'finite-degree-drop', curve, chart, (lambda x, y: y,), (ctx.zero, cutoff), -ctx.one, reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[2] / 'src')
    parser.add_argument('--digits', type=int, nargs='+', default=[20, 40, 80])
    parser.add_argument('--trials', type=int, default=5)
    parser.add_argument('--continuation-repeats', type=int, default=20)
    args = parser.parse_args()
    if args.trials < 1 or args.continuation_repeats < 1:
        parser.error('--trials and --continuation-repeats must be positive')
    sys.path.insert(0, str(args.source.resolve()))
    from mpmath import mp
    from generapy import Curve
    from generapy.curves.continuation import _continue_plane_curve_branch

    results = []
    for digits in args.digits:
        ctx = mp.clone()
        ctx.dps = digits
        for name, curve, chart, forms, path, seed, reference in cases(ctx, Curve):
            def lift():
                return _continue_plane_curve_branch(ctx, chart.curve, path, seed)

            def integrate():
                return curve.chart_integral(chart, forms, path, seed)

            lift()
            integrate()
            timings = {'continuation': [], 'integral': []}
            for unused in range(args.trials):
                for key, operation in (('continuation', lift), ('integral', integrate)):
                    repeats = args.continuation_repeats if key == 'continuation' else 1
                    start = perf_counter()
                    for repeat in range(repeats):
                        operation()
                    timings[key].append((perf_counter()-start)/repeats)
            branch, integral = lift(), integrate()
            error = max(abs(a-b)/max(ctx.one, abs(a), abs(b))
                        for a, b in zip(integral.values, reference))
            row = {
                'case': name, 'digits': digits, 'trials': args.trials,
                'continuation_repeats': args.continuation_repeats,
                'median_seconds': {key: median(values) for key, values in timings.items()},
                'samples_seconds': timings, 'segments': branch.steps,
                'refinements': branch.refinements, 'scaled_error': ctx.nstr(error, 12),
                'accurate': bool(error < ctx.mpf(10)**(-digits+3)),
            }
            results.append(row)
            print(f"{name} {digits} digits: {row['median_seconds']} error={row['scaled_error']}", file=sys.stderr, flush=True)
    print(json.dumps({'python': platform.python_version(), 'platform': platform.platform(),
                      'source': str(args.source.resolve()), 'results': results}, indent=2))


if __name__ == '__main__':
    main()
