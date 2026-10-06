"""Reproduce issue #5 on real and complex genus-two curves.

Run each source in a separate process, sequentially, without competing tests.
Each method starts with a fresh context and curve, so first-call timings
include cold Generapy caches. Construction, imports, and reference quadrature
are excluded. First-call and warm
times include both forms; warm samples reuse setup but evaluate the integral
again. Source hashes identify experimental, uncommitted implementations.
"""

import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
from statistics import median
import sys
from time import perf_counter


def real_case(ctx, Curve):
    curve = Curve({(0, 2): 1, (5, 0): -4, (3, 0): 20, (1, 0): -16}, ctx=ctx)
    chart = curve.chart(
        {(0, 2): 1, (0, 0): -96, (2, 0): -200, (4, 0): -140, (6, 0): -40, (8, 0): -4},
        lambda t, w: (2+t*t, t*w, 2*t))
    end = ctx.sqrt(6)
    point = (ctx.mpf(8), ctx.sqrt(120960))
    forms = (lambda x, y: 1/y, lambda x, y: x/y)

    def direct():
        return tuple(ctx.quad(
            lambda t, k=k: (2+t*t)**k / ctx.sqrt((2+t*t)*((2+t*t)**2-1)*(4+t*t)),
            [0, end]) for k in (0, 1))

    return {
        'direct': direct,
        'chart': lambda: curve.chart_integral(chart, forms, (0, end), 4*ctx.sqrt(6)).values,
        'finite-base-abel': lambda: curve.abel_map_kind_1(point, base_place=(2, 0), reduce=False).value,
    }


def complex_case(ctx, Curve):
    # Algebraic five-mode reconstruction from the standalone issue comment.
    u = [ctx.mpc(a, b) for a, b in [('1', '.1'), ('.7', '-.2'), ('1.2', '.15'), ('.8', '.3'), ('.6', '-.1')]]
    v = [ctx.mpc(a, b) for a, b in [('.9', '-.1'), ('1.1', '.25'), ('.5', '-.2'), ('1.3', '.1'), ('.75', '.2')]]
    xi = [ctx.mpf(a) for a in ('-.7', '-.2', '.3', '.8', '1.4')]
    q = [ctx.mpf(a) for a in ('.2', '-.3', '.4')]
    x0 = ctx.fsum(u[i]*v[i]+xi[i] for i in range(5))/5
    labels = [x0-u[i]*v[i] for i in range(5)]
    q[0] += (ctx.fprod(v)-ctx.fprod(u)-2*ctx.polyval(list(reversed(q)), x0))/2
    coefficients = [ctx.one]
    for root in labels:
        out = [ctx.zero]*(len(coefficients)+1)
        for index, value in enumerate(coefficients):
            out[index] -= root*value
            out[index+1] += value
        coefficients = out
    for i in range(3):
        for j in range(3):
            coefficients[i+j] += q[i]*q[j]
    polynomial = lambda x: ctx.polyval(list(reversed(coefficients)), x)
    y0 = (ctx.fprod(u)+ctx.fprod(v))/2
    root = ctx.sqrt(polynomial(x0))
    sign = 1 if abs(root-y0) < abs(root+y0) else -1
    sheet = lambda x: sign*ctx.sqrt(polynomial(x))
    end = x0+ctx.mpf('.025')+ctx.j*ctx.mpf('.005')
    delta = end-x0
    curve = Curve({(0, 2): ctx.one, **{(i, 0): -4*a for i, a in enumerate(coefficients) if a}}, ctx=ctx)
    base, target = (x0, 2*y0), (end, 2*sheet(end))
    forms = (lambda x, y: 1/y, lambda x, y: x/y)
    path = curve.path(base, target)

    def direct():
        return tuple(ctx.quad(
            lambda t, k=k: delta*(x0+t*delta)**k/(2*sheet(x0+t*delta)), [0, 1])
            for k in (0, 1))

    return {
        'direct': direct,
        'prebuilt-path': lambda: curve.integral(forms, path).values,
        'path-and-integral': lambda: curve.integral(forms, curve.path(base, target)).values,
        'finite-base-abel': lambda: curve.abel_map_kind_1(target, base_place=base, reduce=False).value,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[2]/'src')
    parser.add_argument('--digits', type=int, nargs='+', default=[20, 40, 80])
    parser.add_argument('--trials', type=int, default=5)
    args = parser.parse_args()
    if args.trials < 1 or any(digits < 15 for digits in args.digits):
        parser.error('use at least one trial and at least 15 digits')
    source = args.source.resolve()
    digest = hashlib.sha256()
    for path in sorted((source/'generapy').rglob('*.py')):
        digest.update(str(path.relative_to(source)).encode())
        digest.update(path.read_bytes())
    source_hash = digest.hexdigest()
    sys.path.insert(0, str(source))
    from mpmath import mp
    from generapy import Curve

    rows = []
    for digits in args.digits:
        for name, factory in (('real-branch', real_case), ('complex-regular', complex_case)):
            ctx = mp.clone()
            ctx.dps = digits
            methods = factory(ctx, Curve)
            reference = tuple(methods['direct']())
            for method in methods:
                ctx = mp.clone()
                ctx.dps = digits
                call = factory(ctx, Curve)[method]
                start = perf_counter()
                value = tuple(call())
                first_seconds = perf_counter()-start
                samples, errors = [], []
                for trial in range(args.trials+1):
                    error = max(abs(a-b)/max(ctx.one, abs(a), abs(b)) for a, b in zip(value, reference))
                    errors.append(error)
                    if trial < args.trials:
                        start = perf_counter()
                        value = tuple(call())
                        samples.append(perf_counter()-start)
                row = {
                    'case': name, 'method': method, 'digits': digits,
                    'first_seconds': first_seconds, 'median_seconds': median(samples),
                    'samples_seconds': samples, 'scaled_error': ctx.nstr(max(errors), 12),
                    'accurate': bool(max(errors) < ctx.mpf(10)**(-digits+8)),
                }
                rows.append(row)
                print(f"{name} {method} {digits}: {row['median_seconds']:.6f}s error={row['scaled_error']}", file=sys.stderr, flush=True)
    print(json.dumps({'python': platform.python_version(), 'platform': platform.platform(),
                      'mpmath': version('mpmath'), 'generapy': version('generapy'),
                      'source': str(source), 'source_sha256': source_hash, 'trials': args.trials,
                      'results': rows}, indent=2))
    if not all(row['accurate'] for row in rows):
        sys.exit(1)


if __name__ == '__main__':
    main()
