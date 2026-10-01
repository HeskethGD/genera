"""One fresh process, one curve-to-full-periods workload. JSON stdin/stdout."""
import contextlib
import io
import json
from pathlib import Path
import sys
import time


def encode(matrix):
    rows = matrix.rows() if callable(getattr(matrix, 'rows', None)) else matrix.tolist()
    out = []
    for row in rows:
        values = []
        for z in row:
            real = z.real() if callable(z.real) else z.real
            imag = z.imag() if callable(z.imag) else z.imag
            values.append([str(real), str(imag)])
        out.append(values)
    return out


def run(request):
    case, engine, bits = request['case'], request['engine'], request['bits']
    if engine == 'sage':
        from sage.all import Curve, PolynomialRing, QQ, version
        import mpmath  # Sage's bundled mpmath, not genera
        ring = PolynomialRing(QQ, names=('x','y'))
        x, y = ring.gens()
        def poly(terms):
            return sum(QQ(c)*x**i*y**j for i,j,c in terms)
        # Explicit warmup is outside the measured workload in every process.
        Curve(y*y-x**3+x).riemann_surface(prec=53).period_matrix()
        tick = time.perf_counter()
        polynomial = poly(case['terms'])
        surface = Curve(polynomial).riemann_surface(prec=bits)
        periods = surface.period_matrix()
        seconds = time.perf_counter()-tick
        basis = [[[int(i),int(j),str(c)] for (i,j),c in h.dict().items()]
                 for h in surface.cohomology_basis()]
        return dict(status='ok', seconds=seconds, periods=encode(periods), basis=basis,
                    genus=int(surface.genus), actual_engine='sage', marking='sage',
                    version=version(), mpmath_path=mpmath.__file__)
    # Only the genera worker adds the checkout to its import path.
    sys.path.insert(0, str(Path(request['repo'])/'src'))
    import genera
    from genera import algebraic_curve
    from mpmath import mp
    mp.prec = bits
    genera_path = Path(genera.__file__).resolve()
    if not genera_path.is_relative_to(Path(request['repo']).resolve()):
        raise RuntimeError('genera worker imported the wrong checkout')
    with mp.workprec(53):
        algebraic_curve((0,-1,0,1)).periods_kind_1()
    from genera.curves import _stages
    for stage in vars(_stages).values():
        if callable(stage) and hasattr(stage,'cache_clear'):
            stage.cache_clear()
    def poly(terms, x, y):
        return mp.fsum(mp.mpf(c)*x**i*y**j for i,j,c in terms)
    tick = time.perf_counter()
    terms = {(i,j):mp.mpf(c) for i,j,c in case['terms']}
    curve = algebraic_curve(terms)
    supplied = case.get('numerators')
    forms = None
    if supplied is not None:
        derivative = [[i,j-1,str(mp.mpf(c)*j)] for i,j,c in case['terms'] if j]
        forms = tuple((lambda x,y,h=h: poly(h,x,y)/poly(derivative,x,y)) for h in supplied)
    data = curve.periods_kind_1(forms)
    seconds = time.perf_counter()-tick
    full = mp.matrix([list(a)+list(b) for a,b in zip((2*data.omega).tolist(),(2*data.omega_prime).tolist())])
    if supplied is not None:
        basis = supplied
    elif data.engine == 'hyperelliptic':
        # Specialized integrals are x^k dx/z, z=y+B/(2A). Here A=1,
        # so their exact h/F_y numerators are 2*x^k, even for linear-y input.
        basis = [[[k,0,'2']] for k in range(data.genus)]
    else:
        basis = [[[f.numerator[0],f.numerator[1],'1']] for f in data.differentials]
    return dict(status='ok', seconds=seconds, periods=encode(full), basis=basis,
                genus=data.genus, actual_engine=data.engine, marking=data.marking,
                version=genera.__version__, genera_path=str(genera_path),
                validation=curve.validate(data).passed,
                symmetry_residual=str(data.symmetry_residual),
                minimum_imaginary_eigenvalue=str(min(data.imaginary_eigenvalues)))


if __name__ == '__main__':
    request = json.load(sys.stdin)
    diagnostics = io.StringIO()
    try:
        with contextlib.redirect_stdout(diagnostics):
            result = run(request)
    except Exception as exc:
        import traceback
        result = dict(status='error', error_type=type(exc).__name__, error=str(exc),
                      traceback=traceback.format_exc())
    result['working_bits'] = request['bits']
    result['worker_python'] = sys.version
    result['worker_executable'] = sys.executable
    result['diagnostics'] = diagnostics.getvalue()[-4000:]
    print(json.dumps(result))
