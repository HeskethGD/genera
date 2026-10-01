"""Exact differential-basis alignment and checked period-lattice comparison."""
from fractions import Fraction as Q


def basis_change(source, target):
    """Exact rational C with target_i = sum(C_ij*source_j), or fail."""
    source = [{(i,j):Q(c) for i,j,c in h} for h in source]
    target = [{(i,j):Q(c) for i,j,c in h} for h in target]
    if len(source) != len(target):
        raise ValueError('differential basis dimensions differ')
    monomials = sorted(set().union(*(set(h) for h in source+target)))
    n = len(source)
    result = []
    for h in target:
        rows = [[p.get(m,Q(0)) for p in source]+[h.get(m,Q(0))] for m in monomials]
        pivots = []
        at = 0
        for col in range(n):
            pivot = next((i for i in range(at,len(rows)) if rows[i][col]),None)
            if pivot is None:
                raise ValueError('dependent source differential basis')
            rows[at],rows[pivot] = rows[pivot],rows[at]
            scale = rows[at][col]
            rows[at] = [v/scale for v in rows[at]]
            for i in range(len(rows)):
                if i != at:
                    scale = rows[i][col]
                    rows[i] = [a-scale*b for a,b in zip(rows[i],rows[at])]
            pivots.append((at,col))
            at += 1
        if any(not any(row[:n]) and row[n] for row in rows):
            raise ValueError('differential polynomial spans do not match exactly')
        coefficients = [Q(0)]*n
        for row,col in pivots:
            coefficients[col] = rows[row][-1]
        result.append(coefficients)
    return result


def symplectic(matrix):
    """Use Python integers, not floating-point matrix products."""
    n = len(matrix)
    if n % 2 or any(len(row)!=n for row in matrix):
        return False
    g = n//2
    for i in range(n):
        for j in range(n):
            value = sum(matrix[k][i]*matrix[k+g][j]-matrix[k+g][i]*matrix[k][j]
                        for k in range(g))
            expected = int(i<g and j==i+g)-int(j<g and i==j+g)
            if value != expected:
                return False
    return True


def compare(ctx, sample, reference, digits):
    if sample['genus'] != reference['genus']:
        raise ValueError('genus differs from reference')
    change = basis_change(reference['basis'], sample['basis'])
    def decode(data):
        matrix = ctx.matrix([[ctx.mpc(a,b) for a,b in row] for row in data])
        if any(not ctx.isfinite(z) for z in matrix):
            raise ValueError('nonfinite period matrix')
        return matrix
    def realify(matrix):
        return ctx.matrix([[ctx.re(z) for z in row] for row in matrix.tolist()]+
                          [[ctx.im(z) for z in row] for row in matrix.tolist()])
    def maximum(matrix):
        return max(map(abs,matrix))
    with ctx.workdps(digits+40):
        a = decode(sample['periods'])
        c = ctx.matrix([[ctx.mpf(v.numerator)/v.denominator for v in row] for row in change])
        b = c*decode(reference['periods'])
        g = sample['genus']
        if a.rows!=g or a.cols!=2*g or b.rows!=g or b.cols!=2*g:
            raise ValueError('period matrix is not g by 2g')
        approx = realify(a)**-1*realify(b)
        integral = [[int(ctx.nint(v)) for v in row] for row in approx.tolist()]
        is_symplectic = symplectic(integral)
        error = maximum(a*ctx.matrix(integral)-b)
        scale = maximum(b)
        relative = error/scale
        tolerance = ctx.mpf(10)**(-digits+5)
        return dict(passed=bool(is_symplectic and relative<tolerance),
                    exact_symplectic=is_symplectic, relative_error=str(relative),
                    maximum_absolute_error=str(error), tolerance=str(tolerance),
                    integer_map_error=str(maximum(approx-ctx.matrix(integral))),
                    cycle_transform=integral,
                    differential_transform=[[str(v) for v in row] for row in change])
