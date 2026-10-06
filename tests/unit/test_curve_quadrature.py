from tests._support import make_curve
import pytest

from mpmath import mp
from generapy.curves.continuation import _continue_plane_curve_sheets_adaptive
from generapy.curves.integration import _integrate_plane_curve_path
from generapy.curves.polynomial import _prepare_plane_curve
from generapy.curves.quadrature import (
    _gauss_legendre_rule,
    _geometric_edge_panels, _geometric_quadrature_order, _legendre_edge_rule,
)


def test_gauss_rule_reuses_immutable_data_and_separates_precision_and_context(monkeypatch):
    ctx = mp.clone()
    ctx.dps = 20
    from generapy.curves import quadrature
    original = quadrature._legendre_edge_rule
    calls = []
    exposed = []

    def counted(context, order):
        if context is ctx:
            calls.append((ctx.prec, order))
        result = list(original(context, order))
        exposed.append(result)
        return result

    monkeypatch.setattr(quadrature, '_legendre_edge_rule', counted)
    first = _gauss_legendre_rule(ctx, 8)
    assert _gauss_legendre_rule(ctx, 8) is first
    # The builder output can change without mutating the cached rule.
    exposed[0][0] = (ctx.mpf(999), ctx.one)
    assert _gauss_legendre_rule(ctx, 8)[0][0] < 1
    with ctx.workdps(40):
        higher = _gauss_legendre_rule(ctx, 8)
        assert higher is not first
        assert abs(ctx.fsum(weight*node**15 for node, weight in zip(*higher))-ctx.mpf(1)/16) < 100*ctx.eps
    assert _gauss_legendre_rule(ctx, 8) is first
    assert len(calls) == 2
    other = mp.clone()
    other.dps = 20
    assert _gauss_legendre_rule(other, 8) is not first
    assert _gauss_legendre_rule(ctx, 12) is not first
    assert len(calls) == 3
    assert _gauss_legendre_rule.cache_info(ctx).maxsize == 32


@pytest.mark.parametrize('kind', ['path', 'iterated', 'chart'])
def test_repeated_integrations_reuse_rules_without_caching_form_values(monkeypatch, kind):
    from generapy.curves.continuation import _continue_plane_curve_branch
    from generapy.curves.integration import _integrate_plane_curve_branch, _integrate_plane_curve_path_iterated
    ctx = mp.clone()
    ctx.dps = 20
    curve = _prepare_plane_curve(ctx, {(0, 1): 1, (0, 0): -1})
    if kind == 'chart':
        lift = _continue_plane_curve_branch(ctx, curve, (0, 1), 1)
        integrate = _integrate_plane_curve_branch
    else:
        lift = _continue_plane_curve_sheets_adaptive(ctx, curve, (0, 1))
        integrate = _integrate_plane_curve_path_iterated if kind == 'iterated' else _integrate_plane_curve_path
    from generapy.curves import quadrature
    original = quadrature._legendre_edge_rule
    calls = []

    def counted(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(quadrature, '_legendre_edge_rule', counted)
    for coefficient in (1, 2):
        result = integrate(ctx, curve, lift, (lambda x, y: coefficient*x,), quadrature_order=8)
        assert abs(result.values[0]-ctx.mpf(coefficient)/2) < 100*ctx.eps
    assert len(calls) == 1


@pytest.mark.parametrize("dps", (20, 30, 40))
def test_geometry_quadrature_resolves_nearby_branch_value(dps):
    with mp.workdps(dps):
        branch = mp.mpc(1, "0.2")
        curve = _prepare_plane_curve(mp, {
            (0, 2): 1, (1, 0): -1, (0, 0): branch,
        })
        continuation = _continue_plane_curve_sheets_adaptive(mp, curve, (-1, 1))
        order = _geometric_quadrature_order(mp, -1, 1, (branch,))
        assert order > 12
        assert _geometric_quadrature_order(
            mp, -1, 0, (branch,)) < order
        exact = 2 * (continuation.fibres[-1][0]
                     - continuation.fibres[0][0])
        form = (lambda x, y: 1 / y,)
        coarse = _integrate_plane_curve_path(
            mp, curve, continuation, form, quadrature_order=3)
        local = _integrate_plane_curve_path(
            mp, curve, continuation, form, quadrature_order="geometry",
            branch_values=(branch,))
        assert abs(coarse.values[0] - exact) > mp.mpf("1e-8")
        assert abs(local.values[0] - exact) < 100 * mp.eps


def test_edge_panels_resolve_deep_local_refinement_accurately():
    from generapy.curves.quadrature import _geometric_edge_panels, _legendre_edge_rule
    with mp.workdps(30):
        branch = mp.mpc('0.137', '0.00001')
        panels = _geometric_edge_panels(mp, -1, 1, (branch,))
        assert panels[0][0] == 0 and panels[-1][1] == 1
        assert all(a[1] == b[0] for a, b in zip(panels, panels[1:]))
        assert min(b-a for a, b, n in panels) < mp.mpf(2)**-12
        assert len(panels) < 50
        rules = {n: _legendre_edge_rule(mp, n) for a, b, n in panels}
        actual = mp.fsum(2*(b-a)*weight/mp.sqrt(-1+2*(a+(b-a)*node)-branch)
                         for a, b, n in panels for node, weight in rules[n])
        exact = 2*(mp.sqrt(1-branch)-mp.sqrt(-1-branch))
        assert abs(actual-exact) < mp.mpf('1e-28')


def test_edge_panel_work_and_representability_limits(monkeypatch):
    import generapy.curves.quadrature as quadrature
    calls = []

    def high_order(*args):
        calls.append(1)
        return 128

    monkeypatch.setattr(quadrature, '_geometric_quadrature_order', high_order)
    with mp.workdps(20):
        with pytest.raises(mp.NoConvergence, match='panel budget'):
            quadrature._geometric_edge_panels(mp, 0, 1, (1j,), max_panels=3)
        assert len(calls) == 3
        with pytest.raises(mp.NoConvergence, match='midpoint is unresolved'):
            quadrature._geometric_edge_panels(mp, mp.one, mp.one+mp.eps, (1j,))
        with pytest.raises(ValueError, match='positive integer'):
            quadrature._geometric_edge_panels(mp, 0, 1, (1j,), max_panels=0)


@pytest.mark.parametrize("dps", (15, 40, 80))
@pytest.mark.parametrize("order", (3, 8, 24, 80, 160))
def test_edge_rule_matches_builtin_gaussian_quadrature(dps, order):
    ctx = mp.clone()
    ctx.dps = dps
    # The independent tridiagonal eigensolver supplies a reference for
    # both odd and even orders, including geometry-policy order 8.
    nodes, weights = ctx.gauss_quadrature(order, "legendre01")
    actual = _legendre_edge_rule(ctx, order)
    assert len(actual) == order
    for (node, weight), expected_node, expected_weight in zip(actual, nodes, weights):
        assert abs(node - expected_node) < 100*ctx.eps
        assert abs(weight - expected_weight) < 100*ctx.eps


@pytest.mark.parametrize("left,right,singularities,message", [
    (0, 0, (1j,), "distinct endpoints"),
    (0, 1, (), "known singularity"),
    (0, 1, (mp.mpf('.5'),), "meets a known singularity"),
])
def test_geometric_order_rejects_invalid_segments(left, right, singularities, message):
    with pytest.raises(ValueError, match=message):
        _geometric_quadrature_order(mp, left, right, singularities)


def test_geometric_panels_reject_degenerate_intervals_and_invalid_budget():
    with pytest.raises(ValueError, match="at least eight"):
        _geometric_edge_panels(mp, 0, 1, (1j,), max_order=4)
    with pytest.raises(mp.NoConvergence, match="endpoints are unresolved"):
        _geometric_edge_panels(mp, mp.one, mp.one, (1j,))
    with pytest.raises(ValueError, match="positive integer"):
        _geometric_edge_panels(mp, 0, 1, (1j,), max_panels=mp.mpf('1.5'))


@pytest.mark.parametrize("order", (0, 1, mp.mpf('2.5')))
def test_edge_rule_rejects_invalid_order(order):
    with pytest.raises(ValueError, match="integer at least two"):
        _legendre_edge_rule(mp, order)


def test_edge_rule_bounds_newton_work_after_a_bad_initial_estimate(monkeypatch):
    ctx = mp.clone()
    monkeypatch.setattr(ctx, "cos", lambda value: ctx.mpf("1e100"))
    with pytest.raises(ctx.NoConvergence, match="Legendre node did not converge"):
        _legendre_edge_rule(ctx, 8)
