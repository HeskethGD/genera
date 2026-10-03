"""Focused checks for the specialized and general curve dispatch."""

import pytest
from mpmath import mp

from genera import rtheta
from genera.curves import _operations
from genera.curves.algebraic_curve import CurvePlace
from tests._support import make_curve, with_basis

TERMS = {(0, 3): 1, (4, 0): -1, (1, 0): 1, (0, 0): -1}


def test_geometric_interface_keeps_one_marking(monkeypatch):
    ctx = mp.clone()
    ctx.dps = 18
    curve = make_curve(ctx, TERMS)

    def radial_forbidden(*args, **kwargs):
        raise AssertionError('unexpected radial fallback')

    for name in ('_stage_monodromy',):
        monkeypatch.setattr(_operations, name, radial_forbidden)
    periods = curve.periods_kind_1()
    assert curve.genus == 3
    for result in (periods, curve.homology, curve.riemann_constant()):
        assert result.marking == 'geometric-polygon'
        assert curve.validate(result).passed
    assert curve.validate(curve.genus_data).passed
    assert ctx.norm(periods.tau-curve.riemann_matrix()) == 0
    original = periods.omega[0, 0]
    periods.omega[0, 0] = 123
    periods = curve.periods_kind_1()
    assert periods.omega[0, 0] == original
    target = curve.fibre(ctx.mpc('.3', '.7'))[0]
    base = curve.fibre(ctx.mpc('.4', '.8'))[-1]
    value = curve.abel_map_kind_1([target, target], base_place=base).value
    assert ctx.norm(value - 2*curve.abel_map_kind_1(target, base_place=base).value) < ctx.mpf('1e-16')
    assert ctx.norm(curve.abel_map_kind_1([]).value) == 0
    assert ctx.norm(curve.abel_map_kind_1(base, base_place=base).value) == 0
    reduced = curve.abel_map_kind_1([target, target], base_place=base, reduce=True).value
    assert ctx.norm(reduced-curve.lattice_reduce(value, periods).value) < ctx.mpf('1e-16')
    constant = curve.riemann_constant(base_place=base)
    argument = (2*periods.omega)**-1 * value + constant.value
    argument = curve.lattice_reduce(argument, periods.tau).value
    assert abs(rtheta(argument, periods.tau, ctx=ctx)) < ctx.mpf('1e-12')
    with ctx.workdps(23):
        with pytest.warns(UserWarning, match='context changed'):
            higher = curve.periods_kind_1()
        assert higher.marking == 'geometric-polygon'
        assert curve.homology.marking == higher.marking
        assert ctx.norm(higher.omega-periods.omega) < ctx.mpf('1e-16')
    assert curve.periods_kind_1().omega[0, 0] == original


def test_geometric_unsupported_operations_are_explicit():
    ctx = mp.clone()
    curve = make_curve(ctx, TERMS)
    forms = (lambda x, y: 1,)
    for operation in (
        lambda: with_basis(curve, differentials_kind_1=forms).periods_kind_1(),
        lambda: with_basis(curve, differentials_kind_1=forms).riemann_matrix(),
        lambda: with_basis(curve, differentials_kind_1=forms).riemann_constant(),
        lambda: with_basis(curve, differentials_kind_1=forms).abel_map_kind_1([]).value,
    ):
        with pytest.raises(ValueError, match="one form per positive genus"):
            operation()
    for operation in (
        lambda: curve.periods_kind_2(),
        lambda: curve.abel_map_kind_2([]),
    ):
        with pytest.raises(ValueError, match='Construct Curve.*differentials_kind_2'):
            operation()
    chart_place = CurvePlace(0, 1, object())
    for operation in (
        lambda: curve.abel_map_kind_1(chart_place).value,
        lambda: curve.abel_map_kind_1([], base_place=chart_place).value,
        lambda: curve.riemann_constant(base_place=chart_place),
    ):
        with pytest.raises(ValueError, match='invalid chart description'):
            operation()
    with pytest.raises(TypeError, match='_general_backend'):
        make_curve(ctx, TERMS, _general_backend='unknown')


def test_hyperelliptic_automatic_dispatch_stays_specialized(monkeypatch):
    ctx = mp.clone()
    ctx.dps = 18

    def geometric_forbidden(*args, **kwargs):
        raise AssertionError('hyperelliptic curve entered geometric backend')

    for name in ('_stage_geometric_periods', '_stage_geometric_polygon',
                 '_stage_geometric_riemann_constant'):
        monkeypatch.setattr(_operations, name, geometric_forbidden)
    curve = make_curve(ctx, {(0, 2): 1, (3, 0): -1, (1, 0): 1})
    for result in (curve.periods_kind_1(), curve.periods_kind_2(),
                   curve.homology, curve.riemann_constant()):
        assert result.engine == 'hyperelliptic'
        assert result.marking == 'baker'
    default = make_curve(ctx, {(0, 2): 1, (1, 0): 1, (3, 0): -1})
    point = (ctx.mpf(2), ctx.sqrt(6))
    assert ctx.norm(curve.abel_map_kind_1(point).value-default.abel_map_kind_1(point).value) == 0


def test_constructors_agree_and_monodromy_remains_diagnostic():
    from genera.curves.algebraic_curve import Curve
    ctx = mp.clone()
    ctx.dps = 18
    terms = {(0, 3): 1, (3, 0): 1, (0, 0): -1}
    default = make_curve(ctx, terms)
    direct = Curve(terms, ctx=ctx)
    assert default.periods_kind_1().marking == 'geometric-polygon'
    assert ctx.norm(default.riemann_matrix()-direct.riemann_matrix()) == 0
    assert default.monodromy == direct.monodromy
    assert default.validate(default.monodromy).passed
    constructors = (
        lambda t, **kw: make_curve(ctx, t, **kw),
        lambda t, **kw: Curve(t, **kw, ctx=ctx),
    )
    for constructor in constructors:
        with pytest.raises(TypeError, match='_general_backend'):
            constructor(terms, _general_backend='radial')


def test_geometric_chart_endpoints_cutoffs_and_theta():
    ctx = mp.clone()
    ctx.dps = 18
    curve = make_curve(ctx, {(0, 3): 1, (4, 0): -1, (0, 0): 1})
    infinity = curve.monomial_chart(-3, -4)
    branch = curve.chart(
        {(0, 3): 1, (0, 0): -4, (3, 0): -6, (6, 0): -4, (9, 0): -1},
        lambda t, w: (1+t**3, t*w, 3*t**2))
    periods = curve.periods_kind_1()
    inverse = (2*periods.omega)**-1
    tau = periods.tau
    places = []
    for chart, seed in ((infinity, ctx.one), (branch, ctx.root(4, 3))):
        outer = curve.chart_place(chart, seed, ctx.mpf('.3'))
        inner = curve.chart_place(chart, seed, ctx.mpf('.25'))
        first = curve.abel_map_kind_1(outer).value
        second = curve.abel_map_kind_1(inner).value
        # Independently chosen open paths can differ by a period.
        difference = curve.lattice_reduce(first-second, periods).value
        assert ctx.norm(difference) < ctx.mpf('1e-14')
        local = curve.chart_integral(
            chart, periods.differentials, (0, ctx.mpf('.3')), seed)
        junction = curve.abel_map_kind_1((outer.x, outer.y)).value
        assert ctx.norm(first-junction+ctx.matrix(local.values)) < ctx.mpf('1e-14')
        places.append(outer)
    constant = curve.riemann_constant()
    # g=3: the divisor consists of one ramification point and infinity.
    argument = inverse*curve.abel_map_kind_1(places).value + constant.value
    reduced = curve.lattice_reduce(argument, tau).value
    assert abs(rtheta(reduced, tau, ctx=ctx)) < ctx.mpf('1e-11')
    shifted = curve.riemann_constant(base_place=places[0])
    argument = inverse*curve.abel_map_kind_1(places, base_place=places[0]).value+shifted.value
    reduced = curve.lattice_reduce(argument, tau).value
    assert abs(rtheta(reduced, tau, ctx=ctx)) < ctx.mpf('1e-11')
    assert ctx.norm(curve.abel_map_kind_1(places[0], base_place=places[0]).value) == 0
    other = make_curve(ctx, TERMS)
    with pytest.raises(ValueError, match='different curve or precision'):
        other.abel_map_kind_1(places[0]).value
    with ctx.workdps(23):
        with pytest.warns(UserWarning, match='context changed'):
            with pytest.raises(ValueError, match='different curve or precision'):
                curve.abel_map_kind_1(places[0]).value
    assert ctx.dps == 18


def test_geometric_supplied_basis_is_coherent_across_operations(monkeypatch):
    from genera.curves import _stages
    ctx = mp.clone()
    ctx.dps = 18
    curve = make_curve(ctx, TERMS)
    automatic = curve.periods_kind_1()
    constant = curve.riemann_constant()
    point = curve.fibre(ctx.mpc('.3', '.7'))[0]
    automatic_value = curve.abel_map_kind_1(point).value
    chart = curve.monomial_chart(-3, -4)
    infinity = curve.chart_place(chart, 1, ctx.mpf('.3'))
    automatic_infinity = curve.abel_map_kind_1(infinity).value
    change = ctx.matrix([[2, 1, 0], [0, 1, 1], [1, 0, 1]])

    class Form:
        def __init__(self, row):
            self.row = row

        def __call__(self, x, y):
            return ctx.fsum(change[self.row, j]*f(x, y)
                            for j, f in enumerate(automatic.differentials))

    forms = tuple(Form(i) for i in range(3))

    def automatic_forbidden(*args, **kwargs):
        raise AssertionError('custom forms require no automatic Baker basis')

    monkeypatch.setattr(_stages, '_stage_geometric_periods_working', automatic_forbidden)
    monkeypatch.setattr(_operations, '_stage_monodromy', automatic_forbidden)
    supplied = with_basis(curve, differentials_kind_1=forms).periods_kind_1()
    assert supplied.differentials == forms
    assert supplied.marking == 'geometric-polygon'
    assert ctx.norm(supplied.omega-change*automatic.omega) < ctx.mpf('1e-15')
    assert ctx.norm(supplied.omega_prime-change*automatic.omega_prime) < ctx.mpf('1e-15')
    assert ctx.norm(supplied.tau-automatic.tau) < ctx.mpf('1e-15')
    assert ctx.norm(with_basis(curve, differentials_kind_1=forms).riemann_constant().value-constant.value) < ctx.mpf('1e-14')
    value = with_basis(curve, differentials_kind_1=forms).abel_map_kind_1(point).value
    assert ctx.norm(value-change*automatic_value) < ctx.mpf('1e-15')
    chart_value = with_basis(curve, differentials_kind_1=forms).abel_map_kind_1(infinity).value
    assert ctx.norm(chart_value-change*automatic_infinity) < ctx.mpf('1e-14')
    based = with_basis(curve, differentials_kind_1=forms).abel_map_kind_1(point, base_place=infinity).value
    assert ctx.norm(based-(value-chart_value)) < ctx.mpf('1e-14')
    shifted = with_basis(curve, differentials_kind_1=forms).riemann_constant(base_place=infinity)
    expected = constant.value+2*((2*supplied.omega)**-1*chart_value)
    assert ctx.norm(shifted.value-expected) < ctx.mpf('1e-14')
    reduced = with_basis(curve, differentials_kind_1=forms).abel_map_kind_1(point, reduce=True).value
    assert ctx.norm(reduced-curve.lattice_reduce(value, supplied).value) < ctx.mpf('1e-14')


def test_geometric_supplied_basis_accepts_unhashable_forms():
    ctx = mp.clone()
    ctx.dps = 18
    curve = make_curve(ctx, {(0, 2): 1, (3, 0): -1, (1, 0): 1})

    class Form:
        __hash__ = None

        def __call__(self, x, y):
            return 1/y

    form = Form()
    periods = with_basis(curve, differentials_kind_1=(form,)).periods_kind_1()
    assert periods.differentials == (form,)
    assert periods.marking == 'geometric-polygon'
    assert curve.validate(periods).passed


@pytest.mark.parametrize('linear_y', [False, True])
def test_hyperelliptic_supplied_forms_use_geometric_marking(linear_y):
    ctx = mp.clone()
    ctx.dps = 20
    terms = {(0,2):1, (3,0):-1, (1,0):1}
    if linear_y:
        terms.update({(1,1):2, (2,0):1})
    curve = make_curve(ctx, terms)
    forms = (lambda x,y: 1/(y+x if linear_y else y),)
    baker = curve.periods_kind_1()
    general = with_basis(curve, differentials_kind_1=forms).periods_kind_1()
    assert baker.marking == curve.homology.marking == 'baker'
    assert general.marking == 'geometric-polygon'
    assert curve.validate(general).passed
    def full(p):
        return ctx.matrix([[2*p.omega[0,0], 2*p.omega_prime[0,0]]])
    def realify(p):
        return ctx.matrix([[ctx.re(z) for z in p.tolist()[0]],
                           [ctx.im(z) for z in p.tolist()[0]]])
    a,b = full(general),full(baker)
    change = (realify(a)**-1*realify(b)).apply(ctx.nint)
    assert ctx.det(change) == 1
    assert ctx.norm(a*change-b) < ctx.mpf('1e-17')
    base = (ctx.mpf(2),ctx.sqrt(6)-(2 if linear_y else 0))
    target = (ctx.mpf(3),ctx.sqrt(24)-(3 if linear_y else 0))
    value = with_basis(curve, differentials_kind_1=forms).abel_map_kind_1(target, base_place=base).value
    expected = ctx.quad(lambda x:1/ctx.sqrt(x**3-x),[2,3])
    assert abs(value[0]-expected) < ctx.mpf('1e-17')
    constant = with_basis(curve, differentials_kind_1=forms).riemann_constant(base_place=base)
    assert constant.marking == general.marking
    assert abs(rtheta(
        constant.value, general.tau, ctx=ctx)) < ctx.mpf('1e-16')
    second = (lambda x,y: 1,)
    periods = with_basis(curve, differentials_kind_1=forms, differentials_kind_2=second).periods_kind_2()
    assert periods.marking == general.marking
    assert ctx.norm(periods.eta)+ctx.norm(periods.eta_prime) < ctx.mpf('1e-17')
    integral = with_basis(curve, differentials_kind_1=forms, differentials_kind_2=second).abel_map_kind_2(target, base_place=base, reduce=True)
    assert integral.marking == general.marking
    assert abs(integral.value[0]-1) < ctx.mpf('1e-17')
    assert curve.periods_kind_1().marking == 'baker'
