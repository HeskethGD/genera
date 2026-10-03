"""Public curve contracts that do not depend on theta or Kleinian functions."""

from tests._support import make_curve, with_basis

import warnings

import pytest

import genera.curves._operations as curve_operations
from genera.curves import (
    CurveBranchLocus, CurveGenus, CurveHomology, CurveMonodromy, CurvePlace,
    CurveRiemannConstant, CurveAbelMapKind1, CurveAbelMapKind2,
)
from genera.curves.polynomial import _prepare_plane_curve
from genera.curves._stages import _stage_geometric_periods
from genera.curves.jacobian import _finite_geometric_abel_value
from mpmath import mp


@pytest.mark.parametrize("polynomial,message", [
    (None, "sparse mapping"),
    ((), "sparse mapping"),
    ((0, -1, 0, 1), "sparse mapping"),
    ([(0, 2, 1), (3, 0, -1)], "sparse mapping"),
    ({}, "must be nonzero"),
    ({(0, 2): 0}, "must be nonzero"),
    ({(2, 0): 1}, "must depend on y"),
    ({(0, -1): 1}, "nonnegative integers"),
    ({"y": 1}, "nonnegative integers"),
    ({(0, 2): "invalid"}, "coefficients must be numbers"),
    (((0, 2, "invalid"),), "sparse mapping"),
    ({(0, 2): "inf"}, "coefficients must be finite"),
    ({(0, 2): "nan"}, "coefficients must be finite"),
])
def test_curve_constructor_rejects_invalid_polynomial(polynomial, message):
    with pytest.raises(ValueError, match=message):
        make_curve(mp, polynomial)


def test_curve_materializes_input_and_returns_a_copy_of_the_polynomial():
    ctx = mp.clone()
    source = {(0, 2): 1, (3, 0): -1, (1, 0): 1}
    expected = dict(source)
    curve = make_curve(ctx, source)
    source[(3, 0)] = -2
    exported = curve.polynomial
    exported[(3, 0)] = -3
    assert curve.polynomial == expected
    assert repr(curve) == (
        f"Curve(x_degree=3, y_degree=2, ctx.prec={ctx.prec})")

    with ctx.workdps(25), pytest.warns(UserWarning, match="context changed"):
        assert curve.branch_locus.branch_values == (-1, 0, 1)


@pytest.mark.parametrize("attribute,value", [
    ("rounding", "d"), ("trap_complex", True),
])
def test_curve_context_state_changes_warn_once_and_restore(attribute, value):
    ctx = mp.clone()
    curve = make_curve(ctx, {(0, 2): 1, (1, 0): -1})
    original = curve.branch_locus
    previous = getattr(ctx, attribute)
    setattr(ctx, attribute, value)
    with pytest.warns(UserWarning, match="context changed"):
        assert curve.branch_locus == original
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        assert curve.branch_locus == original
        setattr(ctx, attribute, previous)
        assert curve.branch_locus == original
        setattr(ctx, attribute, value)
        assert curve.branch_locus == original
    assert caught == []


def test_curve_recomputes_exact_inputs_without_inventing_coefficient_precision():
    ctx = mp.clone()
    ctx.dps = 15
    rounded = ctx.mpf("0.1")
    exact = make_curve(ctx, {(0, 1): 10, (0, 0): -1})
    inexact = make_curve(ctx, {(0, 1): 1, (0, 0): -rounded})
    with ctx.workdps(40):
        with pytest.warns(UserWarning, match="context changed"):
            exact_y = exact.fibre(0)[0].y
        with pytest.warns(UserWarning, match="context changed"):
            inexact_y = inexact.fibre(0)[0].y
        assert abs(exact_y - ctx.mpf("0.1")) < ctx.eps
        assert inexact_y == rounded
        assert abs(inexact_y - exact_y) > ctx.mpf("1e-20")


def test_period_results_preserve_basis_and_do_not_share_mutable_matrices():
    ctx = mp.clone()
    ctx.dps = 18
    curve = make_curve(ctx, {(0, 2): 1, (1, 0): 1, (3, 0): -1})
    first = curve.periods_kind_1()
    expected = tuple(+getattr(first, name) for name in ("omega", "omega_prime", "tau"))
    for name in ("omega", "omega_prime", "tau"):
        getattr(first, name)[0, 0] = 0
    # A supplied basis selects a different marking without changing the
    # default basis or contaminating its cached matrices.
    forms = (lambda x, y: 2 / y,)
    supplied = with_basis(curve, differentials_kind_1=forms).periods_kind_1()
    assert supplied.differentials == forms
    assert supplied.marking == "geometric-polygon"
    again = curve.periods_kind_1()
    assert again.marking == "baker"
    for name, value in zip(("omega", "omega_prime", "tau"), expected):
        assert getattr(again, name) == value
    second = curve.periods_kind_2()
    saved = tuple(+getattr(second, name) for name in ("eta", "eta_prime", "kappa"))
    for name in ("eta", "eta_prime", "kappa"):
        getattr(second, name)[0, 0] = 0
    second_again = curve.periods_kind_2()
    for name, value in zip(("eta", "eta_prime", "kappa"), saved):
        assert getattr(second_again, name) == value
    assert curve.periods_kind_1().marking == "baker"
    constant = curve.riemann_constant()
    constant.value[0] = 0
    assert ctx.almosteq(curve.riemann_constant().value[0], (1 + expected[2][0, 0]) / 2)


def test_curve_validation_reports_failed_topological_checks():
    curve = make_curve(mp, {(0, 2): 1, (1, 0): -1})
    cases = [
        (CurveBranchLocus(2, (0, 0), (0, 1)), "branch_values_distinct"),
        (CurveBranchLocus(2, (0,), (1,)), "resultant_nonconstant"),
        (CurveGenus(1, 2, 2), "riemann_hurwitz_balance"),
        (CurveMonodromy(1, (-1, 1), (0,), ((0, 1),), (0, 1),
                       2, 0, False, True, 1), "monodromy_transitive"),
        (CurveMonodromy(1, (-1, 1), (0,), ((1, 0),), (0, 1),
                       2, 0, True, False, 1), "monodromy_product_identity"),
        (CurveHomology(1, 2, 0, 2, 0, ((0, 1), (1, 0)),
                       ((1, 0), (0, 1)), "general", "geometric-polygon"),
         "intersection_form_antisymmetric"),
        (CurveHomology(1, 2, 0, 2, 0, ((0, 0), (0, 0)),
                       ((1, 0), (0, 1)), "general", "geometric-polygon"),
         "intersection_form_rank"),
    ]
    for result, name in cases:
        report = curve.validate(result)
        assert not report.passed
        assert report.maximum_residual is None
        assert not next(check.passed for check in report.checks if check.name == name)
    with pytest.raises(TypeError, match="curve result record"):
        curve.validate(mp.eye(1))


@pytest.mark.parametrize("characteristic", [None, ((), ()), ((mp.inf,), (0,))])
def test_riemann_constant_validation_rejects_invalid_characteristics(characteristic):
    curve = make_curve(mp, {(0, 2): 1, (1, 0): 1, (3, 0): -1})
    result = CurveRiemannConstant(
        mp.matrix([mp.mpc('.5', '.5')]), characteristic, None, None,
        "hyperelliptic", "baker")
    report = curve.validate(result)
    assert not report.passed
    assert not next(check.passed for check in report.checks
                    if check.name == "characteristic_finite")


@pytest.mark.parametrize("value", [mp.matrix([[0, 1]]), mp.matrix([mp.inf])])
def test_riemann_constant_validation_rejects_invalid_values(value):
    curve = make_curve(mp, {(0, 2): 1, (1, 0): 1, (3, 0): -1})
    result = CurveRiemannConstant(
        value, ((0,), (0,)), None, None, "hyperelliptic", "baker")
    assert not curve.validate(result).passed


def test_curve_validation_uses_recorded_residual_and_current_precision():
    ctx = mp.clone()
    ctx.dps = 15
    curve = make_curve(ctx, {(0, 2): 1, (1, 0): 1, (3, 0): -1})
    residual = ctx.mpf("1e-8")
    result = CurveRiemannConstant(
        ctx.matrix([ctx.mpc('.5', '.5')]), ((ctx.mpf('.5'),), (ctx.mpf('.5'),)),
        None, residual, "general", "geometric-polygon")
    assert curve.validate(result).passed
    with ctx.workdps(40), pytest.warns(UserWarning, match="context changed"):
        report = curve.validate(result)
    assert not report.passed
    assert report.maximum_residual == residual
    assert not next(check.passed for check in report.checks
                    if check.name == "max_sheet_residual")


@pytest.mark.parametrize("periods", [None, [[1, 2]], [["invalid"]]])
def test_lattice_reduction_rejects_invalid_periods(periods):
    curve = make_curve(mp, {(0, 2): 1, (1, 0): 1, (3, 0): -1})
    with pytest.raises(ValueError, match="periods must be square"):
        curve.lattice_reduce([0], periods)


@pytest.mark.parametrize("value", [None, [0, 1], [[0, 1]], ["invalid"]])
def test_lattice_reduction_rejects_invalid_vectors(value):
    curve = make_curve(mp, {(0, 2): 1, (1, 0): 1, (3, 0): -1})
    with pytest.raises(ValueError, match="genus-length column vector"):
        curve.lattice_reduce(value, mp.matrix([[1j]]))


@pytest.mark.parametrize("forms", [0, (), (1,)])
def test_public_integrals_reject_invalid_differential_sequences(forms):
    curve = make_curve(mp, {(0, 2): 1, (1, 0): -1})
    path = curve.path((1, 1), (4, 2))
    with pytest.raises(ValueError, match="sequence of callables"):
        curve.integral(forms, path)


def test_public_integral_requires_a_lifted_path():
    curve = make_curve(mp, {(0, 2): 1, (1, 0): -1})
    with pytest.raises(TypeError, match="must be a CurvePath"):
        curve.integral(lambda x, y: 1, [(1, 1), (4, 2)])


@pytest.mark.parametrize("place,message", [
    (None, "finite .* pair"),
    ((mp.inf, 1), "finite .* pair"),
    ((0, 0), "regular for the x projection"),
])
def test_public_path_rejects_nonregular_endpoints(place, message):
    curve = make_curve(mp, {(0, 2): 1, (1, 0): -1})
    with pytest.raises(ValueError, match=message):
        curve.path(place, (1, 1))


def test_second_kind_only_basis_selects_consistent_geometric_marking():
    ctx = mp.clone()
    ctx.dps = 18
    curve = make_curve(ctx, {(0, 2): 1, (1, 0): 1, (3, 0): -1},
                       differentials_kind_2=(lambda x, y: x / y,))
    first = curve.periods_kind_1()
    second = curve.periods_kind_2()
    assert curve.differentials_kind_1 is None
    assert first.marking == second.marking == curve.homology.marking == "geometric-polygon"
    assert curve.riemann_constant().marking == first.marking
    target = (ctx.mpf(2), ctx.sqrt(6))
    assert curve.abel_map_kind_1(target).marking == first.marking
    assert curve.abel_map_kind_2(target).marking == first.marking


def test_hyperelliptic_chart_endpoints_require_supplied_forms():
    curve = make_curve(mp, {(0, 2): 1, (1, 0): 1, (3, 0): -1})
    chart = curve.monomial_chart(-2, -3)
    place = curve.chart_place(chart, 1, mp.mpf('.1'))
    for target in (place, [place]):
        with pytest.raises(ValueError, match="chart-backed places require the general pipeline"):
            curve.abel_map_kind_1(target).value
    with pytest.raises(ValueError, match="base_place must be one affine point"):
        curve.abel_map_kind_1([], base_place=[]).value


def test_hyperelliptic_second_kind_unreduced_result_and_scalar_target():
    ctx = mp.clone()
    ctx.dps = 18
    curve = make_curve(ctx, {(0, 2): 1, (1, 0): 1, (3, 0): -1})
    target = (ctx.mpf(2), ctx.sqrt(6))
    result = curve.abel_map_kind_2(target)
    assert result.reduction_shift is None
    assert result.engine == "hyperelliptic"
    assert all(ctx.isfinite(entry) for entry in result.value)
    with pytest.raises(ValueError, match="affine point"):
        curve.abel_map_kind_1(2).value


def test_geometric_periods_reject_nonpositive_imaginary_part(monkeypatch):
    ctx = mp.clone()
    ctx.dps = 18
    curve = _prepare_plane_curve(ctx, {(0, 2): 1, (3, 0): -1, (1, 0): 1})
    monkeypatch.setattr(
        curve_operations, "_tau_imaginary_eigenvalues",
        lambda context, tau: (-context.one,))
    with pytest.raises(ValueError, match="not positive definite"):
        curve_operations._geometric_first_kind_periods(ctx, curve)


def test_geometric_abel_open_path_rejects_an_unmatched_place():
    ctx = mp.clone()
    ctx.dps = 18
    curve = _prepare_plane_curve(ctx, {(0, 2): 1, (3, 0): -1, (1, 0): 1})
    data = _stage_geometric_periods(ctx, curve)
    with pytest.raises(ValueError, match="could not be matched"):
        _finite_geometric_abel_value(ctx, curve, data, CurvePlace(2, 0), ())


@pytest.mark.parametrize("geometric", [False, True])
@pytest.mark.parametrize("explicit_base", [False, True])
def test_first_kind_abel_record_retains_reduction_shift(geometric, explicit_base):
    ctx = mp.clone()
    ctx.dps = 18
    curve = make_curve(ctx, {(0, 2): 1, (1, 0): 1, (3, 0): -1})
    forms = (lambda x, y: 1 / y,) if geometric else None
    point = (ctx.mpf(2), ctx.sqrt(6))
    base = (ctx.mpf(3), ctx.sqrt(24)) if explicit_base else None
    raw = with_basis(curve, differentials_kind_1=forms).abel_map_kind_1([point] * 5, base_place=base)
    reduced = with_basis(curve, differentials_kind_1=forms).abel_map_kind_1([point] * 5, base_place=base, reduce=True)
    assert isinstance(raw, CurveAbelMapKind1)
    assert raw._fields == CurveAbelMapKind2._fields
    assert raw.reduction_shift is None
    assert reduced.engine == raw.engine == ("general" if geometric else "hyperelliptic")
    assert reduced.marking == raw.marking == ("geometric-polygon" if geometric else "baker")
    periods = with_basis(curve, differentials_kind_1=forms).periods_kind_1()
    m, n = reduced.reduction_shift
    assert isinstance(m, int) and isinstance(n, int)
    assert ctx.norm(raw.value - reduced.value - 2 * periods.omega * m
                    - 2 * periods.omega_prime * n) < ctx.mpf("1e-14")


def test_unreduced_abel_record_does_not_request_periods(monkeypatch):
    ctx = mp.clone()
    ctx.dps = 18
    curve = make_curve(ctx, {(0, 2): 1, (1, 0): 1, (3, 0): -1})

    def periods_forbidden(*args, **kwargs):
        raise AssertionError("an unreduced Abel map must not request period matrices")

    monkeypatch.setattr(curve_operations, "periods", periods_forbidden)
    from genera.curves._hyperelliptic import operations as specialized
    monkeypatch.setattr(specialized, "_first_kind_periods", periods_forbidden)
    for base in (None, (ctx.mpf(3), ctx.sqrt(24))):
        result = curve.abel_map_kind_1((ctx.mpf(2), ctx.sqrt(6)), base_place=base)
        assert result.reduction_shift is None
        assert all(ctx.isfinite(entry) for entry in result.value)


@pytest.mark.parametrize("argument", ["differentials_kind_1", "differentials_kind_2"])
@pytest.mark.parametrize("forms", [0, (), (1,)])
def test_constructor_rejects_invalid_bound_bases(argument, forms):
    with pytest.raises(ValueError, match=argument + " must be a sequence of callables"):
        make_curve(mp, {(0, 2): 1, (1, 0): -1}, **{argument: forms})


def test_constructor_materializes_bases_without_running_numerical_stages(monkeypatch):
    def numerical_stage_forbidden(*args, **kwargs):
        raise AssertionError("constructing a bound basis must remain lazy")

    monkeypatch.setattr(curve_operations, "periods", numerical_stage_forbidden)
    monkeypatch.setattr(curve_operations, "homology", numerical_stage_forbidden)
    du = lambda x, y: 1 / y
    dr = lambda x, y: x / y
    source = [du]
    second_source = (form for form in (dr,))
    curve = make_curve(mp, {(0, 2): 1, (1, 0): 1, (3, 0): -1},
                       differentials_kind_1=source, differentials_kind_2=second_source)
    source.clear()
    assert curve.differentials_kind_1 == (du,)
    assert curve.differentials_kind_2 == (dr,)
    assert tuple(second_source) == ()
    with pytest.raises(AttributeError):
        curve.differentials_kind_1 = (dr,)
    with pytest.raises(AttributeError):
        curve.differentials_kind_2 = (du,)


@pytest.mark.parametrize("polynomial", [
    {(0, 2): 1, (1, 0): 1, (3, 0): -1},
    {(0, 3): 1, (4, 0): -1, (0, 0): 1},
])
def test_missing_second_kind_basis_explains_constructor_remedy(polynomial, monkeypatch):
    curve = make_curve(mp, polynomial, differentials_kind_1=(lambda x, y: 1 / y,))

    def integration_forbidden(*args, **kwargs):
        raise AssertionError("missing second-kind basis should fail before integration")

    monkeypatch.setattr(curve_operations, "periods", integration_forbidden)
    monkeypatch.setattr(curve_operations, "abel_map", integration_forbidden)
    for method in (curve.periods_kind_2, lambda: curve.abel_map_kind_2([])):
        with pytest.raises(ValueError, match="Construct Curve.*differentials_kind_2"):
            method()


@pytest.mark.parametrize("method,args,kwargs", [
    ("periods_kind_1", (), {"differentials_kind_1": ()}),
    ("periods_kind_2", (), {"differentials_kind_2": ()}),
    ("riemann_matrix", ((),), {}),
    ("riemann_constant", (), {"differentials_kind_1": ()}),
    ("abel_map_kind_1", ([], ()), {}),
    ("abel_map_kind_2", ([],), {"differentials_kind_2": ()}),
])
def test_bound_basis_methods_reject_per_call_overrides(method, args, kwargs):
    curve = make_curve(mp, {(0, 2): 1, (1, 0): 1, (3, 0): -1})
    with pytest.raises(TypeError):
        getattr(curve, method)(*args, **kwargs)


def test_bound_custom_basis_reuses_its_period_cache_without_mutable_aliases():
    ctx = mp.clone()
    ctx.dps = 18
    forms = (lambda x, y: 1 / y,)
    curve = make_curve(ctx, {(0, 2): 1, (1, 0): 1, (3, 0): -1},
                       differentials_kind_1=forms, differentials_kind_2=(lambda x, y: x / y,))
    curve.periods_kind_2()
    first = curve.periods_kind_1()
    expected = +first.omega
    first.omega[0, 0] = 0
    assert curve.periods_kind_1().omega == expected
    assert curve.homology.marking == curve.periods_kind_1().marking
