"""Curve-based Kleinian inputs preserve explicit-period conventions."""

import pytest
from mpmath import mp

from generapy import (
    Curve, kleinian_p, kleinian_sigma, kleinian_sigma_jet,
    kleinian_sigma_normalization, kleinian_zeta,
)
from generapy import kleinian as implementation


POLYNOMIAL = {(0, 2): 1, (3, 0): -4, (1, 0): 4}


@pytest.fixture
def curve():
    ctx = mp.clone()
    ctx.dps = 22
    return Curve(POLYNOMIAL, ctx=ctx)


def test_curve_and_mapping_agree_with_explicit_periods(curve):
    ctx = curve.ctx
    first = curve.periods_kind_1()
    second = curve.periods_kind_2()
    characteristic = curve.riemann_constant().characteristic
    args = ([ctx.mpf('0.23')], first.omega, first.tau, second.kappa)
    expected = kleinian_sigma(*args, characteristic,
                              normalization="hyperelliptic", ctx=ctx)
    assert ctx.almosteq(kleinian_sigma(args[0], curve=curve), expected)
    assert ctx.almosteq(kleinian_sigma(args[0], curve=POLYNOMIAL, ctx=ctx), expected)
    expected_jet = kleinian_sigma_jet(*args, 2, characteristic,
                                    normalization="hyperelliptic", ctx=ctx)
    jet = kleinian_sigma_jet(args[0], curve=curve, order=2)
    assert all(ctx.almosteq(jet[index], value)
               for index, value in expected_jet.items())
    assert ctx.almosteq(kleinian_zeta(args[0], curve=curve)[0],
                       kleinian_zeta(*args, characteristic, ctx=ctx)[0])
    assert ctx.almosteq(kleinian_p(args[0], curve=curve, indices=(0, 0)),
                       kleinian_p(*args, (0, 0), characteristic, ctx=ctx))
    constant = kleinian_sigma_normalization(curve=curve)
    assert ctx.almosteq(constant, kleinian_sigma_normalization(
        first.omega, first.tau, characteristic, "hyperelliptic", ctx=ctx))
    raw = kleinian_sigma_jet([0], curve=curve, order=1, normalization="theta")
    assert ctx.almosteq(constant * raw[(1,)], 1)
    assert kleinian_sigma_normalization(first.omega, first.tau, ctx=ctx) == 1


def test_logarithmic_functions_skip_scalar_normalization(curve, monkeypatch):
    def unnecessary(*args, **kwargs):
        raise AssertionError("scalar normalization must not be computed")
    monkeypatch.setattr(implementation, 'kleinian_sigma_normalization', unnecessary)
    monkeypatch.setattr(implementation, '_hyperelliptic_sigma_normalization', unnecessary)
    assert curve.ctx.isfinite(kleinian_zeta([.2], curve=curve)[0])
    assert curve.ctx.isfinite(kleinian_p([.2], curve=curve, indices=(0, 0)))


def test_mapping_cache_order_context_and_precision(curve):
    ctx = curve.ctx
    implementation._mapped_curve_function_data.cache_clear(ctx)
    kleinian_sigma([.2], curve=POLYNOMIAL, ctx=ctx)
    info = implementation._mapped_curve_function_data.cache_info(ctx)
    kleinian_sigma([.3], curve=dict(reversed(list(POLYNOMIAL.items()))), ctx=ctx)
    assert implementation._mapped_curve_function_data.cache_info(ctx).misses == info.misses
    assert implementation._mapped_curve_function_data.cache_info(ctx).hits == info.hits + 1
    ctx.dps += 3
    kleinian_sigma([.2], curve=POLYNOMIAL, ctx=ctx)
    assert implementation._mapped_curve_function_data.cache_info(ctx).misses == info.misses + 1
    other = mp.clone()
    other.dps = ctx.dps
    implementation._mapped_curve_function_data.cache_clear(other)
    kleinian_sigma([.2], curve=POLYNOMIAL, ctx=other)
    assert implementation._mapped_curve_function_data.cache_info(other).misses == 1


def test_curve_cache_does_not_expose_mutable_matrices(curve, monkeypatch):
    value = kleinian_sigma([.2], curve=curve)
    first = curve.periods_kind_1()
    first.omega[0, 0] = 999
    def repeated(*args, **kwargs):
        raise AssertionError("point-independent setup should be cached")
    monkeypatch.setattr(curve, 'periods_kind_2', repeated)
    monkeypatch.setattr(curve, 'riemann_constant', repeated)
    assert curve.ctx.almosteq(kleinian_sigma([.2], curve=curve), value)


@pytest.mark.parametrize('kwargs, message', [
    ({'omega': [[1]]}, 'cannot be combined'),
    ({'characteristic': ([0], [0])}, 'cannot be combined'),
    ({'ctx': mp}, 'numerical context'),
])
def test_curve_input_conflicts(curve, kwargs, message):
    with pytest.raises(ValueError, match=message):
        kleinian_sigma([.2], curve=curve, **kwargs)


def test_missing_explicit_data_and_invalid_orders():
    with pytest.raises(ValueError, match='complete explicit period'):
        kleinian_sigma([.2])
    with pytest.raises(ValueError, match='sparse mapping'):
        kleinian_sigma([.2], curve=[1, 2, 3])
    with pytest.raises(ValueError, match='same size'):
        kleinian_sigma_normalization([[1]], [[1j, 0], [0, 1j]])
    with pytest.raises(ValueError, match='normalization must'):
        kleinian_sigma_normalization([[1]], [[1j]], normalization='unknown')


def test_custom_basis_normalization_is_theta_without_integration(curve, monkeypatch):
    custom = Curve(POLYNOMIAL, differentials_kind_1=[lambda x, y: 1 / y],
                   ctx=curve.ctx)
    def unnecessary(*args, **kwargs):
        raise AssertionError("theta multiplier requires no period integration")
    monkeypatch.setattr(custom, 'periods_kind_1', unnecessary)
    assert kleinian_sigma_normalization(curve=custom) == 1
    assert kleinian_sigma_normalization(curve={(0, 3): 1, (4, 0): -1, (0, 0): 1}) == 1
    assert kleinian_sigma_normalization(curve=curve, normalization='theta') == 1
    assert kleinian_sigma_normalization(curve=POLYNOMIAL, normalization='theta') == 1
    with pytest.raises(ValueError, match='cannot be combined'):
        kleinian_sigma_normalization(curve=curve, omega=[[1]], normalization='theta')
    with pytest.raises(ValueError, match='numerical context'):
        kleinian_sigma_normalization(curve=curve, ctx=mp, normalization='theta')
    # Normalization can be inspected without second-kind data; evaluation cannot.
    with pytest.raises(ValueError, match='second-kind basis'):
        kleinian_sigma([.2], curve=custom)


def test_ba_curve_route_matches_integral_route(curve):
    ctx = curve.ctx
    target = (2, ctx.sqrt(24))
    first = curve.periods_kind_1()
    second = curve.periods_kind_2()
    characteristic = curve.riemann_constant().characteristic
    abel = curve.abel_map_kind_1(target).value
    integral = curve.abel_map_kind_2(target).value
    expected = implementation.kleinian_baker_akhiezer(
        [.2], abel, integral, first.omega, first.tau, second.kappa,
        characteristic, ctx=ctx)
    assert ctx.almosteq(implementation.kleinian_baker_akhiezer(
        [.2], curve=curve, target=target), expected)
    assert ctx.almosteq(implementation.kleinian_baker_akhiezer(
        [.2], curve=POLYNOMIAL, target=target, ctx=ctx), expected)
    assert ctx.isfinite(implementation.kleinian_baker_akhiezer(
        [.2], curve=curve, target=target, reduce=True))


@pytest.mark.parametrize('kwargs, message', [
    ({'curve': POLYNOMIAL}, 'requires target'),
    ({'curve': POLYNOMIAL, 'target': 2, 'abel': [1]}, 'requires target'),
    ({'target': 2}, 'require curve'),
    ({'reduce': True}, 'require curve'),
    ({'curve': {(0, 3): 1, (4, 0): -1, (0, 0): 1}, 'target': 2}, 'odd-degree'),
    ({'curve': {(0, 2): 1, (4, 0): -1, (0, 0): 1}, 'target': 2}, 'odd-degree'),
])
def test_ba_curve_contract_errors(kwargs, message):
    with pytest.raises(ValueError, match=message):
        implementation.kleinian_baker_akhiezer([.2], **kwargs)


def test_custom_basis_uses_theta_scaling(curve, monkeypatch):
    # Reuse numerical records to isolate routing from geometric integration.
    first = curve.periods_kind_1()
    second = curve.periods_kind_2()
    constant = curve.riemann_constant()
    custom = Curve(POLYNOMIAL, differentials_kind_1=[lambda x, y: 1 / y],
                   differentials_kind_2=[lambda x, y: x / y], ctx=curve.ctx)
    monkeypatch.setattr(custom, 'periods_kind_1', lambda: first)
    monkeypatch.setattr(custom, 'periods_kind_2', lambda: second)
    monkeypatch.setattr(custom, 'riemann_constant', lambda: constant)
    raw = kleinian_sigma([.2], first.omega, first.tau, second.kappa,
                         constant.characteristic, ctx=curve.ctx)
    assert curve.ctx.almosteq(kleinian_sigma([.2], curve=custom), raw)
    with pytest.raises(ValueError, match='automatic hyperelliptic basis'):
        kleinian_sigma([.2], curve=custom, normalization='hyperelliptic')
    with pytest.raises(ValueError, match='odd-degree'):
        implementation.kleinian_baker_akhiezer([.2], curve=custom, target=(2, 1))


def test_curve_precision_changes_refresh_setup(curve):
    low = kleinian_sigma([.2], curve=curve)
    curve.ctx.dps += 5
    with pytest.warns(UserWarning, match='context changed'):
        high = kleinian_sigma([.2], curve=curve)
    assert abs(high - low) < curve.ctx.mpf('1e-21')
    high_from_map = kleinian_sigma([.2], curve=POLYNOMIAL, ctx=curve.ctx)
    assert curve.ctx.almosteq(high, high_from_map)


def test_normalization_helper_uses_no_second_kind_periods(curve, monkeypatch):
    def unnecessary(*args, **kwargs):
        raise AssertionError("the scalar needs no second-kind periods")
    monkeypatch.setattr(curve, 'periods_kind_2', unnecessary)
    constant = kleinian_sigma_normalization(curve=curve)
    assert curve.ctx.isfinite(constant)
    assert curve.ctx.almosteq(constant, kleinian_sigma_normalization(
        curve=POLYNOMIAL, ctx=curve.ctx))
