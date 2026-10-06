"""Orchestration for the specialized hyperelliptic implementation."""

from ..._context import ctx_lru_cache

from .integration import (
    _branch_target_integrals,
    _hyperelliptic_intervals,
    _second_kind_interval,
)
from .jacobian import (
    _abel_lattice_shift,
    _branch_abel_values,
    _branch_second_kind_values,
    _first_kind_periods,
    _second_kind_periods,
    _symmetrize_period_matrix,
    _validate_legendre_relation,
)
from .model import (
    _admissible_branch_vertex,
    _evaluate_polynomial,
    _hyperelliptic_coefficients,
    _normalize_abel_targets,
    _prepare_hyperelliptic_curve,
    _target_branch_index,
)


@ctx_lru_cache(maxsize=8)
def _hyperelliptic_curve_data(ctx, coefficients, method):
    """Reuse immutable root preparation without computing global integrals."""
    return _prepare_hyperelliptic_curve(ctx, coefficients, method)


def _validate_abel_targets(ctx, coefficients, targets, target_eps):
    for x, y in targets:
        curve_value = _evaluate_polynomial(ctx, coefficients, x)
        evaluation_scale = max(
            ctx.one, abs(y ** 2),
            ctx.fsum(abs(coefficient * x ** degree)
                     for degree, coefficient in enumerate(coefficients)))
        tolerance = 100 * target_eps * evaluation_scale
        if not ctx.almosteq(
                y ** 2, curve_value,
                rel_eps=100 * target_eps, abs_eps=tolerance):
            raise ValueError("target point must satisfy y**2 = P(x)")


def _hyperelliptic_local_abel_difference(ctx, coefficients, targets, base):
    """Cancel a shared Baker branch reference before doing quadrature.

    Return None when the marked references differ, so the caller can use the
    full absolute maps. Subtract the original branch-to-point paths rather
    than introducing a new path that could change the period convention.
    """
    target_eps = +ctx.eps
    with ctx.extraprec(20):
        coefficients = _hyperelliptic_coefficients(ctx, coefficients)
        coefficients, roots, _, _, genus, _ = _hyperelliptic_curve_data(
            ctx, coefficients, "auto")
        points = tuple(targets) + (base,)
        _validate_abel_targets(ctx, coefficients, points, target_eps)
        anchors = []
        at_branch = []
        for x, y in points:
            index = _target_branch_index(ctx, x, y, roots, target_eps)
            at_branch.append(index is not None)
            anchors.append(index if index is not None else
                           _admissible_branch_vertex(ctx, x, roots, target_eps))
        if any(index != anchors[-1] for index in anchors[:-1]):
            return None
        result = ctx.zeros(genus, 1)
        for position, ((x, y), index, branch) in enumerate(
                zip(points, anchors, at_branch)):
            if branch:
                continue
            values = _branch_target_integrals(
                ctx, roots, coefficients[-1], index, x, y, genus, target_eps)
            multiplier = -len(targets) if position == len(targets) else 1
            for row in range(genus):
                result[row] += multiplier * values[row]
    return +result


@ctx_lru_cache(maxsize=8)
def _hyperelliptic_abel_data(ctx, coefficients, method, second_kind):
    """Cache immutable curve and branch data independently of the target.

    The caller enters guard precision before this stage. The context cache
    separates first- and second-kind guards, precision, and curve coefficients;
    each Abel evaluation still validates targets and constructs fresh matrices.
    """
    curve_data = _hyperelliptic_curve_data(ctx, coefficients, method)
    coefficients, roots, unused_tolerance, use_real, genus, even_degree = curve_data
    monomial_count = 2 * genus + 1 if second_kind else genus
    intervals, b_sign = _hyperelliptic_intervals(
        ctx, coefficients, roots, use_real, monomial_count)
    intervals = tuple(tuple(values) for values in intervals)
    branch_values = _branch_abel_values(
        ctx, roots, intervals, coefficients[-1], genus, even_degree)
    second_coefficients = coefficients
    if second_kind and not even_degree:
        second_coefficients += (ctx.zero,)
    second_branch_values = (_branch_second_kind_values(
        ctx, second_coefficients, roots, intervals, genus, even_degree)
        if second_kind else None)
    return (curve_data, intervals, b_sign, branch_values,
            second_coefficients, second_branch_values)


def _hyperelliptic_periods(ctx, coefficients, method="auto",
                           second_kind=False):
    """Compute Baker-marked hyperelliptic period matrices."""
    target_eps = +ctx.eps
    quadrature_guard = 20
    cancellation_guard = 40 if second_kind else 0
    with ctx.extraprec(quadrature_guard + cancellation_guard):
        curve_data = _prepare_hyperelliptic_curve(
            ctx, coefficients, method)
        (coefficients, roots, _unused_root_tolerance, use_real_method,
         genus, even_degree) = curve_data
        monomial_count = 2 * genus + 1 if second_kind else genus
        intervals, b_sign = _hyperelliptic_intervals(
            ctx, coefficients, roots, use_real_method, monomial_count)
        omega, omega_prime, tau, inverse_omega = _first_kind_periods(
            ctx, intervals, genus, even_degree, b_sign, target_eps)
        if second_kind:
            second_coefficients = coefficients
            if not even_degree:
                second_coefficients += (ctx.zero,)
            eta, eta_prime = _second_kind_periods(
                ctx, second_coefficients, intervals, genus, even_degree,
                b_sign)
            kappa = eta * inverse_omega
            _symmetrize_period_matrix(
                ctx, kappa, target_eps, "kappa matrix")
            _validate_legendre_relation(
                ctx, omega, omega_prime, eta, eta_prime, target_eps)
    if second_kind:
        return (+omega, +omega_prime, +eta, +eta_prime,
                +tau, +kappa)
    return +omega, +omega_prime, +tau


def _hyperelliptic_abel_map(
        ctx, coefficients, target, method="auto", reduce=False,
        second_kind=False, _return_shift=False):
    """Evaluate Baker-marked first- and optionally second-kind integrals."""
    target_eps = +ctx.eps
    quadrature_guard = 20
    cancellation_guard = 40 if second_kind else 0
    with ctx.extraprec(quadrature_guard + cancellation_guard):
        targets = _normalize_abel_targets(ctx, target)
        coefficients = _hyperelliptic_coefficients(ctx, coefficients)
        (curve_data, intervals, b_sign, branch_values,
         second_coefficients, second_branch_values) = _hyperelliptic_abel_data(
             ctx, coefficients, method, second_kind)
        (coefficients, roots, _unused_root_tolerance, use_real_method,
         genus, even_degree) = curve_data
        _validate_abel_targets(ctx, coefficients, targets, target_eps)

        monomial_count = 2 * genus + 1 if second_kind else genus
        result = ctx.zeros(genus, 1)
        if second_kind:
            second_result = ctx.zeros(genus, 1)
        for x, y in targets:
            branch_index = _target_branch_index(
                ctx, x, y, roots, target_eps)
            if branch_index is None:
                branch_index = _admissible_branch_vertex(
                    ctx, x, roots, target_eps)
                final_monomials = _branch_target_integrals(
                    ctx, roots, coefficients[-1], branch_index, x, y,
                    monomial_count, target_eps)
            else:
                final_monomials = (ctx.zero,) * monomial_count
            for row in range(genus):
                result[row] += (
                    branch_values[branch_index][row]
                    + final_monomials[row])
                if second_kind:
                    second_result[row] += (
                        second_branch_values[branch_index][row]
                        + _second_kind_interval(
                            ctx, second_coefficients, final_monomials, row,
                            genus))

        if reduce:
            omega, omega_prime, _unused_tau, _unused_inverse = (
                _first_kind_periods(
                    ctx, intervals, genus, even_degree, b_sign, target_eps))
            periods, lattice_shift = _abel_lattice_shift(
                ctx, result, omega, omega_prime, target_eps)
            result -= periods * lattice_shift
            if second_kind:
                eta, eta_prime = _second_kind_periods(
                    ctx, second_coefficients, intervals, genus, even_degree,
                    b_sign)
                second_periods = ctx.matrix(genus, 2 * genus)
                second_periods[:, :genus] = 2 * eta
                second_periods[:, genus:] = 2 * eta_prime
                second_result += second_periods * lattice_shift
    if second_kind and _return_shift:
        shift = (tuple(int(value) for value in lattice_shift)
                 if reduce else None)
        return +result, +second_result, shift
    if second_kind:
        return +result, +second_result
    if _return_shift:
        shift = tuple(int(value) for value in lattice_shift) if reduce else None
        return +result, shift
    return +result
