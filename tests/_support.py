"""Shared helpers for adapting explicit numerical contexts in tests."""

from genera import Curve


def make_curve(ctx, polynomial, **kwargs):
    """Construct a curve through Genera's public context-aware constructor."""
    return Curve(polynomial, ctx=ctx, **kwargs)


def with_basis(curve, *, differentials_kind_1=None, differentials_kind_2=None):
    """Construct an independent curve with the supplied coordinate convention."""
    return Curve(curve.polynomial, ctx=curve.ctx, differentials_kind_1=differentials_kind_1,
                 differentials_kind_2=differentials_kind_2)
