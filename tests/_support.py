"""Shared helpers for adapting explicit numerical contexts in tests."""

from genera import algebraic_curve


def make_curve(ctx, specification, **kwargs):
    """Construct a curve through Genera's public context-aware factory."""
    return algebraic_curve(specification, ctx=ctx, **kwargs)
