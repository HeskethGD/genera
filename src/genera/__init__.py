"""Numerical algebraic curves and Abelian functions."""

from .riemann_theta import rtheta, rtheta_jet

try:
    from ._version import version as __version__
except ImportError:  # pragma: no cover - generated in builds by setuptools-scm
    __version__ = "0+unknown"

__all__ = ["rtheta", "rtheta_jet"]
