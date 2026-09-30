"""Numerical algebraic curves and Abelian functions."""

from .kleinian import (
    kleinian_baker_akhiezer, kleinian_p, kleinian_sigma,
    kleinian_sigma_jet, kleinian_zeta,
)
from .riemann_theta import rtheta, rtheta_jet

try:
    from ._version import version as __version__
except ImportError:  # pragma: no cover - generated in builds by setuptools-scm
    __version__ = "0+unknown"

__all__ = [
    "kleinian_baker_akhiezer", "kleinian_p", "kleinian_sigma",
    "kleinian_sigma_jet", "kleinian_zeta", "rtheta", "rtheta_jet",
]
