"""Numerical algebraic curves and Abelian functions."""

from .curves import (
    Curve, CurveBranchLocus, CurveChart, CurveCheck, CurvePeriodsKind1,
    CurveGenus, CurveHomology, CurveIntegral, CurveLatticeReduction,
    CurveMonodromy, CurvePath, CurvePlace, CurveRiemannConstant,
    CurveAbelMapKind1, CurveAbelMapKind2, CurvePeriodsKind2, CurveValidation,
)
from .kleinian import (
    kleinian_p, kleinian_sigma, kleinian_sigma_jet, kleinian_sigma_normalization, kleinian_zeta,
)
from .riemann_theta import rtheta, rtheta_jet

try:
    from ._version import version as __version__
except ImportError:  # pragma: no cover - generated in builds by setuptools-scm
    __version__ = "0+unknown"

__all__ = [
    "Curve", "CurveBranchLocus", "CurveChart", "CurveCheck",
    "CurvePeriodsKind1", "CurveGenus", "CurveHomology", "CurveIntegral",
    "CurveLatticeReduction", "CurveMonodromy", "CurvePath", "CurvePlace",
    "CurveRiemannConstant", "CurveAbelMapKind1", "CurveAbelMapKind2", "CurvePeriodsKind2",
    "CurveValidation", "kleinian_p", "kleinian_sigma",
    "kleinian_sigma_jet", "kleinian_sigma_normalization", "kleinian_zeta",
    "rtheta", "rtheta_jet",
]
