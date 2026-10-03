"""Numerical algebraic curves.

This module provides tools for computing with smooth plane algebraic curves,
including period matrices, Riemann constants, Abel maps, and integration on
the curve.

The primary interface is the Curve class, which provides lazy
evaluation of the computational pipeline: branch locus, monodromy, genus,
homology, periods, and Riemann constant.

Example:
    >>> from genera import Curve
    >>> from mpmath import mp
    >>> mp.dps = 30
    >>> # Fermat cubic x^3 + y^3 = 1
    >>> curve = Curve({(3, 0): 1, (0, 3): 1, (0, 0): -1})
    >>> curve.genus
    1
    >>> curve.branch_locus.degree
    3
"""

# Public namedtuples for result records
from .algebraic_curve import (
    Curve,
    CurveBranchLocus,
    CurveChart,
    CurveCheck,
    CurvePeriodsKind1,
    CurveGenus,
    CurveHomology,
    CurveIntegral,
    CurveLatticeReduction,
    CurveMonodromy,
    CurvePath,
    CurvePlace,
    CurveRiemannConstant,
    CurveAbelMapKind1,
    CurveAbelMapKind2,
    CurvePeriodsKind2,
    CurveValidation,
)

__all__ = [
    'Curve',
    'CurveBranchLocus',
    'CurveChart',
    'CurveCheck',
    'CurvePeriodsKind1',
    'CurveGenus',
    'CurveHomology',
    'CurveIntegral',
    'CurveLatticeReduction',
    'CurveMonodromy',
    'CurvePath',
    'CurvePlace',
    'CurveRiemannConstant',
    'CurveAbelMapKind1',
    'CurveAbelMapKind2',
    'CurvePeriodsKind2',
    'CurveValidation',
]
