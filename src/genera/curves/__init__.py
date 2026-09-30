"""Numerical algebraic curves.

This module provides tools for computing with smooth plane algebraic curves,
including period matrices, Riemann constants, Abel maps, and integration on
the curve.

The primary interface is the Curve class, which provides lazy
evaluation of the computational pipeline: branch locus, monodromy, genus,
homology, periods, and Riemann constant.

Example:
    >>> from genera import algebraic_curve
    >>> from mpmath import mp
    >>> mp.dps = 30
    >>> # Fermat cubic x^3 + y^3 = 1
    >>> curve = algebraic_curve({(3, 0): 1, (0, 3): 1, (0, 0): -1})
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
    CurveFirstKindPeriods,
    CurveGenus,
    CurveHomology,
    CurveIntegral,
    CurveLatticeReduction,
    CurveMonodromy,
    CurvePath,
    CurvePlace,
    CurveRiemannConstant,
    CurveSecondKindAbelMap,
    CurveSecondKindPeriods,
    CurveValidation,
    algebraic_curve,
)

__all__ = [
    'Curve',
    'algebraic_curve',
    'CurveBranchLocus',
    'CurveChart',
    'CurveCheck',
    'CurveFirstKindPeriods',
    'CurveGenus',
    'CurveHomology',
    'CurveIntegral',
    'CurveLatticeReduction',
    'CurveMonodromy',
    'CurvePath',
    'CurvePlace',
    'CurveRiemannConstant',
    'CurveSecondKindAbelMap',
    'CurveSecondKindPeriods',
    'CurveValidation',
]
