"""Baker's genus-two addition formula at arbitrary Abelian arguments.

Run with: python -m examples.baker_formula.baker_formula_demo
"""

from mpmath import mp

from generapy import Curve
from generapy import kleinian_p as kp
from generapy import kleinian_sigma as ks


def documentation_example():
    """Return the two sides and their residual, evaluated at 30 digits."""
    with mp.workdps(30):
        # y² = 4x(x² - 1)(x² - 4) = 4x⁵ - 20x³ + 16x
        polynomial = {(0, 2): 1, (5, 0): -4, (3, 0): 20, (1, 0): -16}
        curve = Curve(polynomial)

        # Arbitrary vectors in C²; no points on the curve or Abel maps needed.
        u = mp.matrix([mp.mpf("0.3"), mp.mpf("0.2")])
        v = mp.matrix([mp.mpf("0.1"), mp.mpf("0.4")])

        lhs = ks(u + v, curve=curve) * ks(u - v, curve=curve) / (
            ks(u, curve=curve)**2 * ks(v, curve=curve)**2)

        # Baker's indices are one-based; Generapy's indices are zero-based.
        p11_u, p12_u, p22_u = kp(u, curve=curve, indices=((0, 0), (0, 1), (1, 1)))
        p11_v, p12_v, p22_v = kp(v, curve=curve, indices=((0, 0), (0, 1), (1, 1)))
        rhs = p22_u * p12_v - p12_u * p22_v + p11_v - p11_u

        residual = rhs - lhs
        assert abs(residual) < mp.mpf("1e-24"), "Baker addition formula failed"
        return lhs, rhs, residual


if __name__ == "__main__":
    lhs, rhs, residual = documentation_example()
    print("LHS:", mp.nstr(lhs, 12))
    print("RHS:", mp.nstr(rhs, 12))
    print("|RHS - LHS|:", mp.nstr(abs(residual), 6))
