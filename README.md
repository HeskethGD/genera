# Genera

<img src="https://raw.githubusercontent.com/HeskethGD/genera/main/docs/_static/genera-logo.png" alt="Genera logo: a curved, gridded letter g in teal and gold" width="220">

Genera is an early-stage numeric Python package for arbitrary-precision computation
with algebraic curves and Abelian functions, particularly for applications in
integrable systems.

Genera uses [mpmath](https://mpmath.org/) for arbitrary-precision arithmetic.
It is an independent project and does not modify the mpmath namespace.
Genera retains the applicable BSD-3-Clause copyright and licence notice.

## Installation

Install from PyPI:

```sh
pip install genera
```

Requires Python 3.11 or later.

**Documentation:** https://genera.readthedocs.io/en/latest/

## Quick Example

For example, verify Baker's genus-two addition formula using Kleinian sigma
and P functions at arbitrary Abelian arguments. The identity is from
H. F. Baker, *An Introduction to the Theory of Multiply Periodic Functions*
(Cambridge University Press, 1907), Chapter V, p. 100
([digitized edition](https://books.google.com/books?id=0EQLAAAAYAAJ&pg=PA100)):

```python
from mpmath import mp

from genera import Curve
from genera import kleinian_sigma as ks
from genera import kleinian_p as kp

mp.dps = 30

# y² = 4x(x² - 1)(x² - 4) = 4x⁵ - 20x³ + 16x
polynomial = {(0, 2): 1, (5, 0): -4, (3, 0): 20, (1, 0): -16}
curve = Curve(polynomial)

u = mp.matrix([mp.mpf("0.3"), mp.mpf("0.2")])
v = mp.matrix([mp.mpf("0.1"), mp.mpf("0.4")])

lhs = ks(u + v, curve=curve) * ks(u - v, curve=curve) / (
    ks(u, curve=curve)**2 * ks(v, curve=curve)**2)

# Baker's indices are one-based; Genera's indices are zero-based.
p11_u, p12_u, p22_u = kp(u, curve=curve, indices=((0, 0), (0, 1), (1, 1)))
p11_v, p12_v, p22_v = kp(v, curve=curve, indices=((0, 0), (0, 1), (1, 1)))
rhs = p22_u * p12_v - p12_u * p22_v + p11_v - p11_u

print(mp.nstr(abs(rhs - lhs), 6))  # 6.31089e-29
```

The curve supplies the period data and canonical hyperelliptic sigma
normalization automatically. See the [example walkthrough](https://genera.readthedocs.io/en/latest/examples/baker_formula.html)
for the identity and conventions, or run the
[complete example](https://github.com/HeskethGD/genera/blob/v0.1.0/examples/baker_formula/baker_formula_demo.py).

## Development

Create an isolated environment and install the development dependencies:

```sh
python -m venv .venv
.venv/bin/python -m pip install -e '.[develop]'
```

Run the local checks with:

```sh
.venv/bin/ruff check .
.venv/bin/pytest
.venv/bin/pytest --cov=genera --cov-report=term-missing
.venv/bin/sphinx-build -W -b html docs build/sphinx/html
.venv/bin/python -m build
```

Further numerical APIs and documentation will be added incrementally.

## Examples

Literature-based numerical examples live in the repository's `examples/`
directory. They are not installed as part of the `genera` package and do
not require plotting libraries. Run them from a source checkout, for example:

```sh
.venv/bin/python -m examples.manakov.manakov_kleinian_demo
```

The dynamics examples compare their Abelian-function solutions against a
shared arbitrary-precision Runge--Kutta implementation. Lightweight versions
of the examples run in the test suite to keep the published workflows valid.
