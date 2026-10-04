# Generapy

<img src="https://raw.githubusercontent.com/HeskethGD/generapy/main/docs/_static/generapy-logo.png" alt="Generapy logo: a curved, gridded letter g in teal and gold" width="220">

Generapy is an early-stage numeric Python package for arbitrary-precision computation
with algebraic curves and Abelian functions, particularly for applications in
integrable systems.

Generapy uses [mpmath](https://mpmath.org/) for arbitrary-precision arithmetic.
It is an independent project and does not modify the mpmath namespace.
Generapy retains the applicable BSD-3-Clause copyright and licence notice.

## Installation

Install from PyPI:

```sh
pip install generapy
```

Requires Python 3.11 or later.

**Documentation:** https://generapy.readthedocs.io/en/latest/

## Quick Example

For example, verify Baker's genus-two addition formula using Kleinian sigma
and P functions at arbitrary Abelian arguments. The identity is from
H. F. Baker, *An Introduction to the Theory of Multiply Periodic Functions*
(Cambridge University Press, 1907), Chapter V, p. 100
([digitized edition](https://books.google.com/books?id=0EQLAAAAYAAJ&pg=PA100)):

```python
from mpmath import mp

from generapy import Curve
from generapy import kleinian_sigma as ks
from generapy import kleinian_p as kp

mp.dps = 30

# y² = 4x(x² - 1)(x² - 4) = 4x⁵ - 20x³ + 16x
polynomial = {(0, 2): 1, (5, 0): -4, (3, 0): 20, (1, 0): -16}
curve = Curve(polynomial)

u = mp.matrix([mp.mpf("0.3"), mp.mpf("0.2")])
v = mp.matrix([mp.mpf("0.1"), mp.mpf("0.4")])

lhs = ks(u + v, curve=curve) * ks(u - v, curve=curve) / (
    ks(u, curve=curve)**2 * ks(v, curve=curve)**2)

# Baker's indices are one-based; Generapy's indices are zero-based.
p11_u, p12_u, p22_u = kp(u, curve=curve, indices=((0, 0), (0, 1), (1, 1)))
p11_v, p12_v, p22_v = kp(v, curve=curve, indices=((0, 0), (0, 1), (1, 1)))
rhs = p22_u * p12_v - p12_u * p22_v + p11_v - p11_u

print(mp.nstr(abs(rhs - lhs), 6))  # 6.31089e-29
```

The curve supplies the period data and canonical hyperelliptic sigma
normalization automatically. See the [example walkthrough](https://generapy.readthedocs.io/en/latest/examples/baker_formula.html)
for the identity and conventions, or run the
[complete example](https://github.com/HeskethGD/generapy/blob/v0.1.0/examples/baker_formula/baker_formula_demo.py).

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
.venv/bin/pytest --cov=generapy --cov-report=term-missing
.venv/bin/sphinx-build -W -b html docs build/sphinx/html
.venv/bin/python -m build
```

Further numerical APIs and documentation will be added incrementally.

### Publishing releases

Publishing is manual through GitHub Actions. Push a release tag, then open
**Actions → Publish to PyPI → Run workflow**, select `main`, and enter the
tag (for example, `v0.1.0`). Pushing a commit or tag does not trigger publishing.
The workflow tests the tagged commit on supported Python versions, builds and
validates the distributions, checks an installed wheel, and publishes to PyPI.

For the one-time setup, create a GitHub environment named `pypi` and register
a [PyPI Trusted Publisher](https://docs.pypi.org/trusted-publishers/adding-a-publisher/)
for owner `HeskethGD`, repository `generapy`, workflow `release.yml`, and
environment `pypi`. For the first release, register a
[pending publisher](https://docs.pypi.org/trusted-publishers/creating-a-project-through-oidc/)
for project `generapy` in your PyPI account's Publishing settings. No API token
is needed. Publish each version once; use a new tag and version for changes
after publication.

Before the first release, rename the GitHub repository to `generapy` and
configure a Read the Docs project with the slug `generapy` to match the links
and publisher configuration above.
Commit the package rename and move the unpublished `v0.1.0` tag to that commit
before running the release workflow.

## Examples

Literature-based numerical examples live in the repository's `examples/`
directory. They are not installed as part of the `generapy` package and do
not require plotting libraries. Run them from a source checkout, for example:

```sh
.venv/bin/python -m examples.manakov.manakov_kleinian_demo
```

The dynamics examples compare their Abelian-function solutions against a
shared arbitrary-precision Runge--Kutta implementation. Lightweight versions
of the examples run in the test suite to keep the published workflows valid.
