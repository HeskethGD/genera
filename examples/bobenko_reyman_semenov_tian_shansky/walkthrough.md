# Bobenko–Reyman–Semenov-Tian-Shansky genus-three top

This example evaluates the genus-three theta-functional solution of the
classical Kowalewski top given by A. I. Bobenko, A. G. Reyman and
M. A. Semenov-Tian-Shansky in *The Kowalewski top 99 years later: a Lax pair,
generalizations and explicit solutions*.

It intentionally contains only the working genus-three construction. The
separate experimental genus-two degeneration has not been migrated.

## Formula and checks

For the numerical period matrix and marked shifts in the script, equation
(7.42) gives the angular-momentum variables `(ell1, ell2, ell3)` and gravity
variables `(g1, g2, g3)` as quotients of genus-three Riemann theta functions
and their directional derivatives.

The implementation checks:

- the two cases of the theta decomposition in equation (7.64);
- constancy of the quotient in equation (7.60);
- the four algebraic first integrals; and
- all six physical components against an independent Euler–Poisson RK4
  trajectory.

The RK4 implementation is shared with the other dynamics examples through
`examples._rk4`; it does not reuse any part of the theta formula.

## Run

From the Generapy repository root:

```bash
.venv/bin/python -m \
    examples.bobenko_reyman_semenov_tian_shansky.kowalewski_genus_three
```

Use `--stop`, `--steps` and `--samples` to change the comparison interval and
resolution. `steps` must be divisible by `samples - 1`. The script exits with
an error if the component errors, imaginary parts or invariant drift exceed
`--tol`.

The checked phase lies on a numerically real component of the invariant torus.
At 30 decimal digits, a short comparison with 400 RK4 steps gives component
errors around `1e-12`, with imaginary parts and invariant drift of the same
order.
