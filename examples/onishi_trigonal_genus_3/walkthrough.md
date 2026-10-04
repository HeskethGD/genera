# Verifying Lemma 5.1 of "Abelian functions for trigonal curves of genus three"

`lemma_5_1_demo.py` numerically verifies the fourteen four-index relations
of Lemma 5.1 of J. C. Eilbeck, V. Z. Enolskii, S. Matsutani, Y. Onishi,
E. Previato, *Abelian functions for trigonal curves of genus three*
(arXiv:math/0610019v2; the LaTeX source lives in `latex/`), against
Generapy's general plane-curve and Kleinian machinery.  The curves and the
evaluation points are chosen for this example alone; nothing is taken
from Bernatska's or any other example, and no period, characteristic or
function value is taken from the paper.

## The lemma

For the general (3,4)-curve

    f(x, y) = y^3 + (mu1 x + mu4) y^2 + (mu2 x^2 + mu5 x + mu8) y
              - (x^4 + mu3 x^3 + mu6 x^2 + mu9 x + mu12),

with the holomorphic basis (dx/f_y, x dx/f_y, y dx/f_y) and the standard
sigma function of the paper, Lemma 5.1 states fourteen linear relations

    wp_{ijkell} = polynomial in wp_{ij}, wp_{13}, wp_{1333} and the mu_j,

covering every four-index function except wp_{1333}, which is kept as an
independent generator.  Every relation was checked to be homogeneous of
weight -(w_i+w_j+w_k+w_ell) under the weights w(u1, u2, u3) = (5, 2, 1),
w(mu_j) = j before implementation, and each right-hand side is
transcribed verbatim in `lemma_relations()`.

## Pipeline

For each curve the demo runs:

1. `curve_monodromy` -- genus 3, transitive permutation system,
   product identity (a nonsingularity check).
2. `curve_periods` with Onishi's holomorphic basis; for the purely
   trigonal curve also with the corrected second-kind forms (below),
   which yields the symmetric quadratic form kappa directly.
3. `curve_riemann_constant` with a chart-backed place at the unique
   point at infinity (monomial chart x = t^-3, y = t^-4 w); the
   characteristic must come out half-integral.
4. Sigma sanity checks: sigma is odd, and sigma vanishes on the Abel
   image (based at infinity) of a degree-two divisor.  Both pin down the
   theta characteristic independently of the relations.
5. One batched `kleinian_p` call per generic point for every two-index
   function, wp_1333 and all fourteen four-index left-hand sides; the
   residuals of the fourteen relations are reported over three generic
   points.

Convention mapping: Generapy's `curve_periods` returns omega as half the
a-period matrix and eta as minus one half of the a-periods of supplied
second-kind forms, so kappa = eta*omega^-1 = -eta'*(omega')^-1 is exactly
the quadratic form `kleinian_sigma` wants.  Zero-based kleinian indices
(i-1 for the paper's i) follow from the differential order.

## Results (dps 30, run time about 2.5 minutes)

    general curve (all moduli non-zero):
        Riemann characteristic  a = (0, 0, 1/2),  b = (0, 0, 1/2)
        sigma parity residual                 0.0
        sigma divisor residual              1.7e-16
        kappa closure residual              2.9e-15
        worst relation residual (relative)  1.5e-16

    purely trigonal curve y^3 = x^4 + 2x^3 + 3x^2 + 5x + 7:
        Riemann characteristic  a = (1/2, 0, 0),  b = (1/2, 0, 0)
        sigma parity residual                 0.0
        sigma divisor residual              7.8e-11
        worst relation residual (relative)  3.1e-11

All fourteen relations hold on both curves.  On the purely trigonal
curve nothing is fitted, so this is a complete, independent verification
of both Generapy's engine and the relations.  On the general curve the six
entries of the symmetric quadratic form are closed against the fourteen
relations at the first evaluation point (fourteen complex equations in
six complex unknowns; the closure residual 2.9e-15 is itself a
nontrivial consistency check), and all relations are then verified at
two further points.

## Findings about the paper

Working through the lemma surfaced several internal inconsistencies of
arXiv:math/0610019v2.  None of them affect Lemma 5.1, which verifies
cleanly, but they matter for anyone reproducing the paper's sigma
normalization:

1. **The displayed h_j triple (Section 3) has sign errors.**  For the
   purely trigonal specialization, comparing the paper's Appendix-A
   bi-differential with the Section-3 realization
   Omega = (d Sigma + sum_k omega_k(z) eta_k(x)) dx dz determines
        h_1 = (5 x^2 + 3 mu3 x + mu6) y,   h_2 = 2 x y,   h_3 = x^2,
   to roundoff (verified in development at 1e-31 relative).  These are
   the negatives of the displayed h_1, h_2, h_3 = -(5x^2 - 3mu3 x +
   mu6) y, -2xy, -x^2, i.e. the display additionally flips the sign of
   the mu3 x term relative to its other terms.  The displayed triple
   fails the Appendix-A realization by O(1) and produces an
   asymmetric (hence invalid) kappa.

2. **The Appendix-A polynomial F does not work for general moduli, at
   least as we read it.**  Our transcription passes the symmetry and
   diagonal F((x,y),(x,y)) = f_y^2 checks in the purely trigonal
   specialization (where it reproduces the paper's own Sigma) but
   violates the diagonal condition when mu1, mu4 are non-zero, so the
   general-curve h_j cannot be recovered from it.  We could not rule
   out a misreading of the F display.

3. **The Section-7 sigma expansion disagrees with the Lemma 5.1
   normalization.**  On the purely trigonal curve the expansion
   reduces to a few hand-checkable terms (C_5, one C_8 term, four C_11
   terms); comparing theta with the series and subtracting the
   second-kind quadratic form leaves a t^2 drift (about 0.013 between
   t = 0.1 and t = 0.2 on the u3 axis, reported by the demo) instead
   of the t^16-suppressed constant expected from the omitted C_13+.
   Since no transcription error can mimic a pure quadratic-form shift
   there, the expansion as displayed corresponds to sigma multiplied
   by a non-trivial exp(transp(u) A u / 2).

4. **Precision note (Generapy side).**  `curve_riemann_constant` on the
   general curve (whose monodromy paths have minimum branch clearance
   3e-4) returned an O(1)-wrong value at dps 20 -- silently, since the
   wrong value still reduced to half-integral characteristic entries.
   At dps 25 and 30 the computation is exact and is confirmed by both
   the sigma-divisor test and the relations.  Use dps >= 25; the demo
   defaults to 30 and its sigma checks would catch a recurrence.

For these reasons the demo's two tracks differ: the purely trigonal
curve is verified parameter-free using the Appendix-A-derived
second-kind forms, while the general curve closes the six kappa
entries against the relations themselves and verifies everything else.

## Running

From the Generapy repository root:

    .venv/bin/python -m examples.onishi_trigonal_genus_3.lemma_5_1_demo
    .venv/bin/python -m examples.onishi_trigonal_genus_3.lemma_5_1_demo \
        --dps 35 --curves purely-trigonal

The script exits non-zero and lists the failing stages if any check
fails.  Generapy and its mpmath dependency are sufficient.