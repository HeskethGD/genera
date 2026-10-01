# Riemann theta baseline: mpmath `rtheta` vs Wolfram `SiegelTheta`

This compares the current generic mpmath direct sum with Wolfram Engine's multidimensional `SiegelTheta`. It includes genus 1, 2 and 3, zero and nonzero half-characteristics, diagonal factorisation controls, and genuinely coupled Riemann matrices.

Timed rows are warmed fixed-tau evaluations over varying z vectors. Wolfram process startup and result formatting are excluded.
Derivative rows are correctness-only because Wolfram 14.3 computes them here through `NumericalCalculus`ND`, not a directly comparable native `SiegelTheta` derivative call.

- Generated: `2026-10-01T16:44:17.825129+01:00`
- mpmath revision: `c311b91 (dirty)`
- Python: `3.13.4`
- Platform: `macOS-27.0.1-arm64-arm-64bit-Mach-O`
- Decimal precisions: `50, 100`
- mpmath timing: median of `3` repeats, `1` batch per repeat
- Wolfram timing: `RepeatedTiming` target `0.2` seconds, `1` batch
- Arguments per batch: up to `4`

## Summary

Across all timed rows, the geometric-mean `mpmath / Wolfram` time ratio was **0.40x**. The range was **0.03x to 1.23x**.

The largest scaled numerical difference was `6.34354e-50`.

- At 50 dps: geometric-mean ratio **0.42x**; largest scaled difference `6.34354e-50`.
- At 100 dps: geometric-mean ratio **0.38x**; largest scaled difference `2.66322e-100`.

## Timing

| dps | Case | Genus | Characteristic | mpmath (ms/call) | Wolfram (ms/call) | mpmath / Wolfram |
|---:|---|---:|---|---:|---:|---:|
| 50 | g1 | 1 | zero | 0.258 | 0.401 | 0.64x |
| 50 | g1 | 1 | half | 0.260 | 0.428 | 0.61x |
| 50 | g2-diagonal | 2 | zero | 2.027 | 1.886 | 1.07x |
| 50 | g2-diagonal | 2 | half | 2.300 | 2.230 | 1.03x |
| 50 | g2-coupled | 2 | zero | 2.288 | 2.499 | 0.92x |
| 50 | g2-coupled | 2 | half | 2.300 | 2.770 | 0.83x |
| 50 | g2-coupled | 2 | dz0 | — | — | — |
| 50 | g2-coupled | 2 | half-dz0 | — | — | — |
| 50 | g3-diagonal | 3 | zero | 19.817 | 16.125 | 1.23x |
| 50 | g3-diagonal | 3 | half | 20.810 | 18.757 | 1.11x |
| 50 | g3-coupled | 3 | zero | 24.295 | 23.286 | 1.04x |
| 50 | g3-coupled | 3 | half | 24.732 | 25.902 | 0.95x |
| 50 | g3-coupled | 3 | dz0 | — | — | — |
| 50 | g3-coupled | 3 | half-dz0 | — | — | — |
| 50 | g2-near-boundary | 2 | zero | 2.286 | 2.372 | 0.96x |
| 50 | g2-near-boundary | 2 | half | 2.340 | 2.736 | 0.86x |
| 50 | g2-unreduced | 2 | zero | 1.848 | 13.833 | 0.13x |
| 50 | g2-unreduced | 2 | half | 1.921 | 14.183 | 0.14x |
| 50 | g2-unreduced | 2 | real | 2.249 | 14.258 | 0.16x |
| 50 | g3-unreduced | 3 | zero | 12.535 | 331.894 | 0.04x |
| 50 | g3-unreduced | 3 | half | 12.069 | 373.616 | 0.03x |
| 50 | g3-unreduced | 3 | real | 12.637 | 376.124 | 0.03x |
| 50 | g3-block | 3 | zero | 22.918 | 24.050 | 0.95x |
| 50 | g3-block | 3 | half | 23.193 | 27.100 | 0.86x |
| 100 | g1 | 1 | zero | 0.347 | 0.628 | 0.55x |
| 100 | g1 | 1 | half | 0.353 | 0.644 | 0.55x |
| 100 | g2-diagonal | 2 | zero | 3.866 | 3.978 | 0.97x |
| 100 | g2-diagonal | 2 | half | 3.836 | 4.439 | 0.86x |
| 100 | g2-coupled | 2 | zero | 4.227 | 4.499 | 0.94x |
| 100 | g2-coupled | 2 | half | 4.223 | 4.980 | 0.85x |
| 100 | g2-coupled | 2 | dz0 | — | — | — |
| 100 | g2-coupled | 2 | half-dz0 | — | — | — |
| 100 | g3-diagonal | 3 | zero | 49.537 | 42.379 | 1.17x |
| 100 | g3-diagonal | 3 | half | 50.418 | 48.580 | 1.04x |
| 100 | g3-coupled | 3 | zero | 59.642 | 61.449 | 0.97x |
| 100 | g3-coupled | 3 | half | 60.471 | 68.545 | 0.88x |
| 100 | g2-near-boundary | 2 | zero | 4.123 | 4.561 | 0.90x |
| 100 | g2-near-boundary | 2 | half | 4.126 | 4.773 | 0.86x |
| 100 | g2-unreduced | 2 | zero | 2.902 | 25.226 | 0.12x |
| 100 | g2-unreduced | 2 | half | 3.005 | 25.912 | 0.12x |
| 100 | g2-unreduced | 2 | real | 2.953 | 26.306 | 0.11x |
| 100 | g3-unreduced | 3 | zero | 28.839 | 907.203 | 0.03x |
| 100 | g3-unreduced | 3 | half | 28.978 | 1028.684 | 0.03x |
| 100 | g3-unreduced | 3 | real | 29.433 | 1024.042 | 0.03x |
| 100 | g3-block | 3 | zero | 56.723 | 72.873 | 0.78x |
| 100 | g3-block | 3 | half | 57.880 | 71.102 | 0.81x |

A ratio above 1 means Wolfram was faster; below 1 means mpmath was faster. These are steady-state call costs, not command-line latency.

## Numerical agreement

Error is `abs(mpmath - Wolfram) / max(1, abs(Wolfram))`, maximised over the varying z vectors.

| dps | Case | Characteristic | Maximum scaled difference | Diagonal factorisation error |
|---:|---|---|---:|---:|
| 50 | g1 | zero | `1.044049e-52` | — |
| 50 | g1 | half | `6.887557e-52` | — |
| 50 | g2-diagonal | zero | `1.346782e-51` | `2.42447e-51` |
| 50 | g2-diagonal | half | `1.336382e-51` | `1.377511e-51` |
| 50 | g2-coupled | zero | `3.09377e-52` | — |
| 50 | g2-coupled | half | `8.35239e-53` | — |
| 50 | g2-coupled | dz0 | `1.336382e-51` | — |
| 50 | g2-coupled | half-dz0 | `1.55438e-51` | — |
| 50 | g3-diagonal | zero | `2.39602e-51` | `2.584997e-51` |
| 50 | g3-diagonal | half | `3.340956e-52` | `3.340956e-52` |
| 50 | g3-coupled | zero | `2.172027e-52` | — |
| 50 | g3-coupled | half | `3.340956e-52` | — |
| 50 | g3-coupled | dz0 | `5.011434e-52` | — |
| 50 | g3-coupled | half-dz0 | `1.336382e-51` | — |
| 50 | g2-near-boundary | zero | `2.532076e-51` | — |
| 50 | g2-near-boundary | half | `1.377511e-51` | — |
| 50 | g2-unreduced | zero | `6.343541e-50` | — |
| 50 | g2-unreduced | half | `5.542208e-50` | — |
| 50 | g2-unreduced | real | `4.940529e-50` | — |
| 50 | g3-unreduced | zero | `1.909158e-51` | — |
| 50 | g3-unreduced | half | `1.670478e-51` | — |
| 50 | g3-unreduced | real | `4.226012e-51` | — |
| 50 | g3-block | zero | `3.246185e-52` | `2.125883e-51` |
| 50 | g3-block | half | `4.724825e-52` | `3.347475e-52` |
| 100 | g1 | zero | `1.418256e-101` | — |
| 100 | g1 | half | `7.143671e-102` | — |
| 100 | g2-diagonal | zero | `1.277227e-102` | `1.363041e-101` |
| 100 | g2-diagonal | half | `7.143671e-102` | `7.363528e-102` |
| 100 | g2-coupled | zero | `8.929589e-103` | — |
| 100 | g2-coupled | half | `7.363528e-102` | — |
| 100 | g2-coupled | dz0 | `7.363528e-102` | — |
| 100 | g2-coupled | half-dz0 | `7.7924e-102` | — |
| 100 | g3-diagonal | zero | `1.383838e-101` | `1.381818e-101` |
| 100 | g3-diagonal | half | `1.785918e-102` | `1.799816e-102` |
| 100 | g3-coupled | zero | `1.342184e-101` | — |
| 100 | g3-coupled | half | `2.525669e-102` | — |
| 100 | g2-near-boundary | zero | `1.345697e-101` | — |
| 100 | g2-near-boundary | half | `7.986867e-102` | — |
| 100 | g2-unreduced | zero | `2.663222e-100` | — |
| 100 | g2-unreduced | half | `2.172614e-100` | — |
| 100 | g2-unreduced | real | `2.222519e-100` | — |
| 100 | g3-unreduced | zero | `1.000126e-101` | — |
| 100 | g3-unreduced | half | `7.363528e-102` | — |
| 100 | g3-unreduced | real | `3.29307e-101` | — |
| 100 | g3-block | zero | `8.676288e-103` | `1.135843e-101` |
| 100 | g3-block | half | `1.785918e-102` | `1.785918e-102` |

## Representative coupled-case values

Values below are included to make convention or parsing errors visible rather than hiding them behind aggregate error figures.

| dps | Case | Characteristic | mpmath first value | Wolfram first value |
|---:|---|---|---:|---:|
| 50 | g2-coupled | zero | `1.174894515 + 0.05306062102j` | `1.174894515 + 0.05306062102j` |
| 50 | g2-coupled | half | `0.9370524852 + 0.05319331675j` | `0.9370524852 + 0.05319331675j` |
| 50 | g2-coupled | dz0 | `-0.4212883156 - 0.2974781291j` | `-0.4212883156 - 0.2974781291j` |
| 50 | g2-coupled | half-dz0 | `-0.9837401351 - 0.2877326014j` | `-0.9837401351 - 0.2877326014j` |
| 50 | g3-coupled | zero | `1.248338655 + 0.006475364864j` | `1.248338655 + 0.006475364864j` |
| 50 | g3-coupled | half | `1.380160949e-53 + 5.46626337e-53j` | `0.0 + 0.0j` |
| 50 | g3-coupled | dz0 | `-0.3493647391 - 0.2184857273j` | `-0.3493647391 - 0.2184857273j` |
| 50 | g3-coupled | half-dz0 | `0.3247585964 - 0.2202690092j` | `0.3247585964 - 0.2202690092j` |
| 100 | g2-coupled | zero | `1.174894515 + 0.05306062102j` | `1.174894515 + 0.05306062102j` |
| 100 | g2-coupled | half | `0.9370524852 + 0.05319331675j` | `0.9370524852 + 0.05319331675j` |
| 100 | g2-coupled | dz0 | `-0.4212883156 - 0.2974781291j` | `-0.4212883156 - 0.2974781291j` |
| 100 | g2-coupled | half-dz0 | `-0.9837401351 - 0.2877326014j` | `-0.9837401351 - 0.2877326014j` |
| 100 | g3-coupled | zero | `1.248338655 + 0.006475364864j` | `1.248338655 + 0.006475364864j` |
| 100 | g3-coupled | half | `-7.236026087e-104 + 4.473478316e-103j` | `0.0 + 0.0j` |

## Cases

| Case | Genus | Type |
|---|---:|---|
| g1 | 1 | genus-one control |
| g2-diagonal | 2 | factorisation control |
| g2-coupled | 2 | coupled matrix |
| g3-diagonal | 3 | factorisation control |
| g3-coupled | 3 | coupled matrix |
| g2-near-boundary | 2 | near-boundary direct-path control |
| g2-unreduced | 2 | strongly unreduced matrix |
| g3-unreduced | 3 | strongly unreduced matrix |
| g3-block | 3 | genus-one plus genus-two block control |

The diagonal cases independently verify that a genus-g theta value equals the product of its genus-one components, for both tested characteristics. The coupled cases have nonzero off-diagonal real and imaginary parts and therefore do not factor.

Wolfram calls use `SiegelTheta[tau, z]` and `SiegelTheta[{a, b}, tau, z]`; no `EllipticTheta` calls are made.

## Incomplete workloads

A failure here does not discard successful rows from the same run.

| dps | Case | Workload | Reason |
|---:|---|---|---|
| 100 | g3-coupled | dz0 | skipped: Wolfram contour differentiation is impractical above 50 dps |
| 100 | g3-coupled | half-dz0 | skipped: Wolfram contour differentiation is impractical above 50 dps |
