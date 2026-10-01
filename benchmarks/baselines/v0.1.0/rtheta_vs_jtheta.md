# Riemann theta genus-one baseline: `rtheta` vs `jtheta`

This compares the generic Riemann theta evaluator, including its cost-selected Siegel reduction, with mpmath's specialised genus-one `jtheta` implementation. It is not expected to outperform `jtheta` on ordinary genus-one inputs.

Each timing evaluates 6 deterministic complex arguments. Times are median wall-clock microseconds per individual function call; lower is better.

- Generated: `2026-10-01T15:11:53.402153+01:00`
- mpmath revision: `c311b91 (dirty)`
- Python: `3.13.4`
- Platform: `macOS-27.0.1-arm64-arm-64bit-Mach-O`
- Decimal precisions: `50, 100`
- Timing repeats: `5` (median reported)
- Batches per repeat: `1`
- Arguments per batch: `6`
- Fixed-tau sweep arguments per batch: `50`
- Fixed-tau sweep repeats: `3` (median reported)

## Summary

- At 50 dps, the geometric-mean `rtheta / jtheta` ratio was **1.32x warm** and **8.59x cold**. The largest scaled error was `7.6722e-49`.
- At 100 dps, the geometric-mean `rtheta / jtheta` ratio was **1.23x warm** and **7.66x cold**. The largest scaled error was `1.4922e-98`.

## Timing results

| dps | Case | tau | Workload | rtheta warm (us) | rtheta cold (us) | jtheta (us) | Warm / jtheta | Cold / jtheta |
|---:|---|---:|---|---:|---:|---:|---:|---:|
| 50 | comfortable | `(0.2 + 1.2j)` | theta1 | 561.03 | 2463.58 | 136.66 | 4.11x | 18.03x |
| 50 | comfortable | `(0.2 + 1.2j)` | theta2 | 302.31 | 1769.88 | 131.76 | 2.29x | 13.43x |
| 50 | comfortable | `(0.2 + 1.2j)` | theta3 | 327.56 | 1761.18 | 111.01 | 2.95x | 15.86x |
| 50 | comfortable | `(0.2 + 1.2j)` | theta4 | 332.15 | 1816.51 | 110.71 | 3.00x | 16.41x |
| 50 | comfortable | `(0.2 + 1.2j)` | theta3 d1 | 336.95 | 3455.89 | 102.63 | 3.28x | 33.67x |
| 50 | comfortable | `(0.2 + 1.2j)` | theta3 d2 | 343.84 | 3937.84 | 110.28 | 3.12x | 35.71x |
| 50 | generic | `(0.35 + 0.75j)` | theta1 | 300.92 | 1624.19 | 316.38 | 0.95x | 5.13x |
| 50 | generic | `(0.35 + 0.75j)` | theta2 | 314.17 | 1672.16 | 333.75 | 0.94x | 5.01x |
| 50 | generic | `(0.35 + 0.75j)` | theta3 | 284.74 | 1626.29 | 308.09 | 0.92x | 5.28x |
| 50 | generic | `(0.35 + 0.75j)` | theta4 | 289.74 | 1618.17 | 316.16 | 0.92x | 5.12x |
| 50 | generic | `(0.35 + 0.75j)` | theta3 d1 | 364.17 | 3490.21 | 435.91 | 0.84x | 8.01x |
| 50 | generic | `(0.35 + 0.75j)` | theta3 d2 | 724.99 | 3983.56 | 583.15 | 1.24x | 6.83x |
| 50 | near-boundary | `(0.2 + 0.3j)` | theta1 | 338.88 | 1664.44 | 329.36 | 1.03x | 5.05x |
| 50 | near-boundary | `(0.2 + 0.3j)` | theta2 | 343.87 | 1677.07 | 292.85 | 1.17x | 5.73x |
| 50 | near-boundary | `(0.2 + 0.3j)` | theta3 | 333.12 | 1677.48 | 287.42 | 1.16x | 5.84x |
| 50 | near-boundary | `(0.2 + 0.3j)` | theta4 | 334.19 | 1656.60 | 320.84 | 1.04x | 5.16x |
| 50 | near-boundary | `(0.2 + 0.3j)` | theta3 d1 | 441.76 | 3628.39 | 409.55 | 1.08x | 8.86x |
| 50 | near-boundary | `(0.2 + 0.3j)` | theta3 d2 | 470.97 | 5163.87 | 705.53 | 0.67x | 7.32x |
| 50 | strongly-unreduced | `(0.17 + 0.003j)` | theta1 | 498.51 | 3274.79 | 456.17 | 1.09x | 7.18x |
| 50 | strongly-unreduced | `(0.17 + 0.003j)` | theta2 | 510.08 | 3274.32 | 444.50 | 1.15x | 7.37x |
| 50 | strongly-unreduced | `(0.17 + 0.003j)` | theta3 | 483.72 | 3248.77 | 417.10 | 1.16x | 7.79x |
| 50 | strongly-unreduced | `(0.17 + 0.003j)` | theta4 | 477.24 | 3323.26 | 422.69 | 1.13x | 7.86x |
| 50 | strongly-unreduced | `(0.17 + 0.003j)` | theta3 d1 | 673.24 | 5792.53 | 723.56 | 0.93x | 8.01x |
| 50 | strongly-unreduced | `(0.17 + 0.003j)` | theta3 d2 | 751.42 | 6903.33 | 1117.58 | 0.67x | 6.18x |
| 100 | comfortable | `(0.2 + 1.2j)` | theta1 | 358.27 | 1882.47 | 185.10 | 1.94x | 10.17x |
| 100 | comfortable | `(0.2 + 1.2j)` | theta2 | 348.20 | 1896.37 | 174.01 | 2.00x | 10.90x |
| 100 | comfortable | `(0.2 + 1.2j)` | theta3 | 364.86 | 1908.62 | 144.19 | 2.53x | 13.24x |
| 100 | comfortable | `(0.2 + 1.2j)` | theta4 | 349.67 | 1905.83 | 145.33 | 2.41x | 13.11x |
| 100 | comfortable | `(0.2 + 1.2j)` | theta3 d1 | 424.55 | 4027.99 | 144.67 | 2.93x | 27.84x |
| 100 | comfortable | `(0.2 + 1.2j)` | theta3 d2 | 450.35 | 4631.85 | 163.20 | 2.76x | 28.38x |
| 100 | generic | `(0.35 + 0.75j)` | theta1 | 390.94 | 1938.68 | 408.65 | 0.96x | 4.74x |
| 100 | generic | `(0.35 + 0.75j)` | theta2 | 389.33 | 1929.83 | 383.01 | 1.02x | 5.04x |
| 100 | generic | `(0.35 + 0.75j)` | theta3 | 380.52 | 1940.77 | 382.56 | 0.99x | 5.07x |
| 100 | generic | `(0.35 + 0.75j)` | theta4 | 376.24 | 1926.31 | 402.92 | 0.93x | 4.78x |
| 100 | generic | `(0.35 + 0.75j)` | theta3 d1 | 469.44 | 4057.76 | 601.88 | 0.78x | 6.74x |
| 100 | generic | `(0.35 + 0.75j)` | theta3 d2 | 506.80 | 4665.03 | 786.51 | 0.64x | 5.93x |
| 100 | near-boundary | `(0.2 + 0.3j)` | theta1 | 466.88 | 2094.83 | 404.58 | 1.15x | 5.18x |
| 100 | near-boundary | `(0.2 + 0.3j)` | theta2 | 464.50 | 2087.36 | 351.32 | 1.32x | 5.94x |
| 100 | near-boundary | `(0.2 + 0.3j)` | theta3 | 459.69 | 2093.89 | 340.94 | 1.35x | 6.14x |
| 100 | near-boundary | `(0.2 + 0.3j)` | theta4 | 463.76 | 2089.59 | 396.69 | 1.17x | 5.27x |
| 100 | near-boundary | `(0.2 + 0.3j)` | theta3 d1 | 593.42 | 4269.63 | 531.78 | 1.12x | 8.03x |
| 100 | near-boundary | `(0.2 + 0.3j)` | theta3 d2 | 692.77 | 4897.44 | 753.40 | 0.92x | 6.50x |
| 100 | strongly-unreduced | `(0.17 + 0.003j)` | theta1 | 605.99 | 3987.56 | 579.74 | 1.05x | 6.88x |
| 100 | strongly-unreduced | `(0.17 + 0.003j)` | theta2 | 620.94 | 3926.09 | 586.31 | 1.06x | 6.70x |
| 100 | strongly-unreduced | `(0.17 + 0.003j)` | theta3 | 587.83 | 3864.97 | 550.18 | 1.07x | 7.02x |
| 100 | strongly-unreduced | `(0.17 + 0.003j)` | theta4 | 578.90 | 3876.44 | 553.15 | 1.05x | 7.01x |
| 100 | strongly-unreduced | `(0.17 + 0.003j)` | theta3 d1 | 788.51 | 6569.80 | 895.82 | 0.88x | 7.33x |
| 100 | strongly-unreduced | `(0.17 + 0.003j)` | theta3 d2 | 903.72 | 7953.87 | 1438.66 | 0.63x | 5.53x |
| 50 | **geometric mean** | — | all workloads | 403.79 | 2636.21 | 306.88 | **1.32x** | **8.59x** |
| 100 | **geometric mean** | — | all workloads | 484.04 | 3013.49 | 393.23 | **1.23x** | **7.66x** |

The geometric-mean rows aggregate every case and workload at each precision.

A ratio above 1 means `jtheta` was faster. The cold measurement includes clearing and rebuilding `rtheta`'s precision-aware tau-preprocessing and derivative-radius caches before every call. The warm measurement uses one fixed tau while varying z.

## Numerical agreement

The reported error is `abs(rtheta - converted_jtheta) / max(1, abs(converted_jtheta))`. It is measured before timing.

| dps | Case | Workload | Maximum scaled error |
|---:|---|---|---:|
| 50 | comfortable | theta1 | `2.6888611e-52` |
| 50 | comfortable | theta2 | `8.8831798e-52` |
| 50 | comfortable | theta3 | `1.3391273e-51` |
| 50 | comfortable | theta4 | `7.9998092e-52` |
| 50 | comfortable | theta3 d1 | `4.8721672e-52` |
| 50 | comfortable | theta3 d2 | `2.4192909e-51` |
| 50 | generic | theta1 | `5.4935553e-52` |
| 50 | generic | theta2 | `1.5277683e-51` |
| 50 | generic | theta3 | `1.1616607e-51` |
| 50 | generic | theta4 | `8.1500524e-52` |
| 50 | generic | theta3 d1 | `1.3483178e-51` |
| 50 | generic | theta3 d2 | `1.2353363e-51` |
| 50 | near-boundary | theta1 | `5.6437988e-52` |
| 50 | near-boundary | theta2 | `9.0654085e-52` |
| 50 | near-boundary | theta3 | `8.8926201e-52` |
| 50 | near-boundary | theta4 | `1.0066624e-51` |
| 50 | near-boundary | theta3 d1 | `1.3119826e-51` |
| 50 | near-boundary | theta3 d2 | `1.5585209e-51` |
| 50 | strongly-unreduced | theta1 | `7.9500477e-50` |
| 50 | strongly-unreduced | theta2 | `7.8879885e-50` |
| 50 | strongly-unreduced | theta3 | `7.8120218e-50` |
| 50 | strongly-unreduced | theta4 | `9.662492e-50` |
| 50 | strongly-unreduced | theta3 d1 | `7.672248e-49` |
| 50 | strongly-unreduced | theta3 d2 | `1.2872197e-49` |
| 100 | comfortable | theta1 | `3.5925656e-102` |
| 100 | comfortable | theta2 | `9.6955659e-102` |
| 100 | comfortable | theta3 | `9.2021281e-102` |
| 100 | comfortable | theta4 | `8.3632109e-102` |
| 100 | comfortable | theta3 d1 | `4.7860861e-102` |
| 100 | comfortable | theta3 d2 | `2.4038232e-101` |
| 100 | generic | theta1 | `5.2169226e-102` |
| 100 | generic | theta2 | `5.5534255e-102` |
| 100 | generic | theta3 | `7.0694102e-102` |
| 100 | generic | theta4 | `3.3828031e-102` |
| 100 | generic | theta3 d1 | `8.8632263e-102` |
| 100 | generic | theta3 d2 | `8.9528365e-102` |
| 100 | near-boundary | theta1 | `6.0128109e-102` |
| 100 | near-boundary | theta2 | `1.0742377e-101` |
| 100 | near-boundary | theta3 | `1.0632593e-101` |
| 100 | near-boundary | theta4 | `5.3618238e-102` |
| 100 | near-boundary | theta3 d1 | `7.8539288e-102` |
| 100 | near-boundary | theta3 d2 | `1.2008934e-101` |
| 100 | strongly-unreduced | theta1 | `1.3152588e-99` |
| 100 | strongly-unreduced | theta2 | `1.3187435e-99` |
| 100 | strongly-unreduced | theta3 | `1.0926056e-99` |
| 100 | strongly-unreduced | theta4 | `1.5258939e-99` |
| 100 | strongly-unreduced | theta3 d1 | `1.4922245e-98` |
| 100 | strongly-unreduced | theta3 d2 | `2.5033463e-99` |

## Cases and conventions

| Case | tau | Purpose |
|---|---:|---|
| comfortable | `0.20 + 1.20 i` | rapid direct convergence |
| generic | `0.35 + 0.75 i` | moderate direct convergence |
| near-boundary | `0.20 + 0.30 i` | transformation overhead outweighs the small predicted saving |
| strongly-unreduced | `0.17 + 0.003 i` | benefits substantially from modular reduction |

The six angular `jtheta` arguments are:

```text
0.00+0.00i, 0.11+0.03i, -0.23+0.07i, 0.37-0.09i, -0.41-0.12i, 0.52+0.16i
```

They are converted using `z = w/pi` and `q = exp(pi*i*tau)`. The four standard half-characteristics use the agreed sign for theta1. A derivative of order d is compared using `d^d rtheta/dz^d = pi^d d^d jtheta/dw^d`.

The near-boundary case checks that the cost gate retains direct summation when a modular transform would save too little work. The strongly-unreduced case checks the selected transformation path. `jtheta` applies its own specialised genus-one transformations.

## Fixed-tau, varying-argument sweeps

This plot-like workload evaluates a deterministic sequence of 50 arguments while keeping tau and the characteristic fixed. Both implementations receive preconstructed arguments, and `rtheta`'s precision-aware tau and derivative-radius caches are populated before timing.

- At 50 dps, the geometric-mean `rtheta / jtheta` ratio was **1.14x for values** and **1.00x for derivatives**.
- At 100 dps, the geometric-mean `rtheta / jtheta` ratio was **1.16x for values** and **0.99x for derivatives**.

| dps | Case | Workload | rtheta (us) | jtheta (us) | rtheta / jtheta | Maximum scaled error |
|---:|---|---|---:|---:|---:|---:|
| 50 | comfortable | theta1 | 262.23 | 130.09 | 2.02x | `3.5415831e-51` |
| 50 | comfortable | theta2 | 264.74 | 125.05 | 2.12x | `3.9828184e-51` |
| 50 | comfortable | theta3 | 271.60 | 109.72 | 2.48x | `1.5826567e-51` |
| 50 | comfortable | theta4 | 270.20 | 109.12 | 2.48x | `1.3854487e-51` |
| 50 | comfortable | theta3 d1 | 330.12 | 113.18 | 2.92x | `1.5449265e-51` |
| 50 | comfortable | theta3 d2 | 344.04 | 121.04 | 2.84x | `2.3982573e-50` |
| 50 | generic | theta1 | 286.95 | 373.96 | 0.77x | `4.8344746e-51` |
| 50 | generic | theta2 | 285.84 | 354.14 | 0.81x | `5.4700337e-51` |
| 50 | generic | theta3 | 284.37 | 356.67 | 0.80x | `2.9888985e-51` |
| 50 | generic | theta4 | 282.02 | 370.40 | 0.76x | `2.9443469e-51` |
| 50 | generic | theta3 d1 | 371.74 | 530.67 | 0.70x | `4.1911333e-51` |
| 50 | generic | theta3 d2 | 378.43 | 715.73 | 0.53x | `8.9579571e-50` |
| 50 | near-boundary | theta1 | 328.68 | 383.08 | 0.86x | `7.2380245e-51` |
| 50 | near-boundary | theta2 | 331.27 | 329.34 | 1.01x | `9.943735e-51` |
| 50 | near-boundary | theta3 | 333.91 | 336.23 | 0.99x | `1.0027234e-50` |
| 50 | near-boundary | theta4 | 328.03 | 377.47 | 0.87x | `9.324295e-51` |
| 50 | near-boundary | theta3 d1 | 439.19 | 482.28 | 0.91x | `6.544671e-51` |
| 50 | near-boundary | theta3 d2 | 465.74 | 643.18 | 0.72x | `4.3508629e-50` |
| 50 | strongly-unreduced | theta1 | 519.40 | 505.50 | 1.03x | `1.5459831e-49` |
| 50 | strongly-unreduced | theta2 | 507.62 | 498.66 | 1.02x | `1.5361202e-49` |
| 50 | strongly-unreduced | theta3 | 500.03 | 471.45 | 1.06x | `3.1418007e-49` |
| 50 | strongly-unreduced | theta4 | 481.30 | 448.69 | 1.07x | `2.559615e-49` |
| 50 | strongly-unreduced | theta3 d1 | 679.80 | 815.45 | 0.83x | `3.582689e-49` |
| 50 | strongly-unreduced | theta3 d2 | 748.64 | 1282.90 | 0.58x | `4.6942711e-49` |
| 100 | comfortable | theta1 | 342.02 | 193.65 | 1.77x | `2.3659336e-101` |
| 100 | comfortable | theta2 | 344.76 | 187.43 | 1.84x | `2.3906155e-101` |
| 100 | comfortable | theta3 | 350.06 | 160.49 | 2.18x | `1.1060932e-101` |
| 100 | comfortable | theta4 | 352.20 | 164.82 | 2.14x | `1.3214063e-101` |
| 100 | comfortable | theta3 d1 | 421.19 | 161.41 | 2.61x | `8.5147281e-102` |
| 100 | comfortable | theta3 d2 | 466.31 | 177.13 | 2.63x | `1.4585897e-100` |
| 100 | generic | theta1 | 380.45 | 467.20 | 0.81x | `3.1155517e-101` |
| 100 | generic | theta2 | 389.76 | 449.03 | 0.87x | `3.4403831e-101` |
| 100 | generic | theta3 | 368.97 | 447.82 | 0.82x | `1.9777893e-101` |
| 100 | generic | theta4 | 373.30 | 462.76 | 0.81x | `1.9085689e-101` |
| 100 | generic | theta3 d1 | 467.22 | 714.76 | 0.65x | `3.044883e-101` |
| 100 | generic | theta3 d2 | 519.77 | 942.18 | 0.55x | `5.795141e-100` |
| 100 | near-boundary | theta1 | 457.14 | 470.03 | 0.97x | `5.3487671e-101` |
| 100 | near-boundary | theta2 | 454.36 | 395.04 | 1.15x | `5.5013408e-101` |
| 100 | near-boundary | theta3 | 467.23 | 398.39 | 1.17x | `6.4424024e-101` |
| 100 | near-boundary | theta4 | 464.50 | 445.45 | 1.04x | `6.4318387e-101` |
| 100 | near-boundary | theta3 d1 | 584.99 | 630.11 | 0.93x | `3.5728164e-101` |
| 100 | near-boundary | theta3 d2 | 694.61 | 841.23 | 0.83x | `2.577385e-100` |
| 100 | strongly-unreduced | theta1 | 615.85 | 626.54 | 0.98x | `2.0336791e-99` |
| 100 | strongly-unreduced | theta2 | 616.52 | 609.28 | 1.01x | `1.3475134e-99` |
| 100 | strongly-unreduced | theta3 | 587.04 | 566.57 | 1.04x | `2.548374e-99` |
| 100 | strongly-unreduced | theta4 | 582.01 | 554.97 | 1.05x | `1.7537724e-99` |
| 100 | strongly-unreduced | theta3 d1 | 816.38 | 1003.37 | 0.81x | `2.4496621e-99` |
| 100 | strongly-unreduced | theta3 d2 | 941.24 | 1601.16 | 0.59x | `4.0322772e-99` |
| 50 | **geometric mean** | values | 335.17 | 292.74 | **1.14x** | `3.1418007e-49` |
| 50 | **geometric mean** | derivatives | 449.39 | 450.23 | **1.00x** | `4.6942711e-49` |
| 100 | **geometric mean** | values | 436.49 | 376.80 | **1.16x** | `2.548374e-99` |
| 100 | **geometric mean** | derivatives | 591.04 | 598.19 | **0.99x** | `4.0322772e-99` |

Sweep geometric means are separated into values and derivatives at each precision.

The sweep arguments follow the same deterministic distribution as the existing `benchmark_jtheta_boundaries_pyperf.py` benchmark: their real parts cover 0.1 to 0.9 and their small imaginary parts are distributed using the golden-ratio conjugate. This represents repeated evaluation along a plot or trajectory without adding vector-valued semantics to either function.
