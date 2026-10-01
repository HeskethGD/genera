# Hyperelliptic period comparison with Sage

This compares mpmath with Sage's arbitrary-precision Riemann-surface period calculation and abelfunctions' machine-precision calculation. The former is the period engine used by the nbruin/RiemannTheta README; RiemannTheta itself begins with an already computed period matrix.

Sage and abelfunctions represent the standard algebraic differential as `x^k dx/(df/dy) = x^k dx/(2y)`. Their matrices are therefore multiplied by the known factor two before recovering the nearest integral change of cycle basis. No numerical scale is fitted.

## Reproduction

```console
.venv/bin/python my_data/riemann_theta/benchmark_hyperelliptic_vs_sage.py --cases g1-symmetric g2-irregular
```

- Date: 2026-10-01T14:05:56.722801+00:00
- Platform: macOS-27.0.1-arm64-arm-64bit-Mach-O
- Python: 3.13.4
- mpmath: 0.1.dev5+g5f0e3bed1.d20261001
- Sage: `/usr/local/bin/sage`
- mpmath timing: median of 3 trials
- Sage timing: median of 1 trials inside one Sage process; interpreter startup and the initial warmup are excluded

## Accuracy

| Implementation | Case | Genus | dps | Relative lattice residual | Integer-map error | Orientation | Form error |
|---|---|---:|---:|---:|---:|---|---:|
| sage | g1-symmetric | 1 | 30 | 3.47e-30 | 6.165e-30 | same | 0.0 |
| sage | g2-irregular | 2 | 30 | 2.845e-29 | 5.425e-29 | same | 0.0 |
| sage | g1-symmetric | 1 | 50 | 1.447e-49 | 2.5e-49 | same | 0.0 |
| sage | g2-irregular | 2 | 50 | 2.385e-49 | 5.071e-49 | same | 0.0 |
| abelfunctions | g1-symmetric | 1 | 15 | 1.383e-16 | 3.14e-16 | same | 0.0 |
| abelfunctions | g2-irregular | 2 | 15 | 2.941e-16 | 7.881e-16 | same | 0.0 |

abelfunctions stores periods in NumPy complex128 arrays, so its dps column denotes the approximately 15-digit comparison target rather than user-selectable arbitrary precision.

## Performance

The timed workload is the complete curve-to-period calculation, including topology/root work and integration but excluding imports and construction of the polynomial input.

| Implementation | Case | Genus | dps | mpmath (ms) | Other (ms) | mpmath / other |
|---|---|---:|---:|---:|---:|---:|
| sage | g1-symmetric | 1 | 30 | 10.801 | 124.331 | 0.0869x |
| sage | g2-irregular | 2 | 30 | 53.135 | 581.420 | 0.0914x |
| sage | g1-symmetric | 1 | 50 | 24.933 | 202.789 | 0.123x |
| sage | g2-irregular | 2 | 50 | 99.838 | 842.669 | 0.118x |
| abelfunctions | g1-symmetric | 1 | 15 | 5.104 | 77.508 | 0.0658x |
| abelfunctions | g2-irregular | 2 | 15 | 25.447 | 1197.332 | 0.0213x |

### Geometric means

| Implementation | dps | mpmath (ms) | Other (ms) | mpmath / other |
|---|---:|---:|---:|---:|
| abelfunctions | 15 | 11.396 | 304.636 | 0.0374x |
| sage | 30 | 23.956 | 268.865 | 0.0891x |
| sage | 50 | 49.893 | 413.381 | 0.121x |
