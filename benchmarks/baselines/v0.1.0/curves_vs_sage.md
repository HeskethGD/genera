# Curve-to-full-periods benchmark

Timing excludes imports, fixed warmup, serialization and validation. Fresh process per sample.
Sage discovers its own basis. Bases are aligned exactly; no fitted factor-of-two correction.
Target digits label the accuracy request. Both engines receive 15 extra working digits, including coefficient conversion; their cost is timed.

| Case | Digits | Implementation | Successes/trials | Median seconds | Accuracy |
|---|---:|---|---:|---:|---|
| hyper-g1-lemniscatic | 30 | hyperelliptic | 3/3 | 0.0179 | passed |
| hyper-g1-lemniscatic | 30 | sage | 3/3 | 0.1682 | passed |
| hyper-g1-lemniscatic | 50 | hyperelliptic | 3/3 | 0.0404 | passed |
| hyper-g1-lemniscatic | 50 | sage | 3/3 | 0.2677 | passed |
| hyper-g2-symmetric | 30 | hyperelliptic | 3/3 | 0.0661 | passed |
| hyper-g2-symmetric | 30 | sage | 3/3 | 0.5083 | passed |
| hyper-g2-symmetric | 50 | hyperelliptic | 3/3 | 0.1488 | passed |
| hyper-g2-symmetric | 50 | sage | 3/3 | 0.7678 | passed |
| hyper-g3-symmetric | 30 | hyperelliptic | 3/3 | 0.1617 | passed |
| hyper-g3-symmetric | 30 | sage | 3/3 | 1.1860 | passed |
| hyper-g3-symmetric | 50 | hyperelliptic | 3/3 | 0.3637 | passed |
| hyper-g3-symmetric | 50 | sage | 3/3 | 1.6402 | passed |
| hyper-g2-even | 30 | hyperelliptic | 3/3 | 0.0875 | passed |
| hyper-g2-even | 30 | sage | 3/3 | 0.7713 | passed |
| hyper-g2-even | 50 | hyperelliptic | 3/3 | 0.1969 | passed |
| hyper-g2-even | 50 | sage | 3/3 | 1.1563 | passed |
| hyper-g2-clustered | 30 | hyperelliptic | 3/3 | 0.0995 | passed |
| hyper-g2-clustered | 30 | sage | 3/3 | 0.9525 | passed |
| hyper-g2-clustered | 50 | hyperelliptic | 3/3 | 0.1350 | passed |
| hyper-g2-clustered | 50 | sage | 3/3 | 1.3821 | passed |
| hyper-g3-generic | 30 | hyperelliptic | 3/3 | 0.1831 | passed |
| hyper-g3-generic | 30 | sage | 3/3 | 15.6123 | passed |
| hyper-g3-generic | 50 | hyperelliptic | 3/3 | 0.4825 | passed |
| hyper-g3-generic | 50 | sage | 3/3 | 18.3732 | passed |
| hyper-g4-symmetric | 30 | hyperelliptic | 3/3 | 0.3139 | passed |
| hyper-g4-symmetric | 30 | sage | 3/3 | 3.8939 | passed |
| hyper-g4-symmetric | 50 | hyperelliptic | 3/3 | 0.7150 | passed |
| hyper-g4-symmetric | 50 | sage | 3/3 | 4.6478 | passed |
| general-g1-lemniscatic | 30 | genera | 3/3 | 0.7022 | passed |
| general-g1-lemniscatic | 30 | sage | 3/3 | 0.1678 | passed |
| general-g1-lemniscatic | 50 | genera | 3/3 | 1.2219 | passed |
| general-g1-lemniscatic | 50 | sage | 3/3 | 0.2701 | passed |
| general-trig-g3 | 30 | genera | 3/3 | 2.4747 | passed |
| general-trig-g3 | 30 | sage | 3/3 | 2.5890 | passed |
| general-trig-g3 | 50 | genera | 3/3 | 3.8921 | passed |
| general-trig-g3 | 50 | sage | 3/3 | 3.6729 | passed |
| general-trig-pure-g3 | 30 | genera | 3/3 | 1.2305 | passed |
| general-trig-pure-g3 | 30 | sage | 3/3 | 0.8303 | passed |
| general-trig-pure-g3 | 50 | genera | 3/3 | 2.0086 | passed |
| general-trig-pure-g3 | 50 | sage | 3/3 | 1.2284 | passed |
| general-trig-dense-g3 | 30 | genera | 3/3 | 2.2787 | passed |
| general-trig-dense-g3 | 30 | sage | 3/3 | 3.3166 | passed |
| general-trig-dense-g3 | 50 | genera | 3/3 | 3.3456 | passed |
| general-trig-dense-g3 | 50 | sage | 3/3 | 4.0840 | passed |
| general-klein-quartic-g3 | 30 | genera | 3/3 | 1.8868 | passed |
| general-klein-quartic-g3 | 30 | sage | 3/3 | 0.9493 | passed |
| general-klein-quartic-g3 | 50 | genera | 3/3 | 2.9192 | passed |
| general-klein-quartic-g3 | 50 | sage | 3/3 | 1.4473 | passed |
| general-fermat-quartic-g3 | 30 | genera | 3/3 | 1.6067 | passed |
| general-fermat-quartic-g3 | 30 | sage | 3/3 | 1.1873 | passed |
| general-fermat-quartic-g3 | 50 | genera | 3/3 | 2.6083 | passed |
| general-fermat-quartic-g3 | 50 | sage | 3/3 | 2.0824 | passed |
| general-super-g4 | 30 | genera | 3/3 | 1.5576 | passed |
| general-super-g4 | 30 | sage | 3/3 | 0.9087 | passed |
| general-super-g4 | 50 | genera | 3/3 | 2.5104 | passed |
| general-super-g4 | 50 | sage | 3/3 | 1.3476 | passed |
| general-kovalevskaya-g3 | 30 | genera | 3/3 | 11.5681 | passed |
| general-kovalevskaya-g3 | 30 | sage | 3/3 | 3.4526 | passed |
| general-kovalevskaya-g3 | 50 | genera | 3/3 | 17.3893 | passed |
| general-kovalevskaya-g3 | 50 | sage | 3/3 | 5.1382 | passed |
| nonmonic-elliptic | 30 | genera | 3/3 | 0.4870 | passed |
| nonmonic-elliptic | 30 | sage | 3/3 | 0.1797 | passed |
| nonmonic-elliptic | 50 | genera | 3/3 | 0.7903 | passed |
| nonmonic-elliptic | 50 | sage | 3/3 | 0.2891 | passed |
| general-trig-test-g3 | 30 | genera | 3/3 | 1.2913 | passed |
| general-trig-test-g3 | 30 | sage | 3/3 | 0.5731 | passed |
| general-trig-test-g3 | 50 | genera | 3/3 | 1.9211 | passed |
| general-trig-test-g3 | 50 | sage | 3/3 | 0.8356 | passed |
| general-trig-infinity-g3 | 30 | genera | 3/3 | 1.2954 | passed |
| general-trig-infinity-g3 | 30 | sage | 3/3 | 1.1318 | passed |
| general-trig-infinity-g3 | 50 | genera | 3/3 | 2.0972 | passed |
| general-trig-infinity-g3 | 50 | sage | 3/3 | 2.0045 | passed |
| klein-x-translated | 30 | genera | 3/3 | 2.7261 | passed |
| klein-x-translated | 30 | sage | 3/3 | 2.0931 | passed |
| klein-x-translated | 50 | genera | 3/3 | 3.9036 | passed |
| klein-x-translated | 50 | sage | 3/3 | 2.2719 | passed |
| klein-x-scaled | 30 | genera | 3/3 | 1.8711 | passed |
| klein-x-scaled | 30 | sage | 3/3 | 0.9146 | passed |
| klein-x-scaled | 50 | genera | 3/3 | 2.9047 | passed |
| klein-x-scaled | 50 | sage | 3/3 | 1.4467 | passed |
| general-mu-x-translated | 30 | genera | 3/3 | 3.1592 | passed |
| general-mu-x-translated | 30 | sage | 3/3 | 6.5084 | passed |
| general-mu-x-translated | 50 | genera | 3/3 | 4.7009 | passed |
| general-mu-x-translated | 50 | sage | 3/3 | 6.4515 | passed |
| hyper-g2-close-001 | 30 | hyperelliptic | 3/3 | 0.1117 | passed |
| hyper-g2-close-001 | 30 | sage | 3/3 | 0.7874 | passed |
| hyper-g2-close-001 | 50 | hyperelliptic | 3/3 | 0.2097 | passed |
| hyper-g2-close-001 | 50 | sage | 3/3 | 1.2581 | passed |
| hyper-g2-close-00001 | 30 | hyperelliptic | 3/3 | 0.1584 | passed |
| hyper-g2-close-00001 | 30 | sage | 3/3 | 1.3937 | passed |
| hyper-g2-close-00001 | 50 | hyperelliptic | 3/3 | 0.2722 | passed |
| hyper-g2-close-00001 | 50 | sage | 3/3 | 1.9940 | passed |
| hyper-g2-complex | 30 | hyperelliptic | 3/3 | 0.6065 | passed |
| hyper-g2-complex | 30 | sage | 3/3 | 0.3995 | passed |
| hyper-g2-complex | 50 | hyperelliptic | 3/3 | 0.7095 | passed |
| hyper-g2-complex | 50 | sage | 3/3 | 0.6227 | passed |
| hyper-g1-linear-y | 30 | hyperelliptic | 3/3 | 0.0181 | passed |
| hyper-g1-linear-y | 30 | sage | 3/3 | 0.2153 | passed |
| hyper-g1-linear-y | 50 | hyperelliptic | 3/3 | 0.0399 | passed |
| hyper-g1-linear-y | 50 | sage | 3/3 | 0.3222 | passed |
| singular-g1-double-point | 30 | genera | 0/3 | — | not_checked |
| singular-g1-double-point | 30 | sage | 3/3 | 0.7800 | not_checked |
| singular-g1-double-point | 50 | genera | 0/3 | — | not_checked |
| singular-g1-double-point | 50 | sage | 3/3 | 0.9456 | not_checked |

## Validated timing ratios

| Case | Digits | Sage / genera | Ratio |
|---|---:|---|---:|
| hyper-g1-lemniscatic | 30 | hyperelliptic | 9.414 |
| hyper-g1-lemniscatic | 50 | hyperelliptic | 6.632 |
| hyper-g2-symmetric | 30 | hyperelliptic | 7.685 |
| hyper-g2-symmetric | 50 | hyperelliptic | 5.159 |
| hyper-g3-symmetric | 30 | hyperelliptic | 7.333 |
| hyper-g3-symmetric | 50 | hyperelliptic | 4.509 |
| hyper-g2-even | 30 | hyperelliptic | 8.816 |
| hyper-g2-even | 50 | hyperelliptic | 5.874 |
| hyper-g2-clustered | 30 | hyperelliptic | 9.570 |
| hyper-g2-clustered | 50 | hyperelliptic | 10.239 |
| hyper-g3-generic | 30 | hyperelliptic | 85.265 |
| hyper-g3-generic | 50 | hyperelliptic | 38.077 |
| hyper-g4-symmetric | 30 | hyperelliptic | 12.405 |
| hyper-g4-symmetric | 50 | hyperelliptic | 6.501 |
| general-g1-lemniscatic | 30 | genera | 0.239 |
| general-g1-lemniscatic | 50 | genera | 0.221 |
| general-trig-g3 | 30 | genera | 1.046 |
| general-trig-g3 | 50 | genera | 0.944 |
| general-trig-pure-g3 | 30 | genera | 0.675 |
| general-trig-pure-g3 | 50 | genera | 0.612 |
| general-trig-dense-g3 | 30 | genera | 1.455 |
| general-trig-dense-g3 | 50 | genera | 1.221 |
| general-klein-quartic-g3 | 30 | genera | 0.503 |
| general-klein-quartic-g3 | 50 | genera | 0.496 |
| general-fermat-quartic-g3 | 30 | genera | 0.739 |
| general-fermat-quartic-g3 | 50 | genera | 0.798 |
| general-super-g4 | 30 | genera | 0.583 |
| general-super-g4 | 50 | genera | 0.537 |
| general-kovalevskaya-g3 | 30 | genera | 0.298 |
| general-kovalevskaya-g3 | 50 | genera | 0.295 |
| nonmonic-elliptic | 30 | genera | 0.369 |
| nonmonic-elliptic | 50 | genera | 0.366 |
| general-trig-test-g3 | 30 | genera | 0.444 |
| general-trig-test-g3 | 50 | genera | 0.435 |
| general-trig-infinity-g3 | 30 | genera | 0.874 |
| general-trig-infinity-g3 | 50 | genera | 0.956 |
| klein-x-translated | 30 | genera | 0.768 |
| klein-x-translated | 50 | genera | 0.582 |
| klein-x-scaled | 30 | genera | 0.489 |
| klein-x-scaled | 50 | genera | 0.498 |
| general-mu-x-translated | 30 | genera | 2.060 |
| general-mu-x-translated | 50 | genera | 1.372 |
| hyper-g2-close-001 | 30 | hyperelliptic | 7.052 |
| hyper-g2-close-001 | 50 | hyperelliptic | 6.000 |
| hyper-g2-close-00001 | 30 | hyperelliptic | 8.800 |
| hyper-g2-close-00001 | 50 | hyperelliptic | 7.325 |
| hyper-g2-complex | 30 | hyperelliptic | 0.659 |
| hyper-g2-complex | 50 | hyperelliptic | 0.878 |
| hyper-g1-linear-y | 30 | hyperelliptic | 11.878 |
| hyper-g1-linear-y | 50 | hyperelliptic | 8.077 |

Ratio > 1 means genera was faster. Single trials are observations, not statistical estimates.
Failures, timeouts and budget skips are retained in results.json and samples.jsonl.
Accuracy tolerance is relative 10^(-digits+5), against a separately computed higher-precision Sage result.
Passing these checks is numerical evidence, not a rigorous error bound.

Source unchanged during run: True

- singular-g1-double-point (30 digits, genera): expected_rejection the polynomial must have distinct roots
- singular-g1-double-point (30 digits, genera): expected_rejection the polynomial must have distinct roots
- singular-g1-double-point (30 digits, genera): expected_rejection the polynomial must have distinct roots
- singular-g1-double-point (50 digits, genera): expected_rejection the polynomial must have distinct roots
- singular-g1-double-point (50 digits, genera): expected_rejection the polynomial must have distinct roots
- singular-g1-double-point (50 digits, genera): expected_rejection the polynomial must have distinct roots
