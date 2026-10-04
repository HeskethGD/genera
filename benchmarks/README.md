# Generapy Benchmarks

Performance benchmarks for generapy's algebraic curve and Abelian function implementations. These benchmarks compare generapy against other symbolic/numeric mathematics systems but are **not part of the generapy package** and are **not required for using generapy**.

## Directory Structure

```
benchmarks/
├── README.md           # This file
├── AGENTS.md          # Safety guide for running in different environments
├── PROGRESS.md        # Curated log of performance improvements
├── curves/            # Algebraic curve period benchmarks
├── rtheta/            # Riemann theta function benchmarks
├── results/           # Git-ignored timestamped results
└── baselines/         # Git-tracked release baselines
    └── vX.Y.Z/        # Versioned baseline results (added at releases)
```

## Prerequisites

### Required
- Python 3.10+
- generapy installed in development mode: `.venv/bin/pip install -e '.[develop]'`

### Optional (for specific benchmarks)
- **Wolfram Engine** (free for developers): Required for `benchmark_rtheta_vs_wolfram.py`
  - Download from https://www.wolfram.com/engine/
  - Default path: `/usr/local/bin/wolframscript`
- **SageMath**: Required for `benchmark_*_vs_sage.py` benchmarks
  - Download from https://www.sagemath.org/
  - Default path: `/usr/local/bin/sage`
- **PARI/GP**: Required for `benchmark_hyperelliptic_vs_pari.py`
  - Install via package manager or from https://pari.math.u-bordeaux.fr/
- **FLINT**: Required for `benchmark_rtheta_vs_flint.py`
  - Install via package manager

**Note**: generapy itself does not require any of these tools. They are only needed for comparative benchmarking.

## Quick Start

Verify your setup with a minimal smoke test:

```bash
# From the generapy repository root
cd benchmarks/curves
../../.venv/bin/python benchmark_curves.py --preset smoke --digits 20 --trials 1
```

This runs two representative curve cases and validates the benchmark harness without requiring external tools.

## Running Benchmarks

All commands should be run from the **generapy repository root**.

### Algebraic Curve Benchmarks

The curve benchmarks compare generapy's period computation against SageMath.

**List available test cases:**
```bash
.venv/bin/python benchmarks/curves/benchmark_curves.py --list
```

**Quick validation (6 cases, ~5 minutes):**
```bash
.venv/bin/python benchmarks/curves/benchmark_curves.py \
  --preset quick --digits 30 --trials 2 --budget 300 --timeout 60
```

**Comprehensive benchmark:**
```bash
.venv/bin/python benchmarks/curves/benchmark_curves.py \
  --preset all --digits 30 50 --trials 3 --budget 1800 --timeout 120
```

**Run specific cases:**
```bash
.venv/bin/python benchmarks/curves/benchmark_curves.py \
  --cases hyper-g2-symmetric general-klein-quartic-g3 \
  --digits 30 --trials 3
```

**Available presets:**
- `smoke`: Two quick cases for validation
- `quick`: Six representative cases
- `established`: Standard test suite
- `expanded`: Extended test cases
- `all`: All 26+ test cases including stress tests

**Key options:**
- `--digits 20 30 50`: Target precision(s)
- `--trials 3`: Number of samples per configuration
- `--budget 1800`: Total time budget in seconds
- `--timeout 120`: Per-worker timeout
- `--engines mpmath sage`: Which implementations to test (default: both)
- `--sage /path/to/sage`: Custom Sage path

### Riemann Theta Benchmarks

**Compare against Wolfram Engine:**
```bash
.venv/bin/python benchmarks/rtheta/benchmark_rtheta_vs_wolfram.py \
  --digits 50 100 --trials 3
```

**Compare against jtheta implementation:**
```bash
.venv/bin/python benchmarks/rtheta/benchmark_rtheta_vs_jtheta.py
```

**Compare hyperelliptic periods with Sage:**
```bash
.venv/bin/python benchmarks/rtheta/benchmark_hyperelliptic_vs_sage.py \
  --digits 30 50 --trials 3
```

**Compare hyperelliptic periods with PARI:**
```bash
.venv/bin/python benchmarks/rtheta/benchmark_hyperelliptic_vs_pari.py
```

**Compare against FLINT:**
```bash
.venv/bin/python benchmarks/rtheta/benchmark_rtheta_vs_flint.py
```

**Internal implementation comparison:**
```bash
.venv/bin/python benchmarks/rtheta/benchmark_rtheta_implementations.py
```

## Understanding Results

### Output Location

All benchmark runs create timestamped directories in `benchmarks/results/`:

```
results/20261015T143022.123456Z/
├── report.md         # Human-readable summary
├── results.json      # Full machine-readable data
└── samples.jsonl     # Raw worker outputs
```

The timestamp format is ISO 8601 UTC: `YYYYMMDDTHHMMSSμμμμμμZ`

### Reading Reports

The `report.md` file contains:
- Timing summaries (median seconds per case)
- Success/failure rates per configuration
- Accuracy validation status
- Speed ratios between implementations

Example output:
```
| Case | Digits | Implementation | Successes/trials | Median seconds | Accuracy |
|---|---:|---|---:|---:|---|
| hyper-g2-symmetric | 30 | generapy | 3/3 | 0.0672 | passed |
| hyper-g2-symmetric | 30 | sage | 3/3 | 0.5307 | passed |
```

Speed ratios > 1.0 indicate generapy was faster than the reference implementation.

## Results Tracking

### Development Workflow

During development, all results go to git-ignored `results/TIMESTAMP/` directories:

```bash
# Run benchmark
.venv/bin/python benchmarks/curves/benchmark_curves.py --preset quick

# View latest results
cat benchmarks/results/$(ls -t benchmarks/results/ | head -1)/report.md
```

### Release Baselines

At release milestones, copy the markdown report to `baselines/`:

```bash
# After running comprehensive benchmarks for v0.1.0
mkdir -p benchmarks/baselines/v0.1.0
cp benchmarks/results/LATEST/report.md \
   benchmarks/baselines/v0.1.0/curves_vs_sage.md

# Document the baseline
echo "## v0.1.0: Initial release baseline" >> benchmarks/PROGRESS.md

git add benchmarks/baselines/v0.1.0/ benchmarks/PROGRESS.md
git commit -m "Add v0.1.0 benchmark baselines"
```

**Only commit the markdown reports** - the JSON files are too large and contain low-level data not needed for tracking progress.

### Progress Log

`PROGRESS.md` is a manually curated log documenting significant performance changes:

```markdown
## v0.2.0: Genus-3 optimization
**Improvement:** 2.1x speedup for genus-3 theta functions at 50 digits
**Baseline:** baselines/v0.2.0/rtheta_vs_wolfram.md
**Details:** Optimized derivative radius selection
```

## Safety and Environment Isolation

**IMPORTANT**: Read [AGENTS.md](AGENTS.md) before running benchmarks that involve external tools (Sage, Wolfram Engine, etc.).

Key points:
- Each external tool runs in an isolated subprocess
- Sage workers run in temporary directories outside the generapy checkout
- Never add generapy to Sage's `PYTHONPATH`
- Avoid running benchmarks concurrently from the same directory

## Troubleshooting

**Sage not found:**
```bash
# Specify custom Sage path
.venv/bin/python benchmarks/curves/benchmark_curves.py --sage /path/to/sage
```

**Wolfram Engine not found:**
Edit the `WOLFRAMSCRIPT` path in the benchmark script, or install to the default location.

**Worker timeout:**
Increase the timeout: `--timeout 300` (5 minutes)

**Budget exceeded:**
Some expensive cases may be skipped. Either:
- Increase budget: `--budget 3600`
- Run fewer cases with specific `--cases`
- Reduce `--digits` or `--trials`

**Memory issues:**
Avoid high-genus curves at very high precision. Start with `--preset quick` to validate setup.

## Contributing Benchmarks

When adding new benchmarks:
1. Follow the existing isolation patterns (see `period_worker.py`)
2. Generate timestamped output in `results/`
3. Produce both JSON and markdown reports
4. Add tests to validate the harness (see `test_benchmark.py`)
5. Document new external dependencies in this README
6. Update `AGENTS.md` with any new safety considerations

## Benchmark Validation

Test the benchmark harness itself:

```bash
.venv/bin/python benchmarks/curves/test_benchmark.py
```

This validates basis conversion, symplectic checks, and worker isolation without running full benchmarks.
