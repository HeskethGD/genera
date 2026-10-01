# Benchmark Progress Log

This file tracks significant performance improvements and benchmark baseline updates for genera. Each entry documents changes that affected performance characteristics.

## Purpose

- Record performance milestones at release boundaries
- Document optimization efforts and their impact
- Provide context for comparing baselines across versions
- Track regressions or unexpected slowdowns

## Format

Each entry should include:
- **Version/Date**: Release version or date of change
- **Benchmark(s)**: Which benchmark(s) were affected
- **Baseline**: Reference to committed baseline file(s)
- **Summary**: What changed and the measured impact
- **Details**: Implementation changes, optimization techniques, or relevant context

## Entries

<!-- Entries will be added as releases occur and optimizations are made -->

### v0.1.0: Initial baseline (2026-10-01)

**Baselines:**
- `baselines/v0.1.0/curves_vs_sage.md` - Algebraic curve period computation
- `baselines/v0.1.0/hyperelliptic_vs_sage.md` - Hyperelliptic periods (isolated)
- `baselines/v0.1.0/rtheta_vs_jtheta.md` - Riemann theta vs jtheta
- `baselines/v0.1.0/rtheta_vs_wolfram.md` - Riemann theta vs Wolfram Engine

**Summary**: Established initial performance baseline for genera

**Algebraic curve periods (curves_vs_sage.md):**
- **Hyperelliptic curves**: genera is 5-10x faster than Sage across all genera (1-4) and precision levels (30-50 digits)
  - Genus 1, 30 digits: 0.018s vs 0.168s (9.4x faster)
  - Genus 2, 30 digits: 0.066s vs 0.508s (7.7x faster)
  - Genus 3, 30 digits: 0.162s vs 1.186s (7.3x faster)
  - Genus 3, 50 digits: 0.364s vs 1.640s (4.5x faster)
  - Genus 3 generic case (expensive): 0.183s vs 15.6s (85x faster!)

- **General curves**: Performance competitive with Sage
  - Trigonal genus-3: genera 2.5s vs Sage 2.6s (comparable)
  - Klein quartic: genera slower (1.9s vs 0.9s) - possible optimization target
  - Kovalevskaya: genera much slower (11.3s vs 3.4s) - uses user-supplied basis

**Accuracy**: All test cases passed accuracy validation at both 30 and 50 digit precision

**Hyperelliptic periods isolated benchmark (hyperelliptic_vs_sage.md):**
- Genus 1, 30 digits: 10ms (genera) vs 110ms (Sage) - 11x faster
- Genus 2, 30 digits: 63ms (genera) vs 539ms (Sage) - 8.5x faster
- vs abelfunctions (fixed 15-digit precision): genera 18-55x faster

**Riemann theta vs jtheta (rtheta_vs_jtheta.md):**
- Warmed steady-state: genera competitive with jtheta (0.67-1.24x ratio at 50 dps)
- Cold start: genera slower (5-18x slower) due to Python overhead vs C implementation
- Accuracy: Maximum scaled error < 1e-98 at 100 digits
- Note: jtheta is a specialized genus-1 C library; genera is pure Python and handles arbitrary genus

**Riemann theta vs Wolfram Engine (rtheta_vs_wolfram.md):**
- Overall geometric mean: genera 0.40x (2.5x faster than Wolfram)
- Unreduced cases: genera 0.03-0.16x (up to 33x faster!) - Wolfram struggles with unreduced tau
- Regular cases: genera competitive to slightly slower (0.55-1.23x)
- Genus 1, 50 digits: 0.26ms vs 0.40ms (1.5x faster)
- Genus 3 unreduced, 100 digits: 29ms vs 1028ms (35x faster!)
- Accuracy: Maximum scaled error < 6.3e-50 at 50 digits, < 2.7e-100 at 100 digits

**Infrastructure established:**
- 26+ algebraic curve test cases covering hyperelliptic and general curves
- Automated basis alignment and symplectic validation
- Isolated worker processes for fair comparison
- Riemann theta benchmarks with multiple reference implementations

The benchmark suite includes:
- **Curve benchmarks**: Comparison with SageMath for algebraic curve period computation
  - 26+ test cases covering hyperelliptic and general curves
  - Validation of differential basis, period matrices, and accuracy
- **Riemann theta benchmarks**: Comparisons with multiple reference implementations
  - Wolfram Engine (SiegelTheta)
  - jtheta implementation
  - FLINT acb_theta
  - SageMath hyperelliptic periods
  - PARI hyperelliptic periods
  - Internal implementation variants

Results tracking:
- Development runs: Git-ignored `results/TIMESTAMP/` directories
- Release baselines: Committed markdown reports in `baselines/vX.Y.Z/`
- This progress log: Manually curated performance narrative

**Next steps**: Establish v0.1.0 baseline after first release.

---

<!-- Future entries go here. Template:

### vX.Y.Z: Title of change (YYYY-MM-DD)

**Benchmark**: benchmark_name.py
**Baseline**: baselines/vX.Y.Z/benchmark_name.md
**Improvement**: 2.1x speedup for genus-3 cases at 50 digits
**Details**:

Optimized the derivative radius selection algorithm in src/genera/curves/_radius.py.
The new approach reduces the number of terms in the theta series expansion by ~40%
for high-genus cases, with minimal impact on low-genus performance.

Specific improvements:
- Genus 3, 50 digits: 8.2s → 3.9s (2.1x)
- Genus 2, 50 digits: 2.1s → 2.0s (1.05x)
- Genus 1: No change (different code path)

-->
