# Safe Benchmark Execution Guide

This document explains how to safely run generapy benchmarks in environments with multiple mathematical tools (Sage, Wolfram Engine, PARI, FLINT) while avoiding version conflicts and ensuring proper isolation.

## Core Principle: Process Isolation

**Every external tool runs in a separate, isolated process** to prevent:
- Python namespace conflicts
- mpmath/generapy version mismatches
- Shared state corruption
- Cache interference

## Critical Safety Rules

### 1. Never Mix Imports

**WRONG** - Never do this:
```python
import generapy
import sage.all  # DANGER: Namespace collision!
```

**RIGHT** - Use subprocess isolation:
```python
# generapy runs in this process
import generapy

# Sage runs in a separate worker process
result = subprocess.run(['sage', '-python', 'worker.py'], ...)
```

All benchmarks follow this pattern. The main process uses generapy; external tools run in isolated workers.

### 2. Sage-Specific Isolation

**Critical**: Sage must NEVER import generapy or mpmath from the development checkout.

**Why?** Sage bundles its own mpmath version. Loading a different mpmath version causes:
- Unpredictable numerical behavior
- Import errors
- Silent computation errors

**How benchmarks handle this:**

```python
# In benchmark_curves.py (main process)
import generapy  # Uses development version

# Worker launched with clean environment
env = dict(os.environ)
env.pop('PYTHONPATH', None)  # Remove generapy from path
env.pop('PYTHONHOME', None)

# Sage worker runs in temporary directory OUTSIDE generapy checkout
with tempfile.TemporaryDirectory(prefix='curve-benchmark-') as cwd:
    subprocess.Popen(['sage', '-python', 'worker.py'], cwd=cwd, env=env, ...)
```

**Never:**
- Add generapy to Sage's `PYTHONPATH`
- Run Sage workers from inside the generapy repository
- Import generapy inside a Sage session

### 3. Working Directory Isolation

**Sage workers must start in a temporary directory outside the generapy checkout.**

From `benchmark_curves.py`:
```python
# Sage MUST start outside the development checkout. Use the same isolated
# cwd for both workers, and leave Sage's normal HOME/cache configuration alone.
with tempfile.TemporaryDirectory(prefix='curve-benchmark-') as cwd:
    process = subprocess.Popen(command, cwd=cwd, ...)
```

**Why?** Prevents Sage from accidentally finding and importing generapy/mpmath from the current directory.

### 4. Wolfram Engine Isolation

Wolfram Engine communicates via subprocess and JSON:

```python
# Send input as Mathematica code
process = subprocess.Popen(['/usr/local/bin/wolframscript', '-code', wolfram_code],
                           stdout=subprocess.PIPE, ...)
# Parse output
result = parse_wolfram_output(process.stdout)
```

No shared Python namespace, so isolation is straightforward.

### 5. PARI/GP and FLINT Isolation

Similar to Wolfram: external processes communicating via stdin/stdout.

## Environment Variables

### Clean Environment for Workers

All workers receive a sanitized environment:

```python
env = dict(os.environ)
env.pop('PYTHONPATH', None)    # Don't inherit Python paths
env.pop('PYTHONHOME', None)    # Don't inherit Python home
# Sage keeps its normal HOME for ~/.sage/cache access
```

### Sage Cache Permissions

Sage may need write access to `~/.sage/cache`. If running in a sandbox:
- Grant permission to write to the user's home directory, OR
- Set `SAGE_CACHE` to a writable location

## Verifying Correct Isolation

### Check generapy Version

```python
# In main process
import generapy
print(generapy.__version__)  # Should be your development version
```

### Check Sage is Using Its Own mpmath

```python
# Inside a Sage worker
from sage.all import *
import mpmath
print(mpmath.__file__)  # Should be inside Sage installation, NOT generapy
```

If the path points to your generapy checkout, **isolation has failed**.

## Concurrent Execution Prevention

The curve benchmarks use a POSIX advisory lock to prevent multiple instances:

```python
# In benchmark_curves.py
lockfile = Path('.run.lock')
handle = lockfile.open('a+')
fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)  # Fail if already locked
```

**Never run multiple benchmark instances from the same directory concurrently.** This prevents:
- Race conditions writing results
- Confused result attribution
- Lock contention

## Version Tracking

Each benchmark records:
- generapy/mpmath version and git commit
- External tool versions (Sage, Wolfram Engine, etc.)
- System information (platform, Python version)

From `results.json`:
```json
{
  "versions": {
    "generapy": "0.1.0.dev0+g1a2b3c4",
    "sage": "10.1",
    "python": "3.11.5"
  },
  "git_state": {
    "commit": "1a2b3c4...",
    "status": "M src/generapy/curves/periods.py"
  }
}
```

**Warning**: Results from modified working trees (`status` not empty) should not be used for official baselines.

## Common Pitfalls

### ❌ Running benchmarks with uncommitted changes
**Risk**: Results don't correspond to any committed version.
**Solution**: Commit changes first, or document as experimental.

### ❌ Modifying code during benchmark runs
**Risk**: Different cases measured with different code.
**Solution**: Don't edit code while benchmarks run. The harness detects source changes and marks results invalid.

### ❌ Putting generapy on Sage's PYTHONPATH
**Risk**: Sage imports wrong mpmath, all Sage results invalid.
**Solution**: Never set PYTHONPATH. Let Sage use its bundled mpmath.

### ❌ Running Sage workers from generapy directory
**Risk**: Python imports from current directory, namespace collision.
**Solution**: Workers automatically use temporary directories. Don't override `cwd`.

### ❌ Comparing results across different systems
**Risk**: Performance differences due to hardware, not code.
**Solution**: Always compare on the same machine. Document system specs in baselines.

### ❌ Running benchmarks in CI without external tools
**Risk**: CI tries to run Sage/Wolfram benchmarks without installations.
**Solution**: Use `--engines generapy` for generapy-only benchmarks in CI.

## Results Management

### What to Commit

**✅ Commit:**
- `baselines/vX.Y.Z/*.md` - Release baseline reports (8-20KB each)
- `PROGRESS.md` - Curated improvement log

**❌ Never commit:**
- `results/TIMESTAMP/` - Timestamped development runs (git-ignored)
- `results.json` - Large files with full period matrices (~1MB)
- `samples.jsonl` - Raw worker outputs (~700KB)
- References to local-only timestamped directories

### Baseline Creation Checklist

Before committing a release baseline:

1. ✅ Clean working tree (no uncommitted changes)
2. ✅ Run comprehensive benchmarks (`--preset all`)
3. ✅ Verify all accuracy checks passed
4. ✅ Review `report.md` for anomalies
5. ✅ Copy ONLY the markdown report to `baselines/vX.Y.Z/`
6. ✅ Update `PROGRESS.md` with summary
7. ✅ Commit with descriptive message

**Never:**
- Commit from modified working tree
- Cherry-pick results (run full suite)
- Edit generated reports manually
- Reference local timestamped directories in committed files

## Timeout and Resource Limits

### Worker Timeouts

```bash
--timeout 120  # Kill worker after 120 seconds (includes process startup)
```

Timeouts prevent:
- Hanging on pathological cases
- Infinite loops
- Resource exhaustion

Killed workers are recorded as timeouts, not successes.

### Budget Limits

```bash
--budget 1800  # Stop launching new workers after 30 minutes total
```

The budget is a **ceiling**, not a guarantee. Expensive cases may be skipped when budget runs out. Skipped cases appear as `skipped_budget` in results.

### Process Cleanup

All workers run in their own process group. On timeout or Ctrl-C:
```python
os.killpg(process.pid, signal.SIGKILL)  # Kill entire process group
process.communicate()  # Clean up pipes
```

This ensures:
- No zombie processes
- No orphaned Sage/Wolfram workers
- Clean shutdown on interrupt

## Platform Notes

### macOS/Linux
All benchmarks target POSIX systems. Advisory locks require POSIX `fcntl`.

### Windows
The benchmark harness has not been tested on Windows. Process group management and file locking may require adaptation.

## Debugging

### Enable verbose output

Most benchmarks accept `--verbose` or similar flags. Check individual script `--help`.

### Inspect worker failures

```bash
# Check stderr from failed worker
cat benchmarks/results/LATEST/samples.jsonl | grep '"status":"worker_exit_error"'
```

### Verify isolation manually

```python
# Launch a Sage worker manually
import subprocess
import tempfile

with tempfile.TemporaryDirectory() as cwd:
    proc = subprocess.run(['sage', '-c', 'import mpmath; print(mpmath.__file__)'],
                          cwd=cwd, capture_output=True, text=True)
    print(proc.stdout)  # Should be inside Sage installation
```

## Summary

**Key takeaways:**
1. External tools run in isolated subprocess workers, never in the main process
2. Sage workers use temporary directories outside generapy checkout
3. Clean environment variables prevent PYTHONPATH contamination
4. One benchmark instance per directory (advisory lock)
5. Only commit aggregated markdown reports, never raw JSON
6. Verify clean working tree before baseline creation
7. Never reference local-only timestamped directories in committed files

Following these guidelines ensures benchmarks are reproducible, safe, and scientifically valid.
