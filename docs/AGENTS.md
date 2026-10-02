# Documentation aims

- Make literature examples readable as mathematical walkthroughs in the
  online documentation. Connect the mathematical problem, conventions,
  Genera's public API and numerical validation, with clear attribution.
- Add the papers used by each walkthrough to `references.rst` and cite those
  shared entries from the example page. Verify bibliographic details against
  the paper, its publisher or its arXiv record; reuse existing citation entries.
- Explain what each example actually computes and checks, including any
  supplied numerical data or parts of the paper that are outside its scope.
- Keep build mechanics, contributor instructions and runtime selection
  policies here rather than in the reader-facing Examples index. Individual
  pages should explain their numerical accuracy and practical limitations.

## Presentation

- Display groups of equations in a single vertical column, with one equation
  per row. Use an `aligned` or `gathered` environment; do not place independent
  equations side by side with `\qquad` or multiple equation columns.
- Preserve the mathematical layout of matrices, vectors and single equations.
- Do not mention temporary local folders or other developer-only workspace
  details in reader-facing documentation.
- Use plots when they help explain the example. For trajectory comparisons,
  use dense RK4 samples drawn as dashed curves and sparse analytic samples
  drawn as markers. Check reality before displaying real components.
- Explain residuals in prose with representative values and the checks being
  enforced. Do not add distracting plots of numerical error or roundoff noise.

## Executable examples

- Reuse tracked runnable code for plots and doctests, preferably through one
  checked helper shared by both. Documentation must build without local
  reference files or temporary data folders.
- Include source excerpts directly from the scripts where useful. Ordinary
  displayed code is explanatory; use executable doctest blocks for the
  calculation and meaningful assertions that keep it in sync with Genera.
- Plots execute during Sphinx HTML builds. Executable `>>>` blocks are checked
  by pytest's documentation doctests and Sphinx's doctest builder. HTML and
  doctest builds execute independently, so account for both costs.

## Runtime and numerical accuracy

- Aim for less than 15 seconds locally per documentation calculation, before
  plot rendering. This is an initial selection budget, not a portable timing
  guarantee. Existing example subprocess tests allow 120 seconds for slower CI.
- Measure analytic evaluation and RK4 costs when choosing an interval,
  precision and sampling grid. Reuse period data and batch evaluations where
  appropriate. Keep expensive variants in standalone examples; defer examples
  that cannot meet a useful accuracy within the documentation budget.
- When RK4 is cheap relative to analytic evaluation, reduce its step size
  independently of the sparse analytic plotting grid. Start RK4 from the
  analytic initial state and then integrate independently.
- Check step refinement and conserved quantities where appropriate. If
  rounded input data impose an accuracy floor, check RK4 convergence between
  numerical grids rather than requiring the analytic comparison error to
  continue shrinking.
- Distinguish arbitrary working precision, RK4 discretization error and the
  precision of external reference data. Choose tolerances for each check from
  those limits; increased precision cannot recover missing reference digits.
  Analytic identities may justify much tighter tolerances than ODE comparisons.
- State whether residuals are absolute or scaled, describe what they validate,
  and make clear that sampled checks do not bound the continuous trajectory.

## Validation workflow

- Run focused doctests and the corresponding standalone example smoke tests
  when changing executable examples. Do not run full test suites for this work
  unless the user explicitly requests them.
- Build HTML with warnings treated as errors, inspect newly added or changed
  plots, and run relevant lint and diff checks. For prose-only changes, use
  the documentation build and diff checks without rerunning numerical tests.
