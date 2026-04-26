# Q1 Publication-Grade Audit (2026-04-26)

## Scope
- Static code audit of core pipeline modules in `src/jafa`.
- Test-suite verification of numerical behavior, unit handling, stitching, reprojection, ratios, CLI, and legacy import compatibility.
- Scientific-method review against implemented integration and uncertainty propagation logic.

## What was checked
1. Full test suite (`pytest -q`).
2. Lint/static checks (`ruff check src tests`).
3. Numerical integration and uncertainty propagation paths in `src/jafa/mapping.py`.
4. Ratio-map masking and uncertainty propagation in `src/jafa/ratios.py`.
5. Reprojection uncertainty behavior in `src/jafa/reprojection.py`.
6. Cube stitching merge behavior in `src/jafa/stitching.py`.

## Findings summary
- **Critical bugs:** none identified in audited paths.
- **High severity:** none identified in audited paths.
- **Medium severity:** one maintainability/forward-compatibility issue fixed (deprecated NumPy integration API usage).
- **Consistency status:** package rename transition appears consistent (`jafa` primary with `jwst_feature_mapper` compatibility layer).
- **Scientific accuracy status (audited paths):** integration and error propagation choices are internally consistent with the package’s documented flux-density assumptions.

## Fixed in this audit
### 1) Deprecated numerical integration function
- Replaced `np.trapz` with `np.trapezoid` in spectral integration helpers.
- Why this matters: avoids NumPy deprecation warnings and future break risk without changing integration semantics.

## Residual risks / recommendations before publication
1. Add an explicit CI gate that fails on new deprecations (e.g., `-W error::DeprecationWarning` in a dedicated job) to prevent regression.
2. Add a short “assumptions and validity limits” section in docs for:
   - expected flux density convention (`F_\nu`) for cgs output integration,
   - assumptions behind uncertainty propagation (independent errors, no covariance term).
3. Add one regression test asserting numerical equivalence of integration output across sorted/unsorted wavelength grids for both native and cgs modes.

## Publication readiness statement (Q1)
Based on this audit pass, the repository appears **publication-ready for Q1** for the audited functionality, with no blocking correctness defects observed. Remaining recommendations are quality-hardening steps rather than release blockers.
