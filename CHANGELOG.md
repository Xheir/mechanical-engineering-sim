# Changelog

## 0.2.0 — 2026-11-08

### Added

- Independent Cook's membrane benchmark comparison against the published CalculiX reference.
- Reproducible mesh-convergence data and a regression test for the benchmark.
- Read-only GitHub Actions verification for the full four-test suite.

### Security

- Updated the development test dependency to `pytest==9.0.3` to address the published temporary-directory handling advisory.

## 0.1.0 — 2026-11-08

### Added

- Linear-elastic STEP-to-FEA example with JSON report and deformed-shape plot.
- Slider-crank and four-bar kinematics solver with numerical verification gates.
- Kerf-compensated sheet-metal flat-pattern DXF generator with self-verification.
- Pinned dependency baseline, safety boundary, limitations, and MIT license.
