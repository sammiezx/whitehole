# Changelog

All notable changes to Whitehole will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned
- Kerr (rotating black hole) metric
- GPU acceleration with Numba
- Interactive Jupyter widgets
- Additional coordinate systems

## [0.1.0] - 2024

### Added

#### Core Framework
- Abstract `SpacetimeMetric` base class with full tensor algebra support
- Metric tensor computation (g_μν)
- Inverse metric computation (g^μν)
- Christoffel symbol calculation (Γ^μ_νσ)

#### Spacetime Metrics
- **Schwarzschild metric** — Static, spherically symmetric black hole
  - Eddington-Finkelstein coordinates support
  - Proper singularity and horizon handling
- **Vaidya metric** — Radiating/collapsing black hole
  - Time-dependent mass function M(v)
  - Dynamical horizon support

#### Causal Structure Computation
- Null geodesic tracer using SciPy's RK45 integrator
- Adaptive step size control
- Null constraint preservation
- Point classification algorithm for J⁻(I⁺) membership
- Event horizon detection from first principles

#### Coordinate Transformations
- Schwarzschild ↔ Kruskal-Szekeres transformation
- Kruskal ↔ Penrose compactification
- Full coordinate chain support

#### Visualization
- Penrose diagram generation
- Multi-coordinate plotting (Schwarzschild, Kruskal, Penrose side-by-side)
- Causal structure overlay on diagrams
- Horizon visualization

#### Documentation
- Comprehensive README with problem statement
- Technical architecture document
- Getting started guide with examples
- Development roadmap (6 phases)
- Stack analysis document

#### Examples
- `demo_schwarzschild.py` with 5 demonstration scenarios:
  1. Basic geodesic tracing
  2. Light cone visualization
  3. Horizon approach behavior
  4. Multi-coordinate comparison
  5. Full Penrose diagram

#### Testing
- pytest test suite setup
- Schwarzschild metric tests
- Geodesic integration tests

### Technical Details
- Python 3.10+ support
- Modern packaging with pyproject.toml and hatchling
- Type hints throughout
- NumPy-style docstrings
