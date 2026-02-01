# Whitehole

**A Computational Framework for Global Causal Structure Visualization in General Relativity**

## What is Whitehole?

Whitehole is a research-grade tool that bridges the gap between theoretical general relativity and computational visualization. It computes and visualizes global causal structure from spacetime data — making event horizons, causal boundaries, and the fundamental structure of spacetime tangible and explorable.

## The Problem We Solve

| Current State | Whitehole Solution |
|---------------|-------------------|
| Penrose diagrams drawn manually, after the fact | Dynamically computed from spacetime data |
| Horizons assumed to be known | Horizons discovered from first principles |
| Numerical relativity outputs metrics without causal insight | Full causal structure exposed intuitively |
| Static, non-interactive visualizations | Interactive, falsifiable exploration |

## Core Principles

### 1. Event Horizons Are Global Objects
- Defined as ∂J⁻(I⁺) — the boundary of the causal past of future null infinity
- Cannot be found locally; must be computed retrospectively
- Our tool makes this nonlocality computationally explicit

### 2. Causality Over Geometry
- Light cones define structure; distances are secondary
- Horizons are causal boundaries, not geometric ones
- Everything flows from null geodesic propagation

### 3. Horizons ≠ Apparent Horizons
- Event horizons depend on the entire future evolution
- Apparent horizons depend on spatial slicing
- We show the difference explicitly

### 4. Visualization Respects Invariance
- No misleading "space at a time" plots by default
- Penrose compactified, conformal diagrams, causal diamonds
- Multiple coordinate systems shown side-by-side

## Architecture Overview

```
┌────────────────────────────────────┐
│  Visualization & Interaction Layer │  ← Penrose, Kruskal, Causal Diamonds
├────────────────────────────────────┤
│  Causal Structure Computation Core │  ← Null Geodesics, J⁻(I⁺), Horizons
├────────────────────────────────────┤
│  Spacetime Input & Evolution Layer │  ← Metrics (Analytic & Numerical)
└────────────────────────────────────┘
```

## Supported Spacetimes

**Analytic Metrics:**
- Schwarzschild (static black hole)
- Vaidya (radiating collapse)
- Oppenheimer-Snyder (dust collapse)
- Kerr (rotating black hole — planned)

**Numerical Input:**
- Metric tensor g_μν(t, x) from numerical relativity codes
- Synthetic collapse models

## Getting Started

See [docs/GETTING_STARTED.md](docs/GETTING_STARTED.md) for installation and first steps.

## Documentation

- [Technical Architecture](docs/ARCHITECTURE.md) — Deep dive into system design
- [Getting Started](docs/GETTING_STARTED.md) — Installation and first visualization
- [Development Roadmap](docs/ROADMAP.md) — Phased execution plan

## Research Extensions

This project is designed to support genuine research:

- **Horizon vs Apparent Horizon Mismatch** — Compare true event horizons with apparent horizons under different slicings
- **Dynamical Horizons** — Marginally trapped surfaces, Hayward horizons, isolated horizons
- **QFT Overlay** — Mode propagation, near-horizon redshift, Hawking radiation precursors
- **Causal Set Approximation** — Discretize spacetime and compare with continuum GR

## Why This Matters

> "This project shows deep conceptual understanding, technical skill, and produces something genuinely useful to the GR community."

Whitehole is designed to be:
- A thesis-level research tool
- A publishable computational physics platform
- An open research framework for causal structure and horizon formation

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to get started.

## Citation

If you use Whitehole in your research, please cite:

```bibtex
@software{whitehole,
  title = {Whitehole: A Computational Framework for Global Causal Structure Visualization in General Relativity},
  author = {Sameer},
  year = {2024},
  url = {https://github.com/sammiezx/whitehole}
}
```
