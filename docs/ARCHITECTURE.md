# Technical Architecture

This document provides a deep dive into Whitehole's three-layer architecture and the mathematical foundations underlying each component.

## System Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                VISUALIZATION & INTERACTION LAYER                │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────────────────┐  │
│  │   Penrose   │  │   Kruskal   │  │   Interactive Controls  │  │
│  │   Diagram   │  │   Diagram   │  │   - Point selection     │  │
│  └─────────────┘  └─────────────┘  │   - Ray launching       │  │
│  ┌─────────────┐  ┌─────────────┐  │   - Parameter tuning    │  │
│  │   Causal    │  │ Schwarzschild│  └─────────────────────────┘  │
│  │   Diamond   │  │  Coordinates │                              │
│  └─────────────┘  └─────────────┘                              │
├─────────────────────────────────────────────────────────────────┤
│                CAUSAL STRUCTURE COMPUTATION CORE                │
│  ┌───────────────────┐  ┌───────────────────────────────────┐  │
│  │  Null Geodesic    │  │  Causal Past/Future Constructor   │  │
│  │  Engine           │  │  J⁻(I⁺) numerical computation     │  │
│  │  - Integration    │  └───────────────────────────────────┘  │
│  │  - Ray tracing    │  ┌───────────────────────────────────┐  │
│  └───────────────────┘  │  Horizon Detector                 │  │
│  ┌───────────────────┐  │  - Event horizon: ∂J⁻(I⁺)         │  │
│  │  Christoffel      │  │  - Apparent horizon (optional)    │  │
│  │  Symbols Γ^μ_νσ   │  └───────────────────────────────────┘  │
│  └───────────────────┘                                         │
├─────────────────────────────────────────────────────────────────┤
│                SPACETIME INPUT & EVOLUTION LAYER                │
│  ┌─────────────────────────┐  ┌─────────────────────────────┐  │
│  │   Analytic Metrics      │  │   Numerical Metrics         │  │
│  │   - Schwarzschild       │  │   - g_μν(t, x) arrays       │  │
│  │   - Vaidya              │  │   - Interpolation engine    │  │
│  │   - Oppenheimer-Snyder  │  │   - NR code interface       │  │
│  │   - Kerr                │  └─────────────────────────────┘  │
│  └─────────────────────────┘                                   │
└─────────────────────────────────────────────────────────────────┘
```

---

## Layer 1: Spacetime Input & Evolution Layer

### Purpose
Define the spacetime manifold (M, g_μν) as raw data that feeds into causal structure computation.

### 1.1 Analytic Metric Modules

Each analytic spacetime is implemented as a module providing:

```
interface SpacetimeMetric {
    // Metric tensor components
    g_μν(x: Coordinates): MetricTensor

    // Inverse metric (for raising indices)
    g^μν(x: Coordinates): MetricTensor

    // Christoffel symbols (for geodesic equation)
    Γ^μ_νσ(x: Coordinates): ChristoffelSymbols

    // Optional: analytic derivatives for accuracy
    ∂_ρ g_μν(x: Coordinates): MetricDerivatives
}
```

#### Schwarzschild Metric
```
ds² = -(1 - 2M/r)dt² + (1 - 2M/r)⁻¹dr² + r²dΩ²

Parameters: M (mass)
Singularity: r = 0
Event horizon: r = 2M
```

#### Vaidya Metric (Ingoing)
```
ds² = -(1 - 2M(v)/r)dv² + 2dvdr + r²dΩ²

Parameters: M(v) (mass function of advanced time)
Key feature: Models radiating/accreting black holes
Use case: Horizon formation during collapse
```

#### Oppenheimer-Snyder Metric
```
Interior: FLRW dust collapse
Exterior: Schwarzschild
Matching: Junction conditions at dust surface

Key feature: Complete collapse model
Use case: Visualize horizon formation from collapsing matter
```

### 1.2 Numerical Metric Interface

For numerical relativity data:

```
interface NumericalSpacetime {
    // Grid structure
    grid: SpacetimeGrid

    // Metric data at grid points
    metric_data: Array<MetricTensor>

    // Interpolation for off-grid points
    interpolate(x: Coordinates): MetricTensor

    // Christoffel computation (numerical derivatives)
    compute_christoffels(x: Coordinates): ChristoffelSymbols
}
```

**Interpolation Methods:**
- Trilinear (fast, lower accuracy)
- Cubic spline (balanced)
- Spectral (high accuracy, computationally intensive)

---

## Layer 2: Causal Structure Computation Core

This is the mathematical heart of the system.

### 2.1 Null Geodesic Engine

#### The Geodesic Equation

For null geodesics (light rays):
```
d²xᵘ/dλ² + Γᵘ_νσ (dxᵛ/dλ)(dxᵠ/dλ) = 0

with null constraint:
g_μν (dxᵘ/dλ)(dxᵛ/dλ) = 0
```

#### Implementation Strategy

**State Vector:**
```
Y = [x⁰, x¹, x², x³, k⁰, k¹, k², k³]

where xᵘ = position, kᵘ = dxᵘ/dλ = tangent vector
```

**First-Order System:**
```
dxᵘ/dλ = kᵘ
dkᵘ/dλ = -Γᵘ_νσ kᵛ kᵠ
```

**Integration Methods:**
- RK4 (4th order Runge-Kutta) — default
- RK45 (adaptive step) — for accuracy-critical regions
- Symplectic integrators — for long-term evolution

**Null Constraint Enforcement:**
After each step, project kᵘ back onto null cone:
```
kᵘ → kᵘ - (g_μν kᵘ kᵛ / 2g_00) δᵘ_0
```

### 2.2 Causal Past/Future Construction

#### Computing J⁻(I⁺) — The Causal Past of Future Null Infinity

**Algorithm:**

1. **Discretize spacetime** into a grid of test points {p_i}

2. **For each point p_i:**
   - Launch null geodesics in multiple directions (angular discretization)
   - Integrate forward in time (toward the future)
   - Check termination conditions:
     - **Escapes:** r → ∞ at late times (finite affine parameter)
     - **Trapped:** r → 0 or r → r_horizon

3. **Classification:**
   - If ANY null ray from p escapes: p ∈ J⁻(I⁺)
   - If ALL null rays are trapped: p ∉ J⁻(I⁺)

4. **Boundary detection:**
   - Event horizon H⁺ = ∂J⁻(I⁺)
   - Points where escape/trap classification changes

```
function classify_point(p: SpacetimePoint) -> CausalClassification:
    rays = generate_null_directions(p, angular_resolution)
    escaping_rays = 0

    for ray in rays:
        trajectory = integrate_null_geodesic(p, ray, max_affine_param)
        if reaches_large_r(trajectory):
            escaping_rays += 1

    if escaping_rays > 0:
        return CAUSAL_PAST_OF_SCRI_PLUS
    else:
        return TRAPPED_REGION
```

### 2.3 Retrodictive Horizon Construction

The event horizon is found by working **backwards from the future.**

**Algorithm:**

1. At late times, identify the **apparent** boundary between escaping and trapped regions

2. Take null generators (null geodesics) that hover at this boundary

3. Trace these generators **backwards in time**

4. The surface they sweep out is the event horizon

**Key Insight:** This makes the global, teleological nature of event horizons computationally explicit. You cannot find the horizon without knowing the future.

### 2.4 Apparent Horizon Detection (Optional)

Unlike event horizons, apparent horizons are **local** — they depend on the choice of spatial slice.

**Definition:** The apparent horizon is the outermost marginally outer trapped surface (MOTS) on a spatial slice.

**Detection:**
- On each spatial slice, find surfaces where outgoing null expansion θ_+ = 0
- Use techniques from numerical relativity (e.g., fast flow methods)

**Comparison Module:**
- Overlay event horizon and apparent horizon
- Show how they differ during dynamical evolution
- Demonstrate slice-dependence of apparent horizons

---

## Layer 3: Visualization & Interaction Layer

### 3.1 Coordinate Transformations

#### Schwarzschild → Kruskal-Szekeres
```
For r > 2M:
  U = √(r/2M - 1) exp(r/4M) cosh(t/4M)
  V = √(r/2M - 1) exp(r/4M) sinh(t/4M)

For r < 2M:
  U = √(1 - r/2M) exp(r/4M) sinh(t/4M)
  V = √(1 - r/2M) exp(r/4M) cosh(t/4M)
```

#### Kruskal → Penrose (Compactification)
```
U' = arctan(U)
V' = arctan(V)

Maps infinite spacetime to finite diagram
45° lines remain null
```

### 3.2 Visualization Modes

#### Penrose Diagram (Default)
- Compactified representation
- Null rays at 45°
- Infinity brought to finite distance
- Ideal for global causal structure

#### Kruskal Diagram
- Non-compact but regular at horizon
- Good for near-horizon physics
- Shows maximal extension

#### Schwarzschild Coordinates
- Traditional (t, r) plot
- Warning: coordinate singularity at horizon
- Useful for comparison with textbooks

#### Causal Diamond View
- Focus on a specific region
- Show J⁺(p) ∩ J⁻(q) for chosen events
- Interactive exploration

### 3.3 Interactive Features

**Point Selection:**
- Click to select spacetime event p
- Display coordinates in all systems
- Show local light cone

**Ray Launching:**
- From selected point, launch null rays
- Visualize trajectory in real-time
- Color-code: escaping (blue) vs trapped (red)

**Parameter Tuning:**
- Adjust black hole mass M
- Modify Vaidya mass function M(v)
- Watch horizon respond dynamically

**Animation:**
- Animate collapse process
- Show horizon formation
- Time-slider for exploration

---

## Data Flow Summary

```
User Input (parameters, interactions)
           │
           ▼
┌─────────────────────────┐
│  Spacetime Layer        │
│  → Compute g_μν, Γ^μ_νσ │
└───────────┬─────────────┘
            │
            ▼
┌─────────────────────────┐
│  Causal Core            │
│  → Trace null geodesics │
│  → Classify points      │
│  → Find horizons        │
└───────────┬─────────────┘
            │
            ▼
┌─────────────────────────┐
│  Visualization Layer    │
│  → Transform coords     │
│  → Render diagrams      │
│  → Handle interaction   │
└─────────────────────────┘
            │
            ▼
      Visual Output
```

---

## Performance Considerations

### Parallelization Opportunities
- Null geodesic integration is embarrassingly parallel
- Each ray is independent
- GPU acceleration via compute shaders

### Caching Strategies
- Cache Christoffel symbols on grid
- Memoize coordinate transformations
- Incremental updates for parameter changes

### Accuracy vs Speed Tradeoffs
- Adaptive integration near horizons
- Coarse-to-fine point classification
- LOD (level of detail) for visualization

---

## Future Architecture Extensions

### Kerr Spacetime Support
- Boyer-Lindquist coordinates
- Ergosphere visualization
- Frame dragging effects

### Numerical Relativity Interface
- HDF5 data import
- Support for common NR codes (Einstein Toolkit, SpEC)
- Time interpolation for animations

### QFT Overlay Module
- Scalar field mode propagation
- Near-horizon redshift visualization
- Hawking temperature computation
