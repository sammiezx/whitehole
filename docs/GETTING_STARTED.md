# Getting Started with Whitehole

This guide walks you through setting up Whitehole and creating your first causal structure visualization.

## Prerequisites

### Mathematical Background
Before diving in, ensure familiarity with:
- General Relativity basics (metrics, geodesics, horizons)
- Coordinate systems (Schwarzschild, Kruskal-Szekeres)
- Penrose diagrams (conformal compactification)

### Technical Requirements
- Python 3.10+ (recommended) OR TypeScript/JavaScript
- NumPy, SciPy (for numerical integration)
- Matplotlib or WebGL-based renderer (for visualization)

## Project Structure

```
whitehole/
├── README.md
├── vision.md
├── docs/
│   ├── ARCHITECTURE.md
│   ├── GETTING_STARTED.md
│   └── ROADMAP.md
├── src/
│   ├── spacetimes/          # Layer 1: Metric definitions
│   │   ├── schwarzschild.py
│   │   ├── vaidya.py
│   │   └── base.py
│   ├── causal/              # Layer 2: Causal computation
│   │   ├── geodesic.py
│   │   ├── classifier.py
│   │   └── horizon.py
│   ├── visualization/       # Layer 3: Rendering
│   │   ├── penrose.py
│   │   ├── kruskal.py
│   │   └── interactive.py
│   └── utils/
│       ├── coordinates.py
│       └── integration.py
├── tests/
├── examples/
└── notebooks/               # Jupyter notebooks for exploration
```

## Phase 1: Core Implementation

### Step 1: Implement the Schwarzschild Metric

Start with the simplest non-trivial spacetime.

**File: `src/spacetimes/base.py`**
```python
from abc import ABC, abstractmethod
import numpy as np

class SpacetimeMetric(ABC):
    """Base class for all spacetime metrics."""

    @abstractmethod
    def metric(self, coords: np.ndarray) -> np.ndarray:
        """Return g_μν at given coordinates."""
        pass

    @abstractmethod
    def inverse_metric(self, coords: np.ndarray) -> np.ndarray:
        """Return g^μν at given coordinates."""
        pass

    @abstractmethod
    def christoffel(self, coords: np.ndarray) -> np.ndarray:
        """Return Γ^μ_νσ at given coordinates."""
        pass
```

**File: `src/spacetimes/schwarzschild.py`**
```python
import numpy as np
from .base import SpacetimeMetric

class Schwarzschild(SpacetimeMetric):
    """
    Schwarzschild metric in (t, r, θ, φ) coordinates.

    ds² = -(1 - 2M/r)dt² + (1 - 2M/r)⁻¹dr² + r²(dθ² + sin²θ dφ²)
    """

    def __init__(self, M: float = 1.0):
        self.M = M
        self.r_s = 2 * M  # Schwarzschild radius

    def metric(self, coords: np.ndarray) -> np.ndarray:
        """
        coords = [t, r, θ, φ]
        Returns 4x4 metric tensor g_μν
        """
        t, r, theta, phi = coords

        if r <= 0:
            raise ValueError("r must be positive")

        f = 1 - self.r_s / r

        g = np.zeros((4, 4))
        g[0, 0] = -f           # g_tt
        g[1, 1] = 1/f if f != 0 else np.inf  # g_rr
        g[2, 2] = r**2         # g_θθ
        g[3, 3] = r**2 * np.sin(theta)**2  # g_φφ

        return g

    def inverse_metric(self, coords: np.ndarray) -> np.ndarray:
        """Returns g^μν"""
        t, r, theta, phi = coords
        f = 1 - self.r_s / r

        g_inv = np.zeros((4, 4))
        g_inv[0, 0] = -1/f if f != 0 else -np.inf
        g_inv[1, 1] = f
        g_inv[2, 2] = 1/r**2
        g_inv[3, 3] = 1/(r**2 * np.sin(theta)**2)

        return g_inv

    def christoffel(self, coords: np.ndarray) -> np.ndarray:
        """
        Returns Christoffel symbols Γ^μ_νσ
        Index convention: Γ[μ, ν, σ]
        """
        t, r, theta, phi = coords
        M = self.M

        Gamma = np.zeros((4, 4, 4))

        f = 1 - 2*M/r

        # Non-zero components (symmetric in lower indices)
        Gamma[0, 0, 1] = Gamma[0, 1, 0] = M / (r**2 * f)
        Gamma[1, 0, 0] = M * f / r**2
        Gamma[1, 1, 1] = -M / (r**2 * f)
        Gamma[1, 2, 2] = -(r - 2*M)
        Gamma[1, 3, 3] = -(r - 2*M) * np.sin(theta)**2
        Gamma[2, 1, 2] = Gamma[2, 2, 1] = 1/r
        Gamma[2, 3, 3] = -np.sin(theta) * np.cos(theta)
        Gamma[3, 1, 3] = Gamma[3, 3, 1] = 1/r
        Gamma[3, 2, 3] = Gamma[3, 3, 2] = 1/np.tan(theta)

        return Gamma
```

### Step 2: Implement the Null Geodesic Engine

**File: `src/causal/geodesic.py`**
```python
import numpy as np
from scipy.integrate import solve_ivp

class NullGeodesicTracer:
    """
    Traces null geodesics through spacetime.

    Solves: d²xᵘ/dλ² + Γᵘ_νσ (dxᵛ/dλ)(dxᵠ/dλ) = 0
    """

    def __init__(self, spacetime):
        self.spacetime = spacetime

    def _geodesic_equation(self, lambda_param, state):
        """
        RHS of the geodesic equation in first-order form.

        state = [x⁰, x¹, x², x³, k⁰, k¹, k², k³]
        where k^μ = dx^μ/dλ
        """
        x = state[:4]  # Position
        k = state[4:]  # Tangent vector

        # Get Christoffel symbols at current position
        Gamma = self.spacetime.christoffel(x)

        # dx^μ/dλ = k^μ
        dx = k

        # dk^μ/dλ = -Γ^μ_νσ k^ν k^σ
        dk = np.zeros(4)
        for mu in range(4):
            for nu in range(4):
                for sigma in range(4):
                    dk[mu] -= Gamma[mu, nu, sigma] * k[nu] * k[sigma]

        return np.concatenate([dx, dk])

    def trace(self, initial_position, initial_direction,
              lambda_range=(0, 100), max_step=0.1):
        """
        Trace a null geodesic from initial conditions.

        Parameters:
            initial_position: [t, r, θ, φ]
            initial_direction: [k^t, k^r, k^θ, k^φ] (will be normalized to null)
            lambda_range: (λ_start, λ_end)

        Returns:
            solution object with .t (λ values) and .y (state vectors)
        """
        # Normalize direction to null
        k = self._normalize_to_null(initial_position, initial_direction)

        # Initial state
        state0 = np.concatenate([initial_position, k])

        # Stopping conditions
        def hit_singularity(lambda_param, state):
            return state[1] - 0.1  # r approaches 0
        hit_singularity.terminal = True

        def escape_to_infinity(lambda_param, state):
            return 1000 - state[1]  # r becomes very large
        escape_to_infinity.terminal = True

        # Integrate
        solution = solve_ivp(
            self._geodesic_equation,
            lambda_range,
            state0,
            method='RK45',
            events=[hit_singularity, escape_to_infinity],
            max_step=max_step,
            dense_output=True
        )

        return solution

    def _normalize_to_null(self, x, k):
        """Adjust k^μ so that g_μν k^μ k^ν = 0"""
        g = self.spacetime.metric(x)

        # Current norm
        norm_sq = 0
        for mu in range(4):
            for nu in range(4):
                norm_sq += g[mu, nu] * k[mu] * k[nu]

        # Solve for k^t such that null condition holds
        # g_tt (k^t)² + 2 g_ti k^t k^i + g_ij k^i k^j = 0
        # For diagonal metric: k^t = sqrt(-g_ij k^i k^j / g_tt)

        spatial_part = 0
        for i in range(1, 4):
            for j in range(1, 4):
                spatial_part += g[i, j] * k[i] * k[j]

        if g[0, 0] != 0:
            k_t_sq = -spatial_part / g[0, 0]
            if k_t_sq >= 0:
                k[0] = np.sqrt(k_t_sq) * np.sign(k[0]) if k[0] != 0 else np.sqrt(k_t_sq)

        return k
```

### Step 3: Implement Point Classification

**File: `src/causal/classifier.py`**
```python
import numpy as np
from .geodesic import NullGeodesicTracer

class CausalClassifier:
    """
    Classifies spacetime points as in J⁻(I⁺) or trapped.
    """

    def __init__(self, spacetime, angular_resolution=8):
        self.spacetime = spacetime
        self.tracer = NullGeodesicTracer(spacetime)
        self.angular_resolution = angular_resolution

    def classify_point(self, position, escape_radius=100):
        """
        Determine if a point can send signals to infinity.

        Returns:
            'escaping' if point ∈ J⁻(I⁺)
            'trapped' if point ∉ J⁻(I⁺)
        """
        directions = self._generate_outgoing_directions(position)

        for direction in directions:
            try:
                solution = self.tracer.trace(position, direction)

                # Check if any ray escaped
                final_r = solution.y[1, -1]
                if final_r > escape_radius:
                    return 'escaping'

            except Exception as e:
                continue  # Ray hit singularity or other issue

        return 'trapped'

    def _generate_outgoing_directions(self, position):
        """Generate null directions pointing outward."""
        directions = []

        # For spherical symmetry, vary the radial component sign
        # and angular directions
        for theta_dir in np.linspace(-1, 1, self.angular_resolution):
            for phi_dir in np.linspace(0, 2*np.pi, self.angular_resolution):
                # Outgoing radial (k^r > 0)
                k_r = 1.0
                k_theta = 0.1 * theta_dir
                k_phi = 0.1 * np.sin(phi_dir)
                k_t = 1.0  # Will be adjusted by null normalization

                directions.append(np.array([k_t, k_r, k_theta, k_phi]))

        return directions

    def classify_grid(self, r_range, t_range, theta=np.pi/2, phi=0):
        """
        Classify a 2D grid of points (useful for Penrose diagram).

        Returns:
            classification: 2D array of 'escaping'/'trapped'
            r_grid, t_grid: coordinate grids
        """
        r_vals = np.linspace(*r_range)
        t_vals = np.linspace(*t_range)

        classification = np.empty((len(t_vals), len(r_vals)), dtype=object)

        for i, t in enumerate(t_vals):
            for j, r in enumerate(r_vals):
                position = np.array([t, r, theta, phi])
                classification[i, j] = self.classify_point(position)

        r_grid, t_grid = np.meshgrid(r_vals, t_vals)
        return classification, r_grid, t_grid
```

### Step 4: Basic Penrose Visualization

**File: `src/visualization/penrose.py`**
```python
import numpy as np
import matplotlib.pyplot as plt

class PenroseVisualizer:
    """
    Visualize spacetime in Penrose (compactified) coordinates.
    """

    def __init__(self, spacetime):
        self.spacetime = spacetime

    def schwarzschild_to_kruskal(self, t, r, M=1.0):
        """Convert Schwarzschild (t, r) to Kruskal (U, V)."""
        r_s = 2 * M

        if r > r_s:
            factor = np.sqrt(r/r_s - 1) * np.exp(r / (2*r_s))
            U = factor * np.cosh(t / (2*r_s))
            V = factor * np.sinh(t / (2*r_s))
        else:
            factor = np.sqrt(1 - r/r_s) * np.exp(r / (2*r_s))
            U = factor * np.sinh(t / (2*r_s))
            V = factor * np.cosh(t / (2*r_s))

        return U, V

    def kruskal_to_penrose(self, U, V):
        """Compactify Kruskal to Penrose coordinates."""
        U_p = np.arctan(U)
        V_p = np.arctan(V)
        return U_p, V_p

    def plot_penrose_diagram(self, classification_data=None,
                              show_horizon=True, show_singularity=True):
        """
        Create a Penrose diagram visualization.
        """
        fig, ax = plt.subplots(figsize=(10, 10))

        # Draw the conformal boundary
        # I⁺ (future null infinity)
        ax.plot([0, np.pi/2], [np.pi/2, 0], 'b-', linewidth=2, label='I⁺')
        # I⁻ (past null infinity)
        ax.plot([0, np.pi/2], [-np.pi/2, 0], 'b-', linewidth=2, label='I⁻')
        # i⁰ (spatial infinity)
        ax.plot(np.pi/2, 0, 'bo', markersize=10, label='i⁰')
        # i⁺ (future timelike infinity)
        ax.plot(0, np.pi/2, 'b^', markersize=10, label='i⁺')
        # i⁻ (past timelike infinity)
        ax.plot(0, -np.pi/2, 'bv', markersize=10, label='i⁻')

        if show_horizon:
            # Event horizon (45° line from origin in Penrose coords)
            ax.plot([0, np.pi/4], [0, np.pi/4], 'r--', linewidth=2, label='Horizon')

        if show_singularity:
            # Singularity (r = 0, spacelike)
            ax.plot([0, np.pi/4], [np.pi/4, np.pi/4], 'k-', linewidth=3, label='Singularity')

        if classification_data is not None:
            # Overlay classification
            classification, r_grid, t_grid = classification_data
            # Convert to Penrose coordinates and plot
            # ... (implementation details)

        ax.set_xlim(-0.1, np.pi/2 + 0.1)
        ax.set_ylim(-np.pi/2 - 0.1, np.pi/2 + 0.1)
        ax.set_xlabel('Compactified radial coordinate')
        ax.set_ylabel('Compactified time coordinate')
        ax.set_title('Penrose Diagram')
        ax.legend()
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)

        return fig, ax
```

## Step 5: Your First Visualization

**File: `examples/first_visualization.py`**
```python
from src.spacetimes.schwarzschild import Schwarzschild
from src.causal.geodesic import NullGeodesicTracer
from src.visualization.penrose import PenroseVisualizer
import numpy as np
import matplotlib.pyplot as plt

# Create spacetime
spacetime = Schwarzschild(M=1.0)

# Create tracer
tracer = NullGeodesicTracer(spacetime)

# Launch a null ray from outside the horizon
# Position: t=0, r=5M, θ=π/2, φ=0
initial_pos = np.array([0.0, 5.0, np.pi/2, 0.0])
# Direction: outgoing radial
initial_dir = np.array([1.0, 1.0, 0.0, 0.0])

# Trace the geodesic
solution = tracer.trace(initial_pos, initial_dir)

# Plot in Schwarzschild coordinates
plt.figure(figsize=(10, 6))
plt.plot(solution.y[1], solution.y[0], 'b-', label='Outgoing ray from r=5M')

# Launch an ingoing ray
initial_dir_in = np.array([1.0, -1.0, 0.0, 0.0])
solution_in = tracer.trace(initial_pos, initial_dir_in)
plt.plot(solution_in.y[1], solution_in.y[0], 'r-', label='Ingoing ray from r=5M')

plt.axvline(x=2.0, color='k', linestyle='--', label='Horizon (r=2M)')
plt.xlabel('r/M')
plt.ylabel('t/M')
plt.title('Null Geodesics in Schwarzschild Spacetime')
plt.legend()
plt.grid(True)
plt.savefig('first_geodesics.png')
plt.show()

print("First visualization complete!")
```

## Running the Example

```bash
# From project root
python examples/first_visualization.py
```

You should see null geodesics traced through Schwarzschild spacetime, with outgoing rays escaping to infinity and ingoing rays falling toward the horizon.

## Next Steps

1. **Implement Vaidya metric** — Add dynamic horizon formation
2. **Build point classifier** — Determine J⁻(I⁺) membership
3. **Create interactive visualization** — Allow user to click and launch rays
4. **Add Penrose coordinate transformation** — Proper conformal diagram

See [ROADMAP.md](ROADMAP.md) for the complete development plan.

## Troubleshooting

### Geodesic Integration Fails Near Horizon
- Reduce `max_step` parameter
- Use adaptive stepping near r = 2M
- Consider tortoise coordinate transformation

### Numerical Instabilities
- Check Christoffel symbol computation
- Verify null normalization
- Ensure r > 0 always

### Coordinate Singularities
- Remember: Schwarzschild coordinates are singular at r = 2M
- Use Kruskal for near-horizon physics
- Implement proper coordinate patches
