"""
Schwarzschild spacetime metric.

The Schwarzschild solution describes the spacetime outside a spherically
symmetric, non-rotating mass M. It is the simplest black hole solution.

Metric (in Schwarzschild coordinates):
    ds² = -(1-2M/r)dt² + (1-2M/r)⁻¹dr² + r²(dθ² + sin²θ dφ²)

Coordinates: [t, r, θ, φ]
    - t: Schwarzschild time (coordinate time at infinity)
    - r: Areal radius
    - θ: Polar angle [0, π]
    - φ: Azimuthal angle [0, 2π]

Key features:
    - Event horizon at r = 2M (coordinate singularity)
    - Physical singularity at r = 0
    - Asymptotically flat as r → ∞
"""

import numpy as np
from numpy.typing import NDArray

from whitehole.spacetimes.base import SpacetimeMetric, CoordinateSystem


class Schwarzschild(SpacetimeMetric):
    """
    Schwarzschild black hole spacetime.

    Parameters:
        M: Black hole mass (in geometric units where G=c=1)
           Default M=1 sets the mass scale.

    Coordinate singularity:
        At r = 2M, the metric components diverge. This is a coordinate
        artifact - use Kruskal coordinates for physics at/across the horizon.

    Physical singularity:
        At r = 0, there is a genuine curvature singularity.
    """

    def __init__(self, M: float = 1.0):
        if M <= 0:
            raise ValueError("Mass M must be positive")
        self.M = M
        self.r_s = 2 * M  # Schwarzschild radius

    @property
    def name(self) -> str:
        return f"Schwarzschild (M={self.M})"

    @property
    def coordinate_system(self) -> CoordinateSystem:
        return CoordinateSystem.SCHWARZSCHILD

    @property
    def dimension(self) -> int:
        return 4

    def _f(self, r: float) -> float:
        """
        The metric function f(r) = 1 - 2M/r.

        This is the key function that determines the causal structure.
        f > 0: outside horizon (r > 2M)
        f = 0: at horizon (r = 2M)
        f < 0: inside horizon (r < 2M)
        """
        if r <= 0:
            return float('-inf')
        return 1.0 - self.r_s / r

    def is_valid_position(self, x: NDArray[np.float64]) -> bool:
        """
        Check if position is valid (not at singularity).

        Returns False for r ≤ 0 or very close to horizon.
        """
        r = x[1]
        if r <= 1e-10:  # At or past singularity
            return False
        if abs(r - self.r_s) < 1e-10:  # At coordinate singularity
            return False
        return True

    def metric(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Compute metric tensor g_μν.

        g_μν = diag(-(1-2M/r), (1-2M/r)⁻¹, r², r²sin²θ)
        """
        t, r, theta, phi = x

        if r <= 0:
            raise ValueError(f"Invalid radius r={r} ≤ 0")

        f = self._f(r)

        g = np.zeros((4, 4), dtype=np.float64)

        # g_tt
        g[0, 0] = -f

        # g_rr (handle horizon carefully)
        if abs(f) > 1e-14:
            g[1, 1] = 1.0 / f
        else:
            # At horizon: return large but finite value
            g[1, 1] = np.sign(f) * 1e14 if f != 0 else 1e14

        # g_θθ
        g[2, 2] = r * r

        # g_φφ
        sin_theta = np.sin(theta)
        g[3, 3] = r * r * sin_theta * sin_theta

        return g

    def inverse_metric(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Compute inverse metric g^μν.

        g^μν = diag(-(1-2M/r)⁻¹, (1-2M/r), r⁻², (r²sin²θ)⁻¹)
        """
        t, r, theta, phi = x

        if r <= 0:
            raise ValueError(f"Invalid radius r={r} ≤ 0")

        f = self._f(r)

        g_inv = np.zeros((4, 4), dtype=np.float64)

        # g^tt
        if abs(f) > 1e-14:
            g_inv[0, 0] = -1.0 / f
        else:
            g_inv[0, 0] = -np.sign(f) * 1e14 if f != 0 else -1e14

        # g^rr
        g_inv[1, 1] = f

        # g^θθ
        g_inv[2, 2] = 1.0 / (r * r)

        # g^φφ
        sin_theta = np.sin(theta)
        if abs(sin_theta) > 1e-14:
            g_inv[3, 3] = 1.0 / (r * r * sin_theta * sin_theta)
        else:
            g_inv[3, 3] = 1e14  # Near pole

        return g_inv

    def christoffel(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Compute Christoffel symbols Γ^μ_νσ.

        Non-zero components for Schwarzschild:
            Γ^t_tr = Γ^t_rt = M/(r²f)
            Γ^r_tt = Mf/r²
            Γ^r_rr = -M/(r²f)
            Γ^r_θθ = -(r-2M)
            Γ^r_φφ = -(r-2M)sin²θ
            Γ^θ_rθ = Γ^θ_θr = 1/r
            Γ^θ_φφ = -sinθ cosθ
            Γ^φ_rφ = Γ^φ_φr = 1/r
            Γ^φ_θφ = Γ^φ_φθ = cotθ

        Returns:
            4x4x4 array with Gamma[μ, ν, σ] = Γ^μ_νσ
        """
        t, r, theta, phi = x

        if r <= 0:
            raise ValueError(f"Invalid radius r={r} ≤ 0")

        M = self.M
        f = self._f(r)

        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)

        Gamma = np.zeros((4, 4, 4), dtype=np.float64)

        # Common factors
        r_sq = r * r

        # Handle near-horizon carefully
        if abs(f) > 1e-14:
            M_over_r2f = M / (r_sq * f)
        else:
            M_over_r2f = np.sign(f) * 1e14 if f != 0 else 1e14

        # Γ^t_tr = Γ^t_rt = M/(r²f)
        Gamma[0, 0, 1] = M_over_r2f
        Gamma[0, 1, 0] = M_over_r2f

        # Γ^r_tt = Mf/r²
        Gamma[1, 0, 0] = M * f / r_sq

        # Γ^r_rr = -M/(r²f)
        Gamma[1, 1, 1] = -M_over_r2f

        # Γ^r_θθ = -(r - 2M) = -r*f
        Gamma[1, 2, 2] = -r * f

        # Γ^r_φφ = -(r - 2M)sin²θ = -r*f*sin²θ
        Gamma[1, 3, 3] = -r * f * sin_theta * sin_theta

        # Γ^θ_rθ = Γ^θ_θr = 1/r
        Gamma[2, 1, 2] = 1.0 / r
        Gamma[2, 2, 1] = 1.0 / r

        # Γ^θ_φφ = -sinθ cosθ
        Gamma[2, 3, 3] = -sin_theta * cos_theta

        # Γ^φ_rφ = Γ^φ_φr = 1/r
        Gamma[3, 1, 3] = 1.0 / r
        Gamma[3, 3, 1] = 1.0 / r

        # Γ^φ_θφ = Γ^φ_φθ = cotθ = cosθ/sinθ
        if abs(sin_theta) > 1e-14:
            cot_theta = cos_theta / sin_theta
            Gamma[3, 2, 3] = cot_theta
            Gamma[3, 3, 2] = cot_theta

        return Gamma

    def is_outside_horizon(self, x: NDArray[np.float64]) -> bool:
        """Check if position is outside the event horizon."""
        return x[1] > self.r_s

    def is_inside_horizon(self, x: NDArray[np.float64]) -> bool:
        """Check if position is inside the event horizon."""
        return 0 < x[1] < self.r_s

    def proper_distance_to_horizon(self, r: float) -> float:
        """
        Compute proper radial distance from r to the horizon.

        For r > 2M:
            ∫_{2M}^{r} dr/√(1-2M/r)

        This integral has an analytic solution.
        """
        if r <= self.r_s:
            return 0.0

        M = self.M
        rs = self.r_s

        # Proper distance: r*√(1-2M/r) + 2M*ln(√(r/2M - 1) + √(r/2M))
        # Simplified for r > 2M
        sqrt_term = np.sqrt(r / rs - 1)
        return (r * np.sqrt(1 - rs/r) +
                rs * np.log(sqrt_term + np.sqrt(r / rs)))

    def surface_gravity(self) -> float:
        """
        Surface gravity κ at the horizon.

        κ = 1/(4M) for Schwarzschild

        This determines the Hawking temperature: T = κ/(2π)
        """
        return 1.0 / (4.0 * self.M)

    def hawking_temperature(self) -> float:
        """
        Hawking temperature T = κ/(2π) = 1/(8πM).

        In geometric units (G=c=ℏ=k_B=1).
        """
        return self.surface_gravity() / (2.0 * np.pi)

    def isco_radius(self) -> float:
        """
        Innermost stable circular orbit (ISCO) radius.

        For Schwarzschild: r_ISCO = 6M
        """
        return 6.0 * self.M

    def photon_sphere_radius(self) -> float:
        """
        Photon sphere radius (unstable circular photon orbits).

        For Schwarzschild: r_photon = 3M
        """
        return 3.0 * self.M
