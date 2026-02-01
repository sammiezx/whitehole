"""
Vaidya spacetime metric.

The Vaidya solution describes a spherically symmetric spacetime with
a radially moving null dust (radiation). It models:
    - Gravitational collapse with radiation
    - Evaporating black holes
    - Horizon formation

Ingoing Vaidya metric (Eddington-Finkelstein coordinates):
    ds² = -(1 - 2M(v)/r) dv² + 2 dv dr + r² dΩ²

where:
    - v is the advanced (ingoing) null coordinate
    - M(v) is the mass function (can increase or decrease)
    - For M(v) = const, reduces to Schwarzschild in EF coordinates

Key physics:
    - If dM/dv > 0: accreting/collapsing (horizon grows)
    - If dM/dv < 0: evaporating (horizon shrinks)
    - If dM/dv = 0: stationary (Schwarzschild limit)

The event horizon location depends on the ENTIRE future evolution M(v),
making Vaidya ideal for demonstrating horizon teleology.
"""

from typing import Callable, Optional
import numpy as np
from numpy.typing import NDArray

from whitehole.spacetimes.base import SpacetimeMetric, CoordinateSystem


class Vaidya(SpacetimeMetric):
    """
    Vaidya spacetime in ingoing Eddington-Finkelstein coordinates.

    The mass function M(v) determines the dynamics. Common choices:
        - M(v) = M_0 : Schwarzschild
        - M(v) = M_0 * v/v_0 for v ∈ [0, v_0] : Linear collapse
        - M(v) = M_0 * (1 - exp(-v/τ)) : Exponential approach to M_0
        - M(v) = M_0 * (1 - v/v_evap) : Linear evaporation

    Coordinates: [v, r, θ, φ]
        - v: Advanced (ingoing) Eddington-Finkelstein time
        - r: Areal radius
        - θ, φ: Angular coordinates

    Example:
        >>> # Collapsing shell reaching M=1 at v=10
        >>> mass_func = lambda v: min(v/10, 1.0)
        >>> spacetime = Vaidya(mass_function=mass_func)
    """

    def __init__(
        self,
        mass_function: Optional[Callable[[float], float]] = None,
        M_final: float = 1.0,
        collapse_time: float = 10.0
    ):
        """
        Initialize Vaidya spacetime.

        Args:
            mass_function: Function M(v) returning mass at advanced time v.
                          If None, uses linear collapse to M_final.
            M_final: Final mass (used if mass_function is None)
            collapse_time: Time to reach M_final (used if mass_function is None)
        """
        if mass_function is not None:
            self._mass_function = mass_function
        else:
            # Default: linear collapse
            self._mass_function = lambda v: min(M_final * max(v, 0) / collapse_time, M_final)

        self.M_final = M_final
        self.collapse_time = collapse_time

    @property
    def name(self) -> str:
        return "Vaidya (Ingoing)"

    @property
    def coordinate_system(self) -> CoordinateSystem:
        return CoordinateSystem.EDDINGTON_FINKELSTEIN

    @property
    def dimension(self) -> int:
        return 4

    def mass(self, v: float) -> float:
        """Evaluate mass function at advanced time v."""
        return self._mass_function(v)

    def mass_derivative(self, v: float, dv: float = 1e-6) -> float:
        """
        Numerical derivative dM/dv.

        Positive: accretion/collapse
        Negative: evaporation
        Zero: stationary
        """
        return (self.mass(v + dv) - self.mass(v - dv)) / (2 * dv)

    def _f(self, v: float, r: float) -> float:
        """Metric function f(v, r) = 1 - 2M(v)/r."""
        if r <= 0:
            return float('-inf')
        M = self.mass(v)
        return 1.0 - 2 * M / r

    def apparent_horizon_radius(self, v: float) -> float:
        """
        Apparent horizon radius at time v.

        For Vaidya: r_AH = 2M(v)

        Note: This is NOT the event horizon in general!
        The event horizon depends on the future.
        """
        return 2 * self.mass(v)

    def is_valid_position(self, x: NDArray[np.float64]) -> bool:
        """Check if position is valid (not at singularity)."""
        r = x[1]
        return r > 1e-10

    def metric(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Compute metric tensor g_μν in (v, r, θ, φ) coordinates.

        ds² = -f dv² + 2 dv dr + r² dΩ²

        g_μν = | -f   1   0   0     |
               |  1   0   0   0     |
               |  0   0   r²  0     |
               |  0   0   0   r²s²θ |
        """
        v, r, theta, phi = x

        if r <= 0:
            raise ValueError(f"Invalid radius r={r} ≤ 0")

        f = self._f(v, r)

        g = np.zeros((4, 4), dtype=np.float64)

        # g_vv
        g[0, 0] = -f

        # g_vr = g_rv = 1 (cross term)
        g[0, 1] = 1.0
        g[1, 0] = 1.0

        # g_rr = 0 (!)
        g[1, 1] = 0.0

        # Angular parts
        g[2, 2] = r * r
        g[3, 3] = r * r * np.sin(theta)**2

        return g

    def inverse_metric(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Compute inverse metric g^μν.

        For Vaidya in ingoing EF:
        g^μν = |  0   1   0   0       |
               |  1   f   0   0       |
               |  0   0  1/r² 0       |
               |  0   0   0  1/(r²s²) |
        """
        v, r, theta, phi = x

        if r <= 0:
            raise ValueError(f"Invalid radius r={r} ≤ 0")

        f = self._f(v, r)

        g_inv = np.zeros((4, 4), dtype=np.float64)

        # g^vv = 0
        g_inv[0, 0] = 0.0

        # g^vr = g^rv = 1
        g_inv[0, 1] = 1.0
        g_inv[1, 0] = 1.0

        # g^rr = f
        g_inv[1, 1] = f

        # Angular parts
        g_inv[2, 2] = 1.0 / (r * r)
        sin_theta = np.sin(theta)
        if abs(sin_theta) > 1e-14:
            g_inv[3, 3] = 1.0 / (r * r * sin_theta**2)
        else:
            g_inv[3, 3] = 1e14

        return g_inv

    def christoffel(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Compute Christoffel symbols Γ^μ_νσ for Vaidya.

        Note: Vaidya has additional terms due to ∂M/∂v ≠ 0.

        The non-zero components involve:
            - Standard Schwarzschild-like terms (with M → M(v))
            - Additional terms from dM/dv
        """
        v, r, theta, phi = x

        if r <= 0:
            raise ValueError(f"Invalid radius r={r} ≤ 0")

        M = self.mass(v)
        dM_dv = self.mass_derivative(v)
        f = self._f(v, r)

        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)

        Gamma = np.zeros((4, 4, 4), dtype=np.float64)

        r_sq = r * r

        # Γ^v_vv = dM/dv / r  (from mass evolution)
        Gamma[0, 0, 0] = dM_dv / r if r > 0 else 0

        # Γ^v_vr = Γ^v_rv = M/r²
        Gamma[0, 0, 1] = M / r_sq
        Gamma[0, 1, 0] = M / r_sq

        # Γ^r_vv = -f * M/r² + dM/dv * f / r  (complex!)
        Gamma[1, 0, 0] = f * M / r_sq - dM_dv / r

        # Γ^r_vr = Γ^r_rv = M/r²
        Gamma[1, 0, 1] = -M / r_sq
        Gamma[1, 1, 0] = -M / r_sq

        # Γ^r_θθ = -r * f (same as Schwarzschild with M(v))
        Gamma[1, 2, 2] = -r * f

        # Γ^r_φφ = -r * f * sin²θ
        Gamma[1, 3, 3] = -r * f * sin_theta * sin_theta

        # Γ^θ_vθ = Γ^θ_θv = 0 for ingoing Vaidya

        # Γ^θ_rθ = Γ^θ_θr = 1/r
        Gamma[2, 1, 2] = 1.0 / r
        Gamma[2, 2, 1] = 1.0 / r

        # Γ^θ_φφ = -sinθ cosθ
        Gamma[2, 3, 3] = -sin_theta * cos_theta

        # Γ^φ_rφ = Γ^φ_φr = 1/r
        Gamma[3, 1, 3] = 1.0 / r
        Gamma[3, 3, 1] = 1.0 / r

        # Γ^φ_θφ = Γ^φ_φθ = cotθ
        if abs(sin_theta) > 1e-14:
            cot_theta = cos_theta / sin_theta
            Gamma[3, 2, 3] = cot_theta
            Gamma[3, 3, 2] = cot_theta

        return Gamma


# Convenience factory functions for common mass profiles

def linear_collapse(M_final: float = 1.0, t_collapse: float = 10.0) -> Vaidya:
    """
    Create Vaidya spacetime with linear mass increase.

    M(v) = M_final * v / t_collapse  for 0 ≤ v ≤ t_collapse
    M(v) = M_final                   for v > t_collapse
    """
    def mass_func(v):
        if v <= 0:
            return 0.0
        elif v >= t_collapse:
            return M_final
        else:
            return M_final * v / t_collapse

    return Vaidya(mass_function=mass_func, M_final=M_final, collapse_time=t_collapse)


def step_collapse(M_final: float = 1.0, t_step: float = 5.0) -> Vaidya:
    """
    Create Vaidya spacetime with sudden (step) mass increase.

    Models a thin shell of matter falling in at v = t_step.
    """
    def mass_func(v):
        return M_final if v >= t_step else 0.0

    return Vaidya(mass_function=mass_func, M_final=M_final, collapse_time=t_step)


def exponential_collapse(M_final: float = 1.0, tau: float = 5.0) -> Vaidya:
    """
    Create Vaidya spacetime with exponential approach to M_final.

    M(v) = M_final * (1 - exp(-v/τ))
    """
    def mass_func(v):
        if v <= 0:
            return 0.0
        return M_final * (1 - np.exp(-v / tau))

    return Vaidya(mass_function=mass_func, M_final=M_final, collapse_time=3*tau)
