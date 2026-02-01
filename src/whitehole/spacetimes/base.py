"""
Base class for spacetime metrics.

All spacetime implementations must provide:
- Metric tensor g_μν
- Inverse metric g^μν
- Christoffel symbols Γ^μ_νσ

Coordinates are always 4-vectors: [x⁰, x¹, x², x³]
For spherically symmetric spacetimes: [t, r, θ, φ]
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum
from typing import Callable

import numpy as np
from numpy.typing import NDArray


class CoordinateSystem(Enum):
    """Supported coordinate systems."""
    SCHWARZSCHILD = "schwarzschild"      # (t, r, θ, φ)
    EDDINGTON_FINKELSTEIN = "eddington"  # (v, r, θ, φ) ingoing
    KRUSKAL = "kruskal"                  # (U, V, θ, φ)
    PAINLEVE_GULLSTRAND = "painleve"     # (T, r, θ, φ)


@dataclass
class GeodesicState:
    """
    State of a geodesic at a given affine parameter.

    Attributes:
        position: 4-position x^μ
        momentum: 4-momentum k^μ = dx^μ/dλ
        affine_param: Current value of affine parameter λ
        is_null: Whether this is a null geodesic
    """
    position: NDArray[np.float64]
    momentum: NDArray[np.float64]
    affine_param: float = 0.0
    is_null: bool = True

    def __post_init__(self):
        self.position = np.asarray(self.position, dtype=np.float64)
        self.momentum = np.asarray(self.momentum, dtype=np.float64)


class SpacetimeMetric(ABC):
    """
    Abstract base class for spacetime metrics.

    All implementations must provide the metric tensor, its inverse,
    and Christoffel symbols. These are the minimum requirements for
    geodesic integration.

    Convention:
        - Signature: (-,+,+,+)
        - Index ordering: μ, ν, σ ∈ {0, 1, 2, 3}
        - Christoffel symbols: Γ^μ_νσ stored as [μ, ν, σ]
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Human-readable name of the spacetime."""
        pass

    @property
    @abstractmethod
    def coordinate_system(self) -> CoordinateSystem:
        """The coordinate system used by this metric."""
        pass

    @property
    @abstractmethod
    def dimension(self) -> int:
        """Spacetime dimension (usually 4)."""
        pass

    @abstractmethod
    def metric(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Compute the metric tensor g_μν at position x.

        Args:
            x: 4-position [x⁰, x¹, x², x³]

        Returns:
            4x4 array representing g_μν
        """
        pass

    @abstractmethod
    def inverse_metric(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Compute the inverse metric g^μν at position x.

        Args:
            x: 4-position

        Returns:
            4x4 array representing g^μν
        """
        pass

    @abstractmethod
    def christoffel(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Compute Christoffel symbols Γ^μ_νσ at position x.

        Args:
            x: 4-position

        Returns:
            4x4x4 array where result[μ, ν, σ] = Γ^μ_νσ
        """
        pass

    def is_valid_position(self, x: NDArray[np.float64]) -> bool:
        """
        Check if a position is valid (not at a singularity, etc).

        Override in subclasses for specific singularity handling.
        """
        return True

    def norm_squared(self, x: NDArray[np.float64], v: NDArray[np.float64]) -> float:
        """
        Compute the squared norm of a vector v at position x.

        g_μν v^μ v^ν

        Returns:
            < 0 for timelike
            = 0 for null
            > 0 for spacelike
        """
        g = self.metric(x)
        return float(np.einsum('i,ij,j', v, g, v))

    def normalize_to_null(
        self,
        x: NDArray[np.float64],
        k: NDArray[np.float64],
        preserve_direction: bool = True
    ) -> NDArray[np.float64]:
        """
        Adjust a vector k to be null at position x.

        For a diagonal metric, solves for k⁰ such that g_μν k^μ k^ν = 0.

        Args:
            x: Position
            k: Initial 4-vector (k⁰ will be adjusted)
            preserve_direction: If True, keep sign of k⁰

        Returns:
            Null-normalized vector
        """
        k = np.asarray(k, dtype=np.float64).copy()
        g = self.metric(x)

        # For diagonal metric: g_00 (k^0)² + g_ii (k^i)² = 0
        # So: k^0 = sqrt(-g_ii k^i k^i / g_00)

        spatial_norm_sq = sum(g[i, i] * k[i]**2 for i in range(1, 4))

        # Add off-diagonal terms if present
        for i in range(1, 4):
            for j in range(i+1, 4):
                spatial_norm_sq += 2 * g[i, j] * k[i] * k[j]

        if g[0, 0] == 0:
            raise ValueError("g_00 = 0, cannot normalize")

        k0_sq = -spatial_norm_sq / g[0, 0]

        if k0_sq < 0:
            raise ValueError(f"Cannot make null: requires imaginary k⁰ (k0_sq={k0_sq})")

        k0_new = np.sqrt(k0_sq)

        if preserve_direction and k[0] < 0:
            k0_new = -k0_new

        k[0] = k0_new
        return k

    def project_to_null(
        self,
        x: NDArray[np.float64],
        k: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """
        Project a vector onto the null cone by scaling.

        Alternative to normalize_to_null: scales the entire vector
        rather than adjusting only k⁰.
        """
        norm_sq = self.norm_squared(x, k)
        if abs(norm_sq) < 1e-14:
            return k  # Already null

        # For small deviations, adjust k⁰
        return self.normalize_to_null(x, k)

    def determinant(self, x: NDArray[np.float64]) -> float:
        """Compute det(g_μν) at position x."""
        return float(np.linalg.det(self.metric(x)))

    def geodesic_equation_rhs(
        self,
        state: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """
        Right-hand side of the geodesic equation in first-order form.

        For state = [x⁰, x¹, x², x³, k⁰, k¹, k², k³]:
            dx^μ/dλ = k^μ
            dk^μ/dλ = -Γ^μ_νσ k^ν k^σ

        Args:
            state: 8-component state vector [x, k]

        Returns:
            8-component derivative [dx/dλ, dk/dλ]
        """
        x = state[:4]
        k = state[4:]

        # Position derivative
        dx = k.copy()

        # Momentum derivative from geodesic equation
        Gamma = self.christoffel(x)
        dk = -np.einsum('ijk,j,k->i', Gamma, k, k)

        return np.concatenate([dx, dk])
