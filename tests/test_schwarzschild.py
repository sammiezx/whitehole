"""
Tests for Schwarzschild spacetime implementation.

Verifies:
    1. Metric properties (signature, determinant)
    2. Christoffel symbol symmetry
    3. Known analytic results (horizons, geodesics)
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from whitehole.spacetimes.schwarzschild import Schwarzschild


class TestSchwarzschildMetric:
    """Tests for the Schwarzschild metric implementation."""

    @pytest.fixture
    def spacetime(self):
        """Create a standard Schwarzschild spacetime with M=1."""
        return Schwarzschild(M=1.0)

    def test_metric_signature(self, spacetime):
        """Verify (-,+,+,+) signature."""
        # Test at r = 5M (outside horizon)
        x = np.array([0.0, 5.0, np.pi/2, 0.0])
        g = spacetime.metric(x)

        # Eigenvalues should have signs (-,+,+,+)
        eigenvalues = np.linalg.eigvalsh(g)
        eigenvalues.sort()

        assert eigenvalues[0] < 0, "Should have one negative eigenvalue"
        assert eigenvalues[1] > 0, "Should have three positive eigenvalues"
        assert eigenvalues[2] > 0
        assert eigenvalues[3] > 0

    def test_metric_diagonal(self, spacetime):
        """Schwarzschild metric should be diagonal."""
        x = np.array([0.0, 5.0, np.pi/2, 0.0])
        g = spacetime.metric(x)

        # Off-diagonal elements should be zero
        for i in range(4):
            for j in range(4):
                if i != j:
                    assert_allclose(g[i, j], 0.0, atol=1e-14)

    def test_metric_components_values(self, spacetime):
        """Verify metric components at known location."""
        M = 1.0
        r = 5.0
        theta = np.pi / 2

        x = np.array([0.0, r, theta, 0.0])
        g = spacetime.metric(x)

        # g_tt = -(1 - 2M/r) = -(1 - 0.4) = -0.6
        expected_gtt = -(1 - 2*M/r)
        assert_allclose(g[0, 0], expected_gtt, rtol=1e-10)

        # g_rr = 1/(1 - 2M/r) = 1/0.6 = 5/3
        expected_grr = 1 / (1 - 2*M/r)
        assert_allclose(g[1, 1], expected_grr, rtol=1e-10)

        # g_θθ = r²
        assert_allclose(g[2, 2], r**2, rtol=1e-10)

        # g_φφ = r² sin²θ (at θ = π/2, sin²θ = 1)
        assert_allclose(g[3, 3], r**2, rtol=1e-10)

    def test_inverse_metric_is_inverse(self, spacetime):
        """g^μρ g_ρν = δ^μ_ν."""
        x = np.array([0.0, 5.0, np.pi/2, 0.0])
        g = spacetime.metric(x)
        g_inv = spacetime.inverse_metric(x)

        # Product should be identity
        product = g_inv @ g
        identity = np.eye(4)

        assert_allclose(product, identity, atol=1e-10)

    def test_metric_at_different_radii(self, spacetime):
        """Test metric at various radii outside horizon."""
        for r in [3.0, 5.0, 10.0, 100.0]:
            x = np.array([0.0, r, np.pi/2, 0.0])
            g = spacetime.metric(x)
            g_inv = spacetime.inverse_metric(x)

            # Check g g^{-1} = I
            product = g @ g_inv
            assert_allclose(product, np.eye(4), atol=1e-10)

    def test_asymptotic_flatness(self, spacetime):
        """At large r, metric should approach Minkowski."""
        r = 1000.0
        x = np.array([0.0, r, np.pi/2, 0.0])
        g = spacetime.metric(x)

        # g_tt → -1
        assert_allclose(g[0, 0], -1.0, rtol=0.01)

        # g_rr → 1
        assert_allclose(g[1, 1], 1.0, rtol=0.01)


class TestChristoffelSymbols:
    """Tests for Christoffel symbol computation."""

    @pytest.fixture
    def spacetime(self):
        return Schwarzschild(M=1.0)

    def test_christoffel_symmetry(self, spacetime):
        """Γ^μ_νσ = Γ^μ_σν (symmetric in lower indices)."""
        x = np.array([0.0, 5.0, np.pi/2, 0.0])
        Gamma = spacetime.christoffel(x)

        for mu in range(4):
            for nu in range(4):
                for sigma in range(4):
                    assert_allclose(
                        Gamma[mu, nu, sigma],
                        Gamma[mu, sigma, nu],
                        atol=1e-14,
                        err_msg=f"Γ^{mu}_{nu}{sigma} ≠ Γ^{mu}_{sigma}{nu}"
                    )

    def test_known_christoffel_values(self, spacetime):
        """Verify specific Christoffel components."""
        M = 1.0
        r = 5.0
        theta = np.pi / 2
        f = 1 - 2*M/r  # = 0.6

        x = np.array([0.0, r, theta, 0.0])
        Gamma = spacetime.christoffel(x)

        # Γ^t_tr = M/(r²f) = 1/(25 * 0.6) = 1/15
        expected = M / (r**2 * f)
        assert_allclose(Gamma[0, 0, 1], expected, rtol=1e-10)
        assert_allclose(Gamma[0, 1, 0], expected, rtol=1e-10)

        # Γ^r_tt = Mf/r² = 0.6/25 = 0.024
        expected = M * f / r**2
        assert_allclose(Gamma[1, 0, 0], expected, rtol=1e-10)

        # Γ^r_rr = -M/(r²f) = -1/15
        expected = -M / (r**2 * f)
        assert_allclose(Gamma[1, 1, 1], expected, rtol=1e-10)

        # Γ^θ_rθ = 1/r = 0.2
        assert_allclose(Gamma[2, 1, 2], 1/r, rtol=1e-10)

        # Γ^φ_rφ = 1/r = 0.2
        assert_allclose(Gamma[3, 1, 3], 1/r, rtol=1e-10)


class TestSpacetimeProperties:
    """Tests for physical properties of the spacetime."""

    @pytest.fixture
    def spacetime(self):
        return Schwarzschild(M=1.0)

    def test_horizon_location(self, spacetime):
        """Horizon should be at r = 2M."""
        assert_allclose(spacetime.r_s, 2.0, rtol=1e-10)

    def test_is_outside_horizon(self, spacetime):
        """Test horizon detection."""
        x_outside = np.array([0.0, 3.0, np.pi/2, 0.0])
        x_at = np.array([0.0, 2.0, np.pi/2, 0.0])
        x_inside = np.array([0.0, 1.5, np.pi/2, 0.0])

        assert spacetime.is_outside_horizon(x_outside) is True
        assert spacetime.is_outside_horizon(x_at) is False
        assert spacetime.is_outside_horizon(x_inside) is False

    def test_isco_radius(self, spacetime):
        """ISCO should be at 6M."""
        assert_allclose(spacetime.isco_radius(), 6.0, rtol=1e-10)

    def test_photon_sphere(self, spacetime):
        """Photon sphere should be at 3M."""
        assert_allclose(spacetime.photon_sphere_radius(), 3.0, rtol=1e-10)

    def test_surface_gravity(self, spacetime):
        """Surface gravity κ = 1/(4M)."""
        expected = 1 / (4 * spacetime.M)
        assert_allclose(spacetime.surface_gravity(), expected, rtol=1e-10)

    def test_norm_squared_timelike(self, spacetime):
        """Timelike vector should have negative norm squared."""
        x = np.array([0.0, 5.0, np.pi/2, 0.0])

        # Stationary observer: u^μ = (1, 0, 0, 0) in coordinate basis
        # Needs normalization: g_tt u^t u^t = -1
        f = 1 - 2/5  # = 0.6
        u = np.array([1/np.sqrt(f), 0, 0, 0])

        norm_sq = spacetime.norm_squared(x, u)
        assert norm_sq < 0, "Timelike vector should have negative norm squared"
        assert_allclose(norm_sq, -1.0, rtol=1e-10)

    def test_norm_squared_null(self, spacetime):
        """Null vector should have zero norm squared."""
        x = np.array([0.0, 5.0, np.pi/2, 0.0])

        # Radial null vector (will be normalized)
        k = np.array([1.0, 1.0, 0.0, 0.0])
        k_null = spacetime.normalize_to_null(x, k)

        norm_sq = spacetime.norm_squared(x, k_null)
        assert_allclose(norm_sq, 0.0, atol=1e-12)


class TestNullNormalization:
    """Tests for null vector normalization."""

    @pytest.fixture
    def spacetime(self):
        return Schwarzschild(M=1.0)

    def test_normalize_outgoing_radial(self, spacetime):
        """Outgoing radial null ray should satisfy null condition."""
        x = np.array([0.0, 5.0, np.pi/2, 0.0])
        k = np.array([1.0, 1.0, 0.0, 0.0])

        k_null = spacetime.normalize_to_null(x, k)

        norm_sq = spacetime.norm_squared(x, k_null)
        assert_allclose(norm_sq, 0.0, atol=1e-12)

    def test_normalize_ingoing_radial(self, spacetime):
        """Ingoing radial null ray should satisfy null condition."""
        x = np.array([0.0, 5.0, np.pi/2, 0.0])
        k = np.array([1.0, -1.0, 0.0, 0.0])

        k_null = spacetime.normalize_to_null(x, k)

        norm_sq = spacetime.norm_squared(x, k_null)
        assert_allclose(norm_sq, 0.0, atol=1e-12)

    def test_normalize_preserves_spatial_direction(self, spacetime):
        """Normalization should only adjust k^t."""
        x = np.array([0.0, 5.0, np.pi/2, 0.0])
        k = np.array([1.0, 1.0, 0.1, 0.05])

        k_null = spacetime.normalize_to_null(x, k)

        # Spatial components should be unchanged
        assert_allclose(k_null[1], k[1], rtol=1e-10)
        assert_allclose(k_null[2], k[2], rtol=1e-10)
        assert_allclose(k_null[3], k[3], rtol=1e-10)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
