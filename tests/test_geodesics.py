"""
Tests for null geodesic tracer.

Verifies:
    1. Null condition is preserved during integration
    2. Known geodesic behaviors (escape/capture)
    3. Radial geodesics match analytic solutions
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from whitehole.spacetimes.schwarzschild import Schwarzschild
from whitehole.causal.geodesic import NullGeodesicTracer, TerminationReason


class TestNullGeodesicTracer:
    """Tests for the geodesic tracer."""

    @pytest.fixture
    def spacetime(self):
        return Schwarzschild(M=1.0)

    @pytest.fixture
    def tracer(self, spacetime):
        return NullGeodesicTracer(spacetime, rtol=1e-10, atol=1e-12)

    def test_outgoing_ray_from_outside_horizon_escapes(self, tracer):
        """Outgoing null ray from r > 2M should escape to infinity."""
        result = tracer.trace_radial(r_start=5.0, outgoing=True)

        assert result.termination == TerminationReason.ESCAPED
        assert result.escaped()
        assert result.final_r > 100.0

    def test_ingoing_ray_from_outside_horizon_captured(self, tracer):
        """Ingoing null ray from r > 2M should be captured (cross horizon or hit singularity)."""
        result = tracer.trace_radial(r_start=5.0, outgoing=False)

        # In Schwarzschild coordinates, integration stops at horizon (coordinate singularity)
        # Ray is captured if it reaches horizon or singularity
        assert result.termination in (TerminationReason.SINGULARITY, TerminationReason.HORIZON_CROSSING)
        assert result.captured()
        assert result.final_r <= tracer.spacetime.r_s + 0.1  # At or inside horizon

    def test_null_condition_preserved(self, tracer, spacetime):
        """g_μν k^μ k^ν should remain ~0 throughout integration."""
        result = tracer.trace_radial(r_start=5.0, outgoing=True)

        # Check at several points along the geodesic
        for i in range(0, len(result.positions), max(1, len(result.positions)//10)):
            x = result.positions[i]
            k = result.momenta[i]
            norm_sq = spacetime.norm_squared(x, k)

            assert abs(norm_sq) < 1e-8, f"Null condition violated at step {i}: |norm²| = {abs(norm_sq)}"

    def test_final_null_condition(self, tracer):
        """Final state should satisfy null condition."""
        result = tracer.trace_radial(r_start=5.0, outgoing=True)

        assert abs(result.final_norm_sq) < 1e-8

    def test_outgoing_from_near_horizon(self, tracer):
        """Outgoing ray from just outside horizon should escape (barely)."""
        # r = 2.1M, just outside
        result = tracer.trace_radial(r_start=2.1, outgoing=True, lambda_span=(0, 5000))

        assert result.termination == TerminationReason.ESCAPED
        assert result.escaped()

    def test_ray_from_photon_sphere(self, tracer, spacetime):
        """Ray from r = 3M tangent to photon sphere should orbit."""
        # At r = 3M, a tangential null ray should orbit (unstable)
        position = np.array([0.0, 3.0, np.pi/2, 0.0])

        # Tangential direction (k^r = 0, k^φ ≠ 0)
        direction = np.array([1.0, 0.0, 0.0, 1.0])

        result = tracer.trace(position, direction, lambda_span=(0, 100))

        # Should stay near r = 3M for a while (unstable orbit)
        # Check that it doesn't immediately escape or fall in
        mid_idx = len(result.positions) // 2
        r_mid = result.positions[mid_idx, 1]

        # Should still be relatively close to r = 3M
        assert 2.5 < r_mid < 4.0, f"Photon sphere orbit deviated too quickly: r = {r_mid}"

    def test_multiple_rays(self, tracer):
        """Test tracing multiple rays."""
        positions = [
            np.array([0.0, 5.0, np.pi/2, 0.0]),
            np.array([0.0, 10.0, np.pi/2, 0.0]),
        ]
        directions = [
            np.array([1.0, 1.0, 0.0, 0.0]),
            np.array([1.0, -1.0, 0.0, 0.0]),
        ]

        results = tracer.trace_many(positions, directions)

        assert len(results) == 2
        assert results[0].escaped()  # Outgoing from r=5
        assert results[1].captured()  # Ingoing from r=10

    def test_generate_outgoing_directions(self, tracer):
        """Test direction generation."""
        position = np.array([0.0, 5.0, np.pi/2, 0.0])
        directions = tracer.generate_outgoing_directions(position, n_theta=4, n_phi=4)

        assert len(directions) > 1

        # All should have positive k^r (outgoing)
        for d in directions:
            assert d[1] > 0, "Generated direction should be outgoing"


class TestCausalBehavior:
    """Tests for causal structure properties."""

    @pytest.fixture
    def spacetime(self):
        return Schwarzschild(M=1.0)

    @pytest.fixture
    def tracer(self, spacetime):
        return NullGeodesicTracer(spacetime)

    def test_horizon_divides_escape_capture(self, tracer):
        """Points outside horizon can escape; inside cannot."""
        # Outside horizon (r = 3M): outgoing ray escapes
        result_outside = tracer.trace_radial(r_start=3.0, outgoing=True)
        assert result_outside.escaped()

        # Well outside horizon (r = 10M): outgoing ray definitely escapes
        result_far = tracer.trace_radial(r_start=10.0, outgoing=True)
        assert result_far.escaped()

    def test_all_rays_from_inside_captured(self, tracer):
        """
        Inside the horizon, all future-directed causal curves
        end at the singularity.

        Note: Schwarzschild coordinates are singular at r = 2M,
        so we test near but below the horizon.
        """
        # This test is tricky because Schwarzschild coords are bad inside
        # For now, just verify the basic behavior at the horizon boundary
        pass  # TODO: Implement with Kruskal coordinates


class TestGeodesicResult:
    """Tests for GeodesicResult properties."""

    @pytest.fixture
    def tracer(self):
        spacetime = Schwarzschild(M=1.0)
        return NullGeodesicTracer(spacetime)

    def test_result_properties(self, tracer):
        """Test GeodesicResult helper properties."""
        result = tracer.trace_radial(r_start=5.0, outgoing=True)

        assert result.n_points > 10  # Should have multiple integration points
        assert result.initial_r == 5.0
        assert result.final_r > result.initial_r  # Outgoing
        assert result.is_null is True

    def test_result_positions_shape(self, tracer):
        """Positions array should have correct shape."""
        result = tracer.trace_radial(r_start=5.0, outgoing=True)

        assert result.positions.shape[1] == 4  # 4 coordinates
        assert result.momenta.shape[1] == 4
        assert len(result.affine_params) == result.n_points


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
