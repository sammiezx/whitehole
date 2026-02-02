"""
Null geodesic tracer.

Integrates the geodesic equation to trace light rays through spacetime.
This is the core computational engine for determining causal structure.

The geodesic equation:
    d²x^μ/dλ² + Γ^μ_νσ (dx^ν/dλ)(dx^σ/dλ) = 0

In first-order form:
    dx^μ/dλ = k^μ
    dk^μ/dλ = -Γ^μ_νσ k^ν k^σ

For null geodesics: g_μν k^μ k^ν = 0
"""

from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.integrate._ivp.ivp import OdeResult

from whitehole.spacetimes.base import SpacetimeMetric


class TerminationReason(Enum):
    """Why a geodesic integration stopped."""
    MAX_PARAMETER = "max_affine_parameter"
    ESCAPED = "escaped_to_infinity"
    SINGULARITY = "hit_singularity"
    HORIZON_CROSSING = "crossed_horizon"
    INTEGRATION_ERROR = "integration_error"
    CUSTOM_EVENT = "custom_termination"


@dataclass
class GeodesicResult:
    """
    Result of geodesic integration.

    Attributes:
        positions: Array of shape (N, 4) with positions along geodesic
        momenta: Array of shape (N, 4) with 4-momenta along geodesic
        affine_params: Array of shape (N,) with affine parameter values
        termination: Why the integration stopped
        is_null: Whether the geodesic is null (should always be True for light)
        final_norm_sq: Final value of g_μν k^μ k^ν (should be ~0 for null)
    """
    positions: NDArray[np.float64]
    momenta: NDArray[np.float64]
    affine_params: NDArray[np.float64]
    termination: TerminationReason
    is_null: bool = True
    final_norm_sq: float = 0.0
    message: str = ""

    @property
    def n_points(self) -> int:
        return len(self.affine_params)

    @property
    def initial_position(self) -> NDArray[np.float64]:
        return self.positions[0]

    @property
    def final_position(self) -> NDArray[np.float64]:
        return self.positions[-1]

    @property
    def initial_r(self) -> float:
        return float(self.positions[0, 1])

    @property
    def final_r(self) -> float:
        return float(self.positions[-1, 1])

    def escaped(self, r_threshold: float = 100.0) -> bool:
        """Did the geodesic escape to large r?"""
        return self.final_r > r_threshold

    def captured(self, r_threshold: float = 0.1) -> bool:
        """Did the geodesic fall into the singularity or cross the horizon?"""
        if self.termination == TerminationReason.HORIZON_CROSSING:
            return True
        if self.termination == TerminationReason.SINGULARITY:
            return True
        return self.final_r < r_threshold


class NullGeodesicTracer:
    """
    Traces null geodesics through a given spacetime.

    This class integrates the geodesic equation for light rays, handling:
    - Null constraint enforcement
    - Singularity detection
    - Escape detection
    - Adaptive stepping near horizons

    Example:
        >>> spacetime = Schwarzschild(M=1.0)
        >>> tracer = NullGeodesicTracer(spacetime)
        >>> result = tracer.trace(
        ...     position=[0, 5, np.pi/2, 0],
        ...     direction=[1, 1, 0, 0]  # outgoing radial
        ... )
        >>> print(f"Escaped: {result.escaped()}")
    """

    def __init__(
        self,
        spacetime: SpacetimeMetric,
        rtol: float = 1e-10,
        atol: float = 1e-12,
        max_step: float = 0.5,
    ):
        """
        Initialize the geodesic tracer.

        Args:
            spacetime: The spacetime metric to trace through
            rtol: Relative tolerance for ODE solver
            atol: Absolute tolerance for ODE solver
            max_step: Maximum step size in affine parameter
        """
        self.spacetime = spacetime
        self.rtol = rtol
        self.atol = atol
        self.max_step = max_step

    def _geodesic_rhs(
        self,
        lambda_param: float,
        state: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """
        Right-hand side of geodesic equation.

        dy/dλ where y = [x^μ, k^μ]
        """
        return self.spacetime.geodesic_equation_rhs(state)

    def _make_singularity_event(
        self,
        r_min: float = 0.01
    ) -> Callable[[float, NDArray], float]:
        """Create event function that triggers at r → 0."""
        def event(lambda_param: float, state: NDArray) -> float:
            r = state[1]
            return r - r_min
        event.terminal = True  # type: ignore[attr-defined]
        event.direction = -1  # type: ignore[attr-defined]
        return event

    def _make_escape_event(
        self,
        r_max: float = 1000.0
    ) -> Callable[[float, NDArray], float]:
        """Create event function that triggers at large r."""
        def event(lambda_param: float, state: NDArray) -> float:
            r = state[1]
            return r_max - r
        event.terminal = True  # type: ignore[attr-defined]
        event.direction = -1  # type: ignore[attr-defined]
        return event

    def _make_horizon_event(
        self,
        r_horizon: float,
        terminal: bool = True
    ) -> Callable[[float, NDArray], float]:
        """Create event function that triggers at horizon crossing."""
        def event(lambda_param: float, state: NDArray) -> float:
            r = state[1]
            # Stop slightly before horizon to avoid coordinate singularity
            return r - (r_horizon + 0.001)
        event.terminal = terminal  # type: ignore[attr-defined]
        event.direction = -1  # type: ignore[attr-defined]
        return event

    def trace(
        self,
        position: NDArray[np.float64],
        direction: NDArray[np.float64],
        lambda_span: tuple[float, float] = (0.0, 1000.0),
        r_escape: float = 500.0,
        r_singularity: float = 0.01,
        r_horizon: float | None = None,
        dense_output: bool = False,
        enforce_null: bool = True,
    ) -> GeodesicResult:
        """
        Trace a null geodesic from initial conditions.

        Args:
            position: Initial 4-position [t, r, θ, φ]
            direction: Initial direction [k^t, k^r, k^θ, k^φ]
                      Will be normalized to null if enforce_null=True
            lambda_span: (λ_start, λ_end) for affine parameter
            r_escape: Consider escaped if r exceeds this
            r_singularity: Consider at singularity if r below this
            r_horizon: Horizon radius (stops integration to avoid coord singularity)
            dense_output: If True, return interpolated solution
            enforce_null: If True, adjust direction to satisfy null condition

        Returns:
            GeodesicResult with trajectory data
        """
        position = np.asarray(position, dtype=np.float64)
        direction = np.asarray(direction, dtype=np.float64)

        # Validate initial position
        if not self.spacetime.is_valid_position(position):
            return GeodesicResult(
                positions=position.reshape(1, 4),
                momenta=direction.reshape(1, 4),
                affine_params=np.array([0.0]),
                termination=TerminationReason.SINGULARITY,
                message="Invalid initial position"
            )

        # Normalize to null if requested
        if enforce_null:
            try:
                direction = self.spacetime.normalize_to_null(position, direction)
            except ValueError as e:
                return GeodesicResult(
                    positions=position.reshape(1, 4),
                    momenta=direction.reshape(1, 4),
                    affine_params=np.array([0.0]),
                    termination=TerminationReason.INTEGRATION_ERROR,
                    message=f"Cannot normalize to null: {e}"
                )

        # Initial state
        state0 = np.concatenate([position, direction])

        # Set up termination events
        events = [
            self._make_singularity_event(r_singularity),
            self._make_escape_event(r_escape),
        ]

        # Add horizon event to avoid coordinate singularity in Schwarzschild coords
        # Auto-detect horizon if spacetime has r_s attribute
        if r_horizon is None and hasattr(self.spacetime, 'r_s'):
            r_horizon = self.spacetime.r_s
        if r_horizon is not None:
            events.append(self._make_horizon_event(r_horizon, terminal=True))

        # Integrate
        try:
            solution = solve_ivp(
                self._geodesic_rhs,
                lambda_span,
                state0,
                method='DOP853',  # High-order method for accuracy
                events=events,
                dense_output=dense_output,
                rtol=self.rtol,
                atol=self.atol,
                max_step=self.max_step,
            )
        except Exception as e:
            return GeodesicResult(
                positions=position.reshape(1, 4),
                momenta=direction.reshape(1, 4),
                affine_params=np.array([0.0]),
                termination=TerminationReason.INTEGRATION_ERROR,
                message=str(e)
            )

        # Extract results
        positions = solution.y[:4, :].T
        momenta = solution.y[4:, :].T
        affine_params = solution.t

        # Determine termination reason
        termination = self._determine_termination(
            solution, r_escape, r_singularity, lambda_span
        )

        # Check final null condition
        final_pos = positions[-1]
        final_mom = momenta[-1]
        final_norm_sq = self.spacetime.norm_squared(final_pos, final_mom)

        return GeodesicResult(
            positions=positions,
            momenta=momenta,
            affine_params=affine_params,
            termination=termination,
            is_null=True,
            final_norm_sq=final_norm_sq,
            message=solution.message if hasattr(solution, 'message') else ""
        )

    def _determine_termination(
        self,
        solution: OdeResult,
        r_escape: float,
        r_singularity: float,
        lambda_span: tuple[float, float]
    ) -> TerminationReason:
        """Determine why integration stopped."""
        # Check events
        if solution.t_events is not None:
            # Event 0: singularity
            if len(solution.t_events[0]) > 0:
                return TerminationReason.SINGULARITY
            # Event 1: escape
            if len(solution.t_events[1]) > 0:
                return TerminationReason.ESCAPED
            # Event 2: horizon crossing (if present)
            if len(solution.t_events) > 2 and len(solution.t_events[2]) > 0:
                return TerminationReason.HORIZON_CROSSING

        # Check final position
        final_r = solution.y[1, -1]
        if final_r < r_singularity:
            return TerminationReason.SINGULARITY
        if final_r > r_escape:
            return TerminationReason.ESCAPED

        # Check if near horizon (coordinate singularity in Schwarzschild coords)
        if hasattr(self.spacetime, 'r_s'):
            if abs(final_r - self.spacetime.r_s) < 0.1:
                return TerminationReason.HORIZON_CROSSING

        # Check if we hit max parameter
        if abs(solution.t[-1] - lambda_span[1]) < 1e-6:
            return TerminationReason.MAX_PARAMETER

        return TerminationReason.INTEGRATION_ERROR

    def trace_radial(
        self,
        r_start: float,
        outgoing: bool = True,
        t_start: float = 0.0,
        theta: float = np.pi / 2,
        phi: float = 0.0,
        **kwargs
    ) -> GeodesicResult:
        """
        Trace a radial null geodesic (convenience method).

        Args:
            r_start: Initial radial coordinate
            outgoing: If True, trace outgoing ray; if False, ingoing
            t_start: Initial time coordinate
            theta, phi: Angular coordinates (default: equatorial plane)
            **kwargs: Additional arguments passed to trace()

        Returns:
            GeodesicResult
        """
        position = np.array([t_start, r_start, theta, phi])

        # Radial direction: k^r = ±1, k^θ = k^φ = 0
        k_r = 1.0 if outgoing else -1.0
        direction = np.array([1.0, k_r, 0.0, 0.0])

        return self.trace(position, direction, **kwargs)

    def trace_many(
        self,
        positions: list[NDArray[np.float64]],
        directions: list[NDArray[np.float64]],
        **kwargs
    ) -> list[GeodesicResult]:
        """
        Trace multiple geodesics (for later parallelization).

        Args:
            positions: List of initial positions
            directions: List of initial directions
            **kwargs: Arguments passed to trace()

        Returns:
            List of GeodesicResult objects
        """
        results = []
        for pos, dir in zip(positions, directions):
            results.append(self.trace(pos, dir, **kwargs))
        return results

    def generate_outgoing_directions(
        self,
        position: NDArray[np.float64],
        n_theta: int = 8,
        n_phi: int = 8,
        radial_component: float = 1.0
    ) -> list[NDArray[np.float64]]:
        """
        Generate a set of outgoing null directions.

        Used for testing whether a point can send signals to infinity.

        Args:
            position: The spacetime point
            n_theta: Number of polar angle samples
            n_phi: Number of azimuthal angle samples
            radial_component: Magnitude of k^r (positive = outgoing)

        Returns:
            List of direction 4-vectors (not yet null-normalized)
        """
        directions = []

        # Purely radial outgoing
        directions.append(np.array([1.0, radial_component, 0.0, 0.0]))

        # Add angular spread
        for theta_dir in np.linspace(-0.5, 0.5, n_theta):
            for phi_dir in np.linspace(0, 2*np.pi, n_phi, endpoint=False):
                k_theta = theta_dir * 0.1
                k_phi = np.sin(phi_dir) * 0.1

                directions.append(np.array([
                    1.0,
                    radial_component,
                    k_theta,
                    k_phi
                ]))

        return directions
