"""
Causal point classifier.

Determines whether spacetime points belong to J⁻(I⁺), the causal past
of future null infinity. This is the key to finding event horizons:

    Event Horizon = ∂J⁻(I⁺)

Algorithm:
    For each point p, launch outgoing null rays.
    If ANY ray escapes to large r at late times: p ∈ J⁻(I⁺)
    If ALL rays are trapped (fall to singularity): p ∉ J⁻(I⁺)
"""

from dataclasses import dataclass
from enum import Enum
from typing import List, Tuple, Optional
import numpy as np
from numpy.typing import NDArray

from whitehole.spacetimes.base import SpacetimeMetric
from whitehole.causal.geodesic import NullGeodesicTracer, GeodesicResult, TerminationReason


class PointClassification(Enum):
    """Classification of a spacetime point's causal nature."""
    ESCAPING = "escaping"    # Point is in J⁻(I⁺) - can send signals to infinity
    TRAPPED = "trapped"      # Point is NOT in J⁻(I⁺) - all signals are trapped
    BOUNDARY = "boundary"    # Point is near the horizon (uncertain)
    INVALID = "invalid"      # Point is at/past singularity


@dataclass
class ClassificationResult:
    """
    Result of classifying a single spacetime point.

    Attributes:
        position: The classified point [t, r, θ, φ]
        classification: The causal classification
        escape_fraction: Fraction of rays that escaped (0 to 1)
        n_rays_traced: Total number of rays traced
        escaping_rays: Results of rays that escaped
        trapped_rays: Results of rays that were trapped
    """
    position: NDArray[np.float64]
    classification: PointClassification
    escape_fraction: float
    n_rays_traced: int
    escaping_rays: List[GeodesicResult]
    trapped_rays: List[GeodesicResult]

    @property
    def r(self) -> float:
        return float(self.position[1])

    @property
    def t(self) -> float:
        return float(self.position[0])


class CausalClassifier:
    """
    Classifies spacetime points as escaping (in J⁻(I⁺)) or trapped.

    This is the core algorithm for horizon detection. For Schwarzschild,
    points with r > 2M should be escaping, r < 2M should be trapped.

    Example:
        >>> spacetime = Schwarzschild(M=1.0)
        >>> classifier = CausalClassifier(spacetime)
        >>> result = classifier.classify_point([0, 3, np.pi/2, 0])
        >>> print(result.classification)  # Should be ESCAPING
    """

    def __init__(
        self,
        spacetime: SpacetimeMetric,
        tracer: Optional[NullGeodesicTracer] = None,
        n_theta: int = 4,
        n_phi: int = 4,
        escape_radius: float = 200.0,
        escape_threshold: float = 0.01,  # Fraction of rays needed to classify as escaping
    ):
        """
        Initialize the classifier.

        Args:
            spacetime: The spacetime metric
            tracer: Geodesic tracer (created if not provided)
            n_theta: Angular resolution in theta
            n_phi: Angular resolution in phi
            escape_radius: Radius beyond which we consider "escaped"
            escape_threshold: Minimum fraction of escaping rays for ESCAPING class
        """
        self.spacetime = spacetime
        self.tracer = tracer or NullGeodesicTracer(spacetime)
        self.n_theta = n_theta
        self.n_phi = n_phi
        self.escape_radius = escape_radius
        self.escape_threshold = escape_threshold

    def classify_point(
        self,
        position: NDArray[np.float64],
        store_rays: bool = False
    ) -> ClassificationResult:
        """
        Classify a single spacetime point.

        Args:
            position: 4-position [t, r, θ, φ]
            store_rays: If True, store individual ray results

        Returns:
            ClassificationResult with classification and metadata
        """
        position = np.asarray(position, dtype=np.float64)

        # Check validity
        if not self.spacetime.is_valid_position(position):
            return ClassificationResult(
                position=position,
                classification=PointClassification.INVALID,
                escape_fraction=0.0,
                n_rays_traced=0,
                escaping_rays=[],
                trapped_rays=[]
            )

        # Generate outgoing directions
        directions = self.tracer.generate_outgoing_directions(
            position, self.n_theta, self.n_phi
        )

        escaping_rays = []
        trapped_rays = []
        n_escaped = 0

        for direction in directions:
            result = self.tracer.trace(
                position, direction,
                r_escape=self.escape_radius,
            )

            if result.termination == TerminationReason.ESCAPED:
                n_escaped += 1
                if store_rays:
                    escaping_rays.append(result)
            else:
                if store_rays:
                    trapped_rays.append(result)

        n_total = len(directions)
        escape_fraction = n_escaped / n_total if n_total > 0 else 0.0

        # Determine classification
        if escape_fraction >= self.escape_threshold:
            classification = PointClassification.ESCAPING
        elif escape_fraction > 0:
            classification = PointClassification.BOUNDARY
        else:
            classification = PointClassification.TRAPPED

        return ClassificationResult(
            position=position,
            classification=classification,
            escape_fraction=escape_fraction,
            n_rays_traced=n_total,
            escaping_rays=escaping_rays if store_rays else [],
            trapped_rays=trapped_rays if store_rays else []
        )

    def classify_radial_slice(
        self,
        r_values: NDArray[np.float64],
        t: float = 0.0,
        theta: float = np.pi / 2,
        phi: float = 0.0,
    ) -> List[ClassificationResult]:
        """
        Classify points along a radial slice (fixed t, θ, φ).

        Useful for finding the horizon location.

        Args:
            r_values: Array of r coordinates to classify
            t: Time coordinate
            theta: Polar angle
            phi: Azimuthal angle

        Returns:
            List of ClassificationResult for each r value
        """
        results = []
        for r in r_values:
            position = np.array([t, r, theta, phi])
            results.append(self.classify_point(position))
        return results

    def find_horizon_radius(
        self,
        r_min: float,
        r_max: float,
        t: float = 0.0,
        theta: float = np.pi / 2,
        phi: float = 0.0,
        tolerance: float = 0.01,
        max_iterations: int = 50
    ) -> Tuple[float, float]:
        """
        Find the horizon radius using bisection.

        Searches for the boundary between escaping and trapped regions.

        Args:
            r_min: Lower bound for search (should be trapped)
            r_max: Upper bound for search (should be escaping)
            t: Time coordinate
            theta, phi: Angular coordinates
            tolerance: Desired precision in r
            max_iterations: Maximum bisection iterations

        Returns:
            (r_horizon, uncertainty) tuple
        """
        # Verify bounds
        pos_min = np.array([t, r_min, theta, phi])
        pos_max = np.array([t, r_max, theta, phi])

        result_min = self.classify_point(pos_min)
        result_max = self.classify_point(pos_max)

        if result_min.classification == PointClassification.ESCAPING:
            raise ValueError(f"r_min={r_min} is escaping, should be trapped")
        if result_max.classification == PointClassification.TRAPPED:
            raise ValueError(f"r_max={r_max} is trapped, should be escaping")

        # Bisection
        low, high = r_min, r_max
        for _ in range(max_iterations):
            mid = (low + high) / 2
            if high - low < tolerance:
                break

            pos_mid = np.array([t, mid, theta, phi])
            result_mid = self.classify_point(pos_mid)

            if result_mid.classification == PointClassification.ESCAPING:
                high = mid
            else:
                low = mid

        horizon_r = (low + high) / 2
        uncertainty = (high - low) / 2

        return horizon_r, uncertainty

    def classify_grid(
        self,
        t_values: NDArray[np.float64],
        r_values: NDArray[np.float64],
        theta: float = np.pi / 2,
        phi: float = 0.0,
    ) -> NDArray[np.int32]:
        """
        Classify a 2D grid of (t, r) points.

        Returns an integer array for visualization:
            1 = ESCAPING
            0 = BOUNDARY
           -1 = TRAPPED
           -2 = INVALID

        Args:
            t_values: Array of t coordinates
            r_values: Array of r coordinates
            theta, phi: Fixed angular coordinates

        Returns:
            2D array of shape (len(t_values), len(r_values))
        """
        result = np.zeros((len(t_values), len(r_values)), dtype=np.int32)

        classification_map = {
            PointClassification.ESCAPING: 1,
            PointClassification.BOUNDARY: 0,
            PointClassification.TRAPPED: -1,
            PointClassification.INVALID: -2,
        }

        for i, t in enumerate(t_values):
            for j, r in enumerate(r_values):
                position = np.array([t, r, theta, phi])
                class_result = self.classify_point(position)
                result[i, j] = classification_map[class_result.classification]

        return result
