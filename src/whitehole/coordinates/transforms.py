"""
Coordinate transformations for Schwarzschild spacetime.

Key coordinate systems:
    1. Schwarzschild (t, r): Standard textbook coordinates
       - Singular at horizon r = 2M
       - Asymptotically Minkowski

    2. Tortoise (t, r*): r* = r + 2M ln|r/2M - 1|
       - Horizon at r* → -∞
       - Light rays at 45° in (t, r*) plane

    3. Eddington-Finkelstein (v, r) or (u, r):
       - Regular at horizon
       - Ingoing: v = t + r*
       - Outgoing: u = t - r*

    4. Kruskal-Szekeres (U, V):
       - Regular everywhere except singularity
       - Maximal extension of spacetime
       - Light rays at 45°

    5. Penrose (U', V'):
       - Compactified Kruskal
       - Entire spacetime in finite region
       - Light rays at 45°
       - I⁺, I⁻, i⁰, i⁺, i⁻ at finite distance

All functions assume geometric units with G = c = 1.
"""


import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq


def schwarzschild_to_tortoise(
    r: float | NDArray,
    M: float = 1.0
) -> float | NDArray:
    """
    Convert Schwarzschild r to tortoise coordinate r*.

    r* = r + 2M ln|r/(2M) - 1|

    Properties:
        - r* → r as r → ∞
        - r* → -∞ as r → 2M (horizon)
        - r* ∈ (-∞, ∞) for r ∈ (2M, ∞)

    Args:
        r: Schwarzschild radial coordinate (must be > 2M)
        M: Black hole mass

    Returns:
        Tortoise coordinate r*
    """
    r = np.asarray(r)
    rs = 2 * M

    # Handle r ≤ 2M (inside horizon)
    mask_outside = r > rs
    r_star = np.zeros_like(r, dtype=np.float64)

    if np.any(mask_outside):
        r_outside = r[mask_outside] if r.ndim > 0 else r
        r_star_outside = r_outside + rs * np.log(r_outside / rs - 1)
        if r.ndim > 0:
            r_star[mask_outside] = r_star_outside
        else:
            r_star = r_star_outside

    # Inside horizon (r < 2M): use modified formula
    mask_inside = (r > 0) & (r < rs)
    if np.any(mask_inside):
        r_inside = r[mask_inside] if r.ndim > 0 else r
        # r* = r + 2M ln|1 - r/(2M)|
        r_star_inside = r_inside + rs * np.log(np.abs(1 - r_inside / rs))
        if r.ndim > 0:
            r_star[mask_inside] = r_star_inside

    return float(r_star) if r.ndim == 0 else r_star


def tortoise_to_schwarzschild(
    r_star: float,
    M: float = 1.0,
    r_guess: float | None = None
) -> float:
    """
    Convert tortoise r* to Schwarzschild r (numerical inversion).

    Uses root-finding since the inverse is transcendental.

    Args:
        r_star: Tortoise coordinate
        M: Black hole mass
        r_guess: Initial guess for r (optional)

    Returns:
        Schwarzschild radial coordinate r
    """
    rs = 2 * M

    def f(r):
        return schwarzschild_to_tortoise(r, M) - r_star

    # Determine search bounds based on r_star
    if r_star > 0:
        # Well outside horizon
        r_low = rs + 1e-6
        r_high = r_star + 10 * M
    else:
        # Near or at horizon (r_star very negative)
        r_low = rs + 1e-10
        r_high = rs + 10 * M * np.exp(r_star / rs)

    try:
        return brentq(f, r_low, r_high)
    except ValueError:
        # Expand search
        return brentq(f, rs + 1e-15, 1000 * M)


def schwarzschild_to_kruskal(
    t: float | NDArray,
    r: float | NDArray,
    M: float = 1.0,
    region: str = "exterior"
) -> tuple[float | NDArray, float | NDArray]:
    """
    Convert Schwarzschild (t, r) to Kruskal-Szekeres (U, V).

    For r > 2M (exterior):
        U = √(r/2M - 1) exp(r/4M) cosh(t/4M)
        V = √(r/2M - 1) exp(r/4M) sinh(t/4M)

    For r < 2M (interior):
        U = √(1 - r/2M) exp(r/4M) sinh(t/4M)
        V = √(1 - r/2M) exp(r/4M) cosh(t/4M)

    Kruskal properties:
        - Metric: ds² = (32M³/r) e^(-r/2M) (-dV² + dU²) + r² dΩ²
        - Horizon at U = 0 or V = 0
        - Singularity at UV = 1
        - Light rays at 45° (dV = ±dU)

    Args:
        t: Schwarzschild time coordinate
        r: Schwarzschild radial coordinate
        M: Black hole mass
        region: "exterior" (r > 2M), "interior" (r < 2M), or "auto"

    Returns:
        (U, V) Kruskal coordinates
    """
    t = np.asarray(t, dtype=np.float64)
    r = np.asarray(r, dtype=np.float64)
    rs = 2 * M

    # Determine region
    if region == "auto":
        is_exterior = r > rs
    else:
        is_exterior = region == "exterior"  # type: ignore[assignment]

    # Broadcast to same shape
    if t.shape != r.shape:
        t, r = np.broadcast_arrays(t, r)

    U = np.zeros_like(r, dtype=np.float64)
    V = np.zeros_like(r, dtype=np.float64)

    if np.isscalar(is_exterior):
        if is_exterior:
            # Exterior: r > 2M
            factor = np.sqrt(r / rs - 1) * np.exp(r / (2 * rs))
            U = factor * np.cosh(t / (2 * rs))
            V = factor * np.sinh(t / (2 * rs))
        else:
            # Interior: r < 2M
            factor = np.sqrt(1 - r / rs) * np.exp(r / (2 * rs))
            U = factor * np.sinh(t / (2 * rs))
            V = factor * np.cosh(t / (2 * rs))
    else:
        # Array with mixed regions
        ext_mask = is_exterior
        int_mask = ~is_exterior & (r > 0)

        if np.any(ext_mask):
            r_ext = r[ext_mask]
            t_ext = t[ext_mask]
            factor = np.sqrt(r_ext / rs - 1) * np.exp(r_ext / (2 * rs))
            U[ext_mask] = factor * np.cosh(t_ext / (2 * rs))
            V[ext_mask] = factor * np.sinh(t_ext / (2 * rs))

        if np.any(int_mask):
            r_int = r[int_mask]
            t_int = t[int_mask]
            factor = np.sqrt(1 - r_int / rs) * np.exp(r_int / (2 * rs))
            U[int_mask] = factor * np.sinh(t_int / (2 * rs))
            V[int_mask] = factor * np.cosh(t_int / (2 * rs))

    # Return scalar if input was scalar
    if t.ndim == 0 and r.ndim == 0:
        return float(U), float(V)

    return U, V


def kruskal_to_schwarzschild(
    U: float,
    V: float,
    M: float = 1.0
) -> tuple[float, float, str]:
    """
    Convert Kruskal (U, V) to Schwarzschild (t, r) and determine region.

    Uses the relations:
        UV = (r/2M - 1) e^(r/2M)  for exterior
        UV = (1 - r/2M) e^(r/2M)  for interior
        V/U = tanh(t/4M)          for exterior

    Args:
        U: Kruskal U coordinate
        V: Kruskal V coordinate
        M: Black hole mass

    Returns:
        (t, r, region) where region is "exterior", "interior", or "white_hole"
    """
    rs = 2 * M
    UV = U * V

    # Determine region
    if V > abs(U):
        # Region II (interior of black hole)
        region = "interior"
    elif V < -abs(U):
        # Region IV (white hole interior)
        region = "white_hole"
    elif U > 0:
        # Region I (our exterior)
        region = "exterior"
    else:
        # Region III (other exterior)
        region = "other_exterior"

    # Find r from UV = (r/2M - 1) e^(r/2M) [exterior] or similar
    # This requires numerical inversion
    def f(r):
        if region in ["exterior", "other_exterior"]:
            return (r / rs - 1) * np.exp(r / rs) - UV
        else:
            return (1 - r / rs) * np.exp(r / rs) - UV

    # Solve for r
    try:
        if region in ["exterior", "other_exterior"]:
            r = brentq(f, rs + 1e-10, 1000 * M)
        else:
            r = brentq(f, 1e-10, rs - 1e-10)
    except ValueError:
        r = rs  # Fallback to horizon

    # Find t from V/U = tanh(t/4M) [exterior] or coth(t/4M) [interior]
    if abs(U) > 1e-14:
        ratio = V / U
        if region in ["exterior", "other_exterior"]:
            t = 2 * rs * np.arctanh(np.clip(ratio, -0.999, 0.999))
        else:
            if abs(ratio) > 1:
                t = 2 * rs * np.arctanh(1 / np.clip(ratio, -100, 100))
            else:
                t = 0.0
    else:
        t = 0.0

    return t, r, region


def kruskal_to_penrose(
    U: float | NDArray,
    V: float | NDArray
) -> tuple[float | NDArray, float | NDArray]:
    """
    Compactify Kruskal coordinates to Penrose coordinates.

    U' = arctan(U)
    V' = arctan(V)

    This maps the infinite Kruskal plane to a finite region:
        U, V ∈ (-∞, ∞) → U', V' ∈ (-π/2, π/2)

    Properties:
        - Light rays still at 45° (null lines preserved)
        - I⁺ at V' = π/2 - U' (for U' < π/2)
        - I⁻ at V' = -π/2 + U'
        - i⁰ at U' = π/2, V' = 0
        - i⁺ at U' = 0, V' = π/2
        - i⁻ at U' = 0, V' = -π/2

    Args:
        U, V: Kruskal coordinates

    Returns:
        (U', V') Penrose (compactified) coordinates
    """
    U_prime = np.arctan(U)
    V_prime = np.arctan(V)
    return U_prime, V_prime


def penrose_to_kruskal(
    U_prime: float | NDArray,
    V_prime: float | NDArray
) -> tuple[float | NDArray, float | NDArray]:
    """
    Inverse of Penrose compactification.

    U = tan(U')
    V = tan(V')
    """
    U = np.tan(U_prime)
    V = np.tan(V_prime)
    return U, V


def schwarzschild_to_penrose(
    t: float | NDArray,
    r: float | NDArray,
    M: float = 1.0,
    region: str = "auto"
) -> tuple[float | NDArray, float | NDArray]:
    """
    Convert Schwarzschild (t, r) directly to Penrose (U', V').

    Composition: Schwarzschild → Kruskal → Penrose

    Args:
        t: Schwarzschild time
        r: Schwarzschild radius
        M: Black hole mass
        region: "exterior", "interior", or "auto"

    Returns:
        (U', V') Penrose coordinates
    """
    U, V = schwarzschild_to_kruskal(t, r, M, region)
    return kruskal_to_penrose(U, V)


def eddington_finkelstein_ingoing(
    t: float | NDArray,
    r: float | NDArray,
    M: float = 1.0
) -> float | NDArray:
    """
    Convert to ingoing Eddington-Finkelstein coordinate v.

    v = t + r*  where r* is the tortoise coordinate

    In (v, r) coordinates:
        - Ingoing null rays are lines of constant v
        - The metric is regular at the horizon
        - ds² = -(1-2M/r)dv² + 2dvdr + r²dΩ²

    Args:
        t: Schwarzschild time
        r: Schwarzschild radius
        M: Black hole mass

    Returns:
        Ingoing EF coordinate v
    """
    r_star = schwarzschild_to_tortoise(r, M)
    return t + r_star


def eddington_finkelstein_outgoing(
    t: float | NDArray,
    r: float | NDArray,
    M: float = 1.0
) -> float | NDArray:
    """
    Convert to outgoing Eddington-Finkelstein coordinate u.

    u = t - r*  where r* is the tortoise coordinate

    In (u, r) coordinates:
        - Outgoing null rays are lines of constant u
        - ds² = -(1-2M/r)du² - 2dudr + r²dΩ²

    Args:
        t: Schwarzschild time
        r: Schwarzschild radius
        M: Black hole mass

    Returns:
        Outgoing EF coordinate u
    """
    r_star = schwarzschild_to_tortoise(r, M)
    return t - r_star


def transform_geodesic_to_penrose(
    positions: NDArray[np.float64],
    M: float = 1.0
) -> NDArray[np.float64]:
    """
    Transform a geodesic trajectory from Schwarzschild to Penrose.

    Args:
        positions: Array of shape (N, 4) with [t, r, θ, φ]
        M: Black hole mass

    Returns:
        Array of shape (N, 4) with [U', V', θ, φ]
    """
    t = positions[:, 0]
    r = positions[:, 1]
    theta = positions[:, 2]
    phi = positions[:, 3]

    # Determine which points are exterior/interior
    rs = 2 * M
    is_exterior = r > rs

    U_prime = np.zeros_like(t)
    V_prime = np.zeros_like(t)

    # Transform exterior points
    if np.any(is_exterior):
        U_ext, V_ext = schwarzschild_to_kruskal(
            t[is_exterior], r[is_exterior], M, "exterior"
        )
        U_prime[is_exterior], V_prime[is_exterior] = kruskal_to_penrose(U_ext, V_ext)

    # Transform interior points
    is_interior = ~is_exterior & (r > 0)
    if np.any(is_interior):
        U_int, V_int = schwarzschild_to_kruskal(
            t[is_interior], r[is_interior], M, "interior"
        )
        U_prime[is_interior], V_prime[is_interior] = kruskal_to_penrose(U_int, V_int)

    result = np.column_stack([U_prime, V_prime, theta, phi])
    return result
