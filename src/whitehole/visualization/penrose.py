"""
Penrose diagram visualization.

Penrose diagrams (conformal diagrams) are the gold standard for visualizing
global causal structure. They compactify infinite spacetime to a finite region
while preserving:
    - Causal relationships (light rays at 45°)
    - The topology of infinity (I⁺, I⁻, i⁰, i⁺, i⁻)

This module provides tools to:
    1. Draw the conformal boundary (infinity)
    2. Draw horizons and singularities
    3. Transform and plot geodesics
    4. Overlay causal classification
"""

from typing import Optional, List, Tuple
import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.patches import Polygon, FancyArrowPatch
from matplotlib.collections import LineCollection

from whitehole.coordinates.transforms import (
    schwarzschild_to_penrose,
    transform_geodesic_to_penrose,
)
from whitehole.causal.geodesic import GeodesicResult


class PenroseDiagram:
    """
    Create and render Penrose diagrams for Schwarzschild spacetime.

    The diagram uses coordinates (U', V') where:
        U' = arctan(U_Kruskal)
        V' = arctan(V_Kruskal)

    Key features:
        - Future null infinity I⁺: V' = π/2 - U' for U' > 0
        - Past null infinity I⁻: V' = -π/2 + U' for U' > 0
        - Spatial infinity i⁰: (U', V') = (π/2, 0)
        - Future timelike infinity i⁺: (0, π/2)
        - Past timelike infinity i⁻: (0, -π/2)
        - Event horizon: U' = 0 (for V' > 0, black hole)
        - Singularity: curve where r = 0

    Example:
        >>> diagram = PenroseDiagram(M=1.0)
        >>> fig, ax = diagram.create_figure()
        >>> diagram.draw_conformal_boundary(ax)
        >>> diagram.draw_horizon(ax)
        >>> diagram.draw_singularity(ax)
        >>> diagram.add_geodesic(ax, geodesic_result)
        >>> plt.show()
    """

    def __init__(self, M: float = 1.0):
        """
        Initialize the Penrose diagram.

        Args:
            M: Black hole mass (determines horizon location in original coords)
        """
        self.M = M
        self.rs = 2 * M

        # Penrose coordinate limits
        self.u_min = -np.pi / 2
        self.u_max = np.pi / 2
        self.v_min = -np.pi / 2
        self.v_max = np.pi / 2

    def create_figure(
        self,
        figsize: Tuple[float, float] = (10, 10),
        title: str = "Penrose Diagram"
    ) -> Tuple[Figure, Axes]:
        """
        Create a figure for the Penrose diagram.

        Returns:
            (fig, ax) tuple
        """
        fig, ax = plt.subplots(figsize=figsize)

        ax.set_xlim(self.u_min - 0.1, self.u_max + 0.1)
        ax.set_ylim(self.v_min - 0.1, self.v_max + 0.1)
        ax.set_aspect('equal')
        ax.set_xlabel(r"$U'$ (compactified)", fontsize=12)
        ax.set_ylabel(r"$V'$ (compactified)", fontsize=12)
        ax.set_title(title, fontsize=14)

        return fig, ax

    def draw_conformal_boundary(
        self,
        ax: Axes,
        color: str = 'blue',
        linewidth: float = 2.0,
        label_infinity: bool = True
    ) -> None:
        """
        Draw the conformal boundary (null and spatial infinity).

        The boundary consists of:
            - I⁺ (future null infinity): where outgoing light escapes
            - I⁻ (past null infinity): where ingoing light comes from
            - i⁰ (spatial infinity): r → ∞ at fixed t
            - i⁺ (future timelike infinity): t → ∞ at fixed r
            - i⁻ (past timelike infinity): t → -∞ at fixed r
        """
        # Right side of diagram (Region I exterior)

        # I⁺ (future null infinity): from i⁰ to i⁺
        # In Penrose coords: U' + V' = π/2, U' ∈ [0, π/2]
        u_plus = np.linspace(0, np.pi/2, 100)
        v_plus = np.pi/2 - u_plus
        ax.plot(u_plus, v_plus, color=color, linewidth=linewidth, label=r'$\mathscr{I}^+$')

        # I⁻ (past null infinity): from i⁻ to i⁰
        # In Penrose coords: V' - U' = -π/2, U' ∈ [0, π/2]
        u_minus = np.linspace(0, np.pi/2, 100)
        v_minus = u_minus - np.pi/2
        ax.plot(u_minus, v_minus, color=color, linewidth=linewidth, label=r'$\mathscr{I}^-$')

        # Mark special points
        if label_infinity:
            # i⁰ (spatial infinity)
            ax.plot(np.pi/2, 0, 'o', color=color, markersize=8)
            ax.annotate(r'$i^0$', (np.pi/2 + 0.05, 0), fontsize=12, color=color)

            # i⁺ (future timelike infinity)
            ax.plot(0, np.pi/2, '^', color=color, markersize=8)
            ax.annotate(r'$i^+$', (0.05, np.pi/2), fontsize=12, color=color)

            # i⁻ (past timelike infinity)
            ax.plot(0, -np.pi/2, 'v', color=color, markersize=8)
            ax.annotate(r'$i^-$', (0.05, -np.pi/2), fontsize=12, color=color)

    def draw_horizon(
        self,
        ax: Axes,
        color: str = 'red',
        linewidth: float = 2.0,
        linestyle: str = '--'
    ) -> None:
        """
        Draw the event horizon.

        The event horizon in Kruskal is at U = 0 (for V > 0, future horizon)
        and V = 0 (for U < 0, past horizon).

        In Penrose coordinates: U' = 0 line for V' > 0.
        """
        # Future event horizon: U' = 0, V' ∈ [0, V'_singularity]
        # The singularity in Penrose is at UV = 1, which maps to a curve
        # For the horizon (U=0), it extends from V=0 to V→∞ (the singularity)
        v_horizon = np.linspace(0, np.arctan(10), 100)  # Approximate
        u_horizon = np.zeros_like(v_horizon)
        ax.plot(u_horizon, v_horizon, color=color, linewidth=linewidth,
                linestyle=linestyle, label='Event Horizon')

        # Past horizon (white hole): V' = 0, U' ∈ [U'_singularity, 0]
        u_past = np.linspace(-np.arctan(10), 0, 100)
        v_past = np.zeros_like(u_past)
        ax.plot(u_past, v_past, color=color, linewidth=linewidth,
                linestyle=linestyle, alpha=0.5)

    def draw_singularity(
        self,
        ax: Axes,
        color: str = 'black',
        linewidth: float = 3.0
    ) -> None:
        """
        Draw the singularity (r = 0).

        In Kruskal, r = 0 corresponds to UV = 1.
        In Penrose: tan(U') tan(V') = 1
        """
        # Future singularity (inside black hole)
        # Parametrize: U = 1/V, then V from small to large
        V_values = np.linspace(1.01, 50, 200)
        U_values = 1.0 / V_values

        U_prime = np.arctan(U_values)
        V_prime = np.arctan(V_values)

        ax.plot(U_prime, V_prime, color=color, linewidth=linewidth,
                label='Singularity (r=0)')

        # Past singularity (white hole)
        ax.plot(-U_prime, -V_prime, color=color, linewidth=linewidth, alpha=0.5)

    def draw_constant_r_curves(
        self,
        ax: Axes,
        r_values: List[float] = [3.0, 4.0, 6.0, 10.0],
        t_range: Tuple[float, float] = (-20, 20),
        n_points: int = 100,
        color: str = 'gray',
        alpha: float = 0.5
    ) -> None:
        """
        Draw curves of constant Schwarzschild r.

        These are timelike curves (inside horizon) or spacelike curves
        (the singularity).
        """
        t_values = np.linspace(*t_range, n_points)

        for r in r_values:
            if r <= 0 or r == self.rs:
                continue

            try:
                U_prime, V_prime = schwarzschild_to_penrose(
                    t_values, np.full_like(t_values, r), self.M
                )
                ax.plot(U_prime, V_prime, color=color, alpha=alpha,
                        linewidth=0.5, label=f'r={r}M' if r == r_values[0] else '')
            except Exception:
                continue

    def draw_constant_t_curves(
        self,
        ax: Axes,
        t_values: List[float] = [-10, -5, 0, 5, 10],
        r_range: Tuple[float, float] = (2.01, 50),
        n_points: int = 100,
        color: str = 'lightgray',
        alpha: float = 0.3
    ) -> None:
        """
        Draw curves of constant Schwarzschild t.

        These are spacelike surfaces.
        """
        r_values = np.linspace(*r_range, n_points)

        for t in t_values:
            try:
                U_prime, V_prime = schwarzschild_to_penrose(
                    np.full_like(r_values, t), r_values, self.M
                )
                ax.plot(U_prime, V_prime, color=color, alpha=alpha,
                        linewidth=0.5)
            except Exception:
                continue

    def add_geodesic(
        self,
        ax: Axes,
        result: GeodesicResult,
        color: str = 'green',
        linewidth: float = 1.5,
        label: Optional[str] = None,
        marker_start: bool = True,
        marker_end: bool = True
    ) -> None:
        """
        Add a geodesic trajectory to the diagram.

        Args:
            ax: Matplotlib axes
            result: GeodesicResult from the tracer
            color: Line color
            linewidth: Line width
            label: Optional label for legend
            marker_start: Mark starting point
            marker_end: Mark ending point
        """
        # Transform to Penrose coordinates
        penrose_coords = transform_geodesic_to_penrose(result.positions, self.M)
        U_prime = penrose_coords[:, 0]
        V_prime = penrose_coords[:, 1]

        # Filter valid points (in the plottable region)
        valid = (np.abs(U_prime) < np.pi/2 - 0.01) & (np.abs(V_prime) < np.pi/2 - 0.01)

        if np.any(valid):
            ax.plot(U_prime[valid], V_prime[valid], color=color,
                    linewidth=linewidth, label=label)

            if marker_start and valid[0]:
                ax.plot(U_prime[0], V_prime[0], 'o', color=color, markersize=6)

            if marker_end and valid[-1]:
                ax.plot(U_prime[valid][-1], V_prime[valid][-1], 's',
                        color=color, markersize=6)

    def add_light_cone(
        self,
        ax: Axes,
        position: NDArray[np.float64],
        size: float = 0.2,
        color: str = 'yellow',
        alpha: float = 0.3
    ) -> None:
        """
        Draw a light cone at a given position.

        In Penrose coordinates, light cones are 45° lines (since null
        rays have dU' = ±dV').

        Args:
            ax: Matplotlib axes
            position: Schwarzschild coordinates [t, r, θ, φ]
            size: Size of light cone in Penrose coordinates
            color: Fill color
            alpha: Transparency
        """
        t, r = position[0], position[1]

        try:
            U_p, V_p = schwarzschild_to_penrose(t, r, self.M)
        except Exception:
            return

        # Light cone vertices (future and past)
        future_cone = np.array([
            [U_p, V_p],
            [U_p + size, V_p + size],
            [U_p - size, V_p + size],
        ])

        past_cone = np.array([
            [U_p, V_p],
            [U_p + size, V_p - size],
            [U_p - size, V_p - size],
        ])

        ax.add_patch(Polygon(future_cone, color=color, alpha=alpha))
        ax.add_patch(Polygon(past_cone, color=color, alpha=alpha * 0.5))
        ax.plot(U_p, V_p, 'ko', markersize=4)

    def shade_trapped_region(
        self,
        ax: Axes,
        color: str = 'lightcoral',
        alpha: float = 0.2
    ) -> None:
        """
        Shade the trapped region (inside event horizon).

        This is the region where all future-directed causal curves
        end at the singularity.
        """
        # The trapped region is bounded by:
        # - The future horizon (U' = 0, V' > 0)
        # - The future singularity (tan(U')tan(V') = 1)
        # - The boundary where U' < 0

        # Create polygon vertices
        n = 100
        V_horizon = np.linspace(0, np.arctan(50), n)
        U_horizon = np.zeros(n)

        # Singularity curve
        V_sing_vals = np.linspace(50, 1.01, n)
        U_sing_vals = 1.0 / V_sing_vals
        U_sing = np.arctan(U_sing_vals)
        V_sing = np.arctan(V_sing_vals)

        # Combine to form polygon
        vertices = np.vstack([
            np.column_stack([U_horizon, V_horizon]),
            np.column_stack([U_sing, V_sing]),
        ])

        polygon = Polygon(vertices, color=color, alpha=alpha, label='Trapped Region')
        ax.add_patch(polygon)

    def create_standard_diagram(
        self,
        geodesics: Optional[List[GeodesicResult]] = None,
        show_grid: bool = True,
        title: str = "Schwarzschild Penrose Diagram"
    ) -> Tuple[Figure, Axes]:
        """
        Create a complete standard Penrose diagram.

        Args:
            geodesics: Optional list of geodesics to plot
            show_grid: Show constant r and t curves
            title: Figure title

        Returns:
            (fig, ax) tuple
        """
        fig, ax = self.create_figure(title=title)

        # Draw structure
        self.shade_trapped_region(ax)
        self.draw_conformal_boundary(ax)
        self.draw_horizon(ax)
        self.draw_singularity(ax)

        if show_grid:
            self.draw_constant_r_curves(ax, r_values=[3, 4, 6, 10, 20])
            self.draw_constant_t_curves(ax, t_values=[-10, -5, 0, 5, 10])

        # Add geodesics
        if geodesics:
            colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(geodesics)))
            for i, geo in enumerate(geodesics):
                escaped = geo.escaped()
                self.add_geodesic(
                    ax, geo,
                    color='blue' if escaped else 'red',
                    label=f'Ray {i+1} ({"escaped" if escaped else "trapped"})'
                )

        ax.legend(loc='upper left', fontsize=9)
        ax.grid(True, alpha=0.2)

        return fig, ax
