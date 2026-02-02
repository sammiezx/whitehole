"""
General spacetime plotter supporting multiple coordinate systems.

Provides side-by-side visualization in:
    - Schwarzschild coordinates (t, r)
    - Kruskal-Szekeres coordinates (U, V)
    - Penrose coordinates (U', V')
"""

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import NDArray

from whitehole.causal.geodesic import GeodesicResult
from whitehole.coordinates.transforms import (
    schwarzschild_to_kruskal,
    transform_geodesic_to_penrose,
)


class SpacetimePlotter:
    """
    Multi-coordinate spacetime visualization.

    Allows plotting geodesics and spacetime structure in multiple
    coordinate systems simultaneously, making the coordinate-independence
    of causal structure apparent.

    Example:
        >>> plotter = SpacetimePlotter(M=1.0)
        >>> plotter.add_geodesic(result, label="Outgoing ray")
        >>> fig = plotter.plot_all_coordinates()
        >>> plt.show()
    """

    def __init__(self, M: float = 1.0):
        """
        Initialize the plotter.

        Args:
            M: Black hole mass
        """
        self.M = M
        self.rs = 2 * M
        self.geodesics: list[tuple[GeodesicResult, dict[str, Any]]] = []

    def add_geodesic(
        self,
        result: GeodesicResult,
        label: str | None = None,
        color: str | None = None,
        **kwargs
    ) -> None:
        """
        Add a geodesic to be plotted.

        Args:
            result: GeodesicResult from tracer
            label: Legend label
            color: Line color (auto-selected if None)
            **kwargs: Additional plotting arguments
        """
        style = {
            'label': label,
            'color': color,
            **kwargs
        }
        self.geodesics.append((result, style))

    def clear_geodesics(self) -> None:
        """Remove all geodesics."""
        self.geodesics = []

    def plot_schwarzschild(
        self,
        ax: Axes,
        show_horizon: bool = True,
        t_range: tuple[float, float] = (-20, 50),
        r_range: tuple[float, float] = (0.5, 15)
    ) -> None:
        """
        Plot in Schwarzschild coordinates (t, r).

        Note: These coordinates are singular at r = 2M.
        """
        ax.set_xlabel(r'$r/M$', fontsize=11)
        ax.set_ylabel(r'$t/M$', fontsize=11)
        ax.set_title('Schwarzschild Coordinates', fontsize=12)
        ax.set_xlim(*r_range)
        ax.set_ylim(*t_range)

        # Draw horizon
        if show_horizon:
            ax.axvline(x=self.rs, color='red', linestyle='--',
                       linewidth=2, label='Horizon (r=2M)')

        # Plot geodesics
        colors = plt.cm.tab10(np.linspace(0, 1, max(len(self.geodesics), 1)))  # type: ignore[attr-defined]
        for i, (geo, style) in enumerate(self.geodesics):
            t = geo.positions[:, 0]
            r = geo.positions[:, 1]

            color = style.get('color') or colors[i]
            label = style.get('label')
            escaped = geo.escaped()

            ax.plot(r, t, color=color, linewidth=1.5, label=label)
            ax.plot(r[0], t[0], 'o', color=color, markersize=5)
            ax.plot(r[-1], t[-1], 's' if escaped else 'x',
                    color=color, markersize=6)

        ax.legend(loc='upper left', fontsize=8)
        ax.grid(True, alpha=0.3)

    def plot_kruskal(
        self,
        ax: Axes,
        show_horizon: bool = True,
        show_singularity: bool = True,
        uv_range: float = 3.0
    ) -> None:
        """
        Plot in Kruskal-Szekeres coordinates (U, V).

        In these coordinates:
            - The horizon is at U = 0 or V = 0
            - Light rays are at 45°
            - The singularity is at UV = 1
        """
        ax.set_xlabel(r'$U$', fontsize=11)
        ax.set_ylabel(r'$V$', fontsize=11)
        ax.set_title('Kruskal-Szekeres Coordinates', fontsize=12)
        ax.set_xlim(-uv_range, uv_range)
        ax.set_ylim(-uv_range, uv_range)
        ax.set_aspect('equal')

        # Draw horizon (U=0 and V=0 lines)
        if show_horizon:
            ax.axhline(y=0, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
            ax.axvline(x=0, color='red', linestyle='--', linewidth=1.5, alpha=0.7)

        # Draw singularity (UV = 1 hyperbola)
        if show_singularity:
            u_pos = np.linspace(0.1, uv_range, 100)
            v_pos = 1 / u_pos
            valid = v_pos < uv_range
            ax.plot(u_pos[valid], v_pos[valid], 'k-', linewidth=2.5, label='Singularity')
            ax.plot(-u_pos[valid], -v_pos[valid], 'k-', linewidth=2.5, alpha=0.5)

        # Plot geodesics
        colors = plt.cm.tab10(np.linspace(0, 1, max(len(self.geodesics), 1)))  # type: ignore[attr-defined]
        for i, (geo, style) in enumerate(self.geodesics):
            t = geo.positions[:, 0]
            r = geo.positions[:, 1]

            # Transform to Kruskal
            try:
                # Handle mixed interior/exterior
                U, V = schwarzschild_to_kruskal(t, r, self.M, region="auto")
            except Exception:
                continue

            color = style.get('color') or colors[i]
            label = style.get('label')

            # Filter valid points
            valid = (np.abs(U) < uv_range) & (np.abs(V) < uv_range)
            if np.any(valid):
                ax.plot(U[valid], V[valid], color=color, linewidth=1.5, label=label)  # type: ignore[index]
                if valid[0]:
                    ax.plot(U[0], V[0], 'o', color=color, markersize=5)  # type: ignore[index]

        ax.legend(loc='upper left', fontsize=8)
        ax.grid(True, alpha=0.3)

    def plot_penrose(
        self,
        ax: Axes,
        show_structure: bool = True
    ) -> None:
        """
        Plot in Penrose (compactified) coordinates.

        The entire spacetime fits in a finite region with light rays at 45°.
        """
        ax.set_xlabel(r"$U'$", fontsize=11)
        ax.set_ylabel(r"$V'$", fontsize=11)
        ax.set_title('Penrose Diagram', fontsize=12)
        ax.set_xlim(-np.pi/2 - 0.1, np.pi/2 + 0.1)
        ax.set_ylim(-np.pi/2 - 0.1, np.pi/2 + 0.1)
        ax.set_aspect('equal')

        if show_structure:
            # Conformal boundary
            u_line = np.linspace(0, np.pi/2, 50)
            ax.plot(u_line, np.pi/2 - u_line, 'b-', linewidth=2, label=r"$\mathscr{I}^+$")
            ax.plot(u_line, u_line - np.pi/2, 'b-', linewidth=2, label=r"$\mathscr{I}^-$")

            # Horizon
            v_h = np.linspace(0, np.arctan(5), 50)
            ax.plot(np.zeros_like(v_h), v_h, 'r--', linewidth=2, label='Horizon')

            # Singularity
            V_sing = np.linspace(1.01, 20, 100)
            U_sing = 1.0 / V_sing
            ax.plot(np.arctan(U_sing), np.arctan(V_sing), 'k-', linewidth=2.5,
                    label='Singularity')

            # Mark infinities
            ax.plot(np.pi/2, 0, 'bo', markersize=6)
            ax.plot(0, np.pi/2, 'b^', markersize=6)
            ax.plot(0, -np.pi/2, 'bv', markersize=6)

        # Plot geodesics
        colors = plt.cm.tab10(np.linspace(0, 1, max(len(self.geodesics), 1)))  # type: ignore[attr-defined]
        for i, (geo, style) in enumerate(self.geodesics):
            try:
                penrose_coords = transform_geodesic_to_penrose(geo.positions, self.M)
                U_p = penrose_coords[:, 0]
                V_p = penrose_coords[:, 1]
            except Exception:
                continue

            color = style.get('color') or colors[i]
            label = style.get('label')

            valid = (np.abs(U_p) < np.pi/2 - 0.05) & (np.abs(V_p) < np.pi/2 - 0.05)
            if np.any(valid):
                ax.plot(U_p[valid], V_p[valid], color=color, linewidth=1.5, label=label)
                if valid[0]:
                    ax.plot(U_p[0], V_p[0], 'o', color=color, markersize=5)

        ax.legend(loc='upper left', fontsize=8)
        ax.grid(True, alpha=0.3)

    def plot_all_coordinates(
        self,
        figsize: tuple[float, float] = (15, 5),
        title: str = "Null Geodesics in Multiple Coordinate Systems"
    ) -> Figure:
        """
        Create a figure with all three coordinate systems side by side.

        Returns:
            Matplotlib Figure
        """
        fig, axes = plt.subplots(1, 3, figsize=figsize)
        fig.suptitle(title, fontsize=14)

        self.plot_schwarzschild(axes[0])
        self.plot_kruskal(axes[1])
        self.plot_penrose(axes[2])

        plt.tight_layout()
        return fig

    def plot_classification_grid(
        self,
        classification: NDArray[np.int32],
        t_values: NDArray[np.float64],
        r_values: NDArray[np.float64],
        ax: Axes | None = None,
        title: str = "Causal Classification"
    ) -> tuple[Figure, Axes]:
        """
        Plot a 2D classification grid.

        Args:
            classification: 2D array from CausalClassifier.classify_grid
            t_values: Time coordinates
            r_values: Radial coordinates
            ax: Optional axes (creates new figure if None)
            title: Plot title

        Returns:
            (fig, ax) tuple
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))
        else:
            fig = ax.get_figure()  # type: ignore[assignment]

        # Custom colormap: red (trapped), yellow (boundary), green (escaping)
        from matplotlib.colors import ListedColormap
        colors = ['darkred', 'gray', 'yellow', 'green']
        cmap = ListedColormap(colors)

        # Plot
        r_grid, t_grid = np.meshgrid(r_values, t_values)
        mesh = ax.pcolormesh(r_grid, t_grid, classification,
                             cmap=cmap, vmin=-2, vmax=1, shading='auto')

        # Colorbar
        cbar = plt.colorbar(mesh, ax=ax, ticks=[-2, -1, 0, 1])
        cbar.ax.set_yticklabels(['Invalid', 'Trapped', 'Boundary', 'Escaping'])

        # Horizon line
        ax.axvline(x=self.rs, color='white', linestyle='--',
                   linewidth=2, label='Horizon')

        ax.set_xlabel(r'$r/M$', fontsize=11)
        ax.set_ylabel(r'$t/M$', fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.legend()

        return fig, ax
