#!/usr/bin/env python3
"""
Whitehole Demo: Schwarzschild Spacetime

This script demonstrates the core functionality of Whitehole:
    1. Creating a Schwarzschild spacetime
    2. Tracing null geodesics (light rays)
    3. Classifying points as escaping or trapped
    4. Visualizing in multiple coordinate systems

Run this script to verify your installation and see causal structure
visualization in action.

Usage:
    python examples/demo_schwarzschild.py
"""

import sys
import numpy as np

# Use non-interactive backend for running without display
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Add src to path for development
sys.path.insert(0, 'src')

from whitehole.spacetimes.schwarzschild import Schwarzschild
from whitehole.causal.geodesic import NullGeodesicTracer
from whitehole.causal.classifier import CausalClassifier
from whitehole.visualization.plotter import SpacetimePlotter
from whitehole.visualization.penrose import PenroseDiagram


def demo_basic_geodesics():
    """
    Demo 1: Trace basic null geodesics and visualize.

    Shows outgoing rays escaping and ingoing rays falling in.
    """
    print("=" * 60)
    print("Demo 1: Basic Null Geodesics")
    print("=" * 60)

    # Create spacetime
    spacetime = Schwarzschild(M=1.0)
    print(f"Created {spacetime.name}")
    print(f"  Horizon radius: r_s = {spacetime.r_s}M")
    print(f"  Photon sphere: r = {spacetime.photon_sphere_radius()}M")
    print(f"  ISCO: r = {spacetime.isco_radius()}M")

    # Create tracer
    tracer = NullGeodesicTracer(spacetime)

    # Trace several geodesics
    print("\nTracing null geodesics...")

    # Outgoing from various radii
    results = []
    for r_start in [3.0, 5.0, 10.0]:
        result = tracer.trace_radial(r_start=r_start, outgoing=True)
        results.append(result)
        print(f"  Outgoing from r={r_start}M: {result.termination.value}, "
              f"final r={result.final_r:.1f}M")

    # Ingoing from various radii
    for r_start in [3.0, 5.0, 10.0]:
        result = tracer.trace_radial(r_start=r_start, outgoing=False)
        results.append(result)
        print(f"  Ingoing from r={r_start}M: {result.termination.value}, "
              f"final r={result.final_r:.4f}M")

    # Visualize
    print("\nCreating visualization...")
    plotter = SpacetimePlotter(M=1.0)
    for i, result in enumerate(results):
        direction = "out" if i < 3 else "in"
        r_start = [3, 5, 10][i % 3]
        plotter.add_geodesic(result, label=f"r₀={r_start}M ({direction})")

    fig = plotter.plot_all_coordinates()
    fig.savefig('examples/output/demo1_geodesics.png', dpi=150, bbox_inches='tight')
    print("  Saved: examples/output/demo1_geodesics.png")

    return results


def demo_penrose_diagram():
    """
    Demo 2: Create a Penrose diagram with geodesics.

    The Penrose diagram compactifies infinite spacetime to show
    global causal structure.
    """
    print("\n" + "=" * 60)
    print("Demo 2: Penrose Diagram")
    print("=" * 60)

    spacetime = Schwarzschild(M=1.0)
    tracer = NullGeodesicTracer(spacetime)

    # Trace geodesics
    geodesics = []

    # Escaping rays
    for r in [3.0, 4.0, 6.0, 10.0]:
        result = tracer.trace_radial(r_start=r, outgoing=True)
        geodesics.append(result)

    # Captured rays
    for r in [3.0, 5.0, 8.0]:
        result = tracer.trace_radial(r_start=r, outgoing=False)
        geodesics.append(result)

    # Create Penrose diagram
    diagram = PenroseDiagram(M=1.0)
    fig, ax = diagram.create_standard_diagram(
        geodesics=geodesics,
        show_grid=True,
        title="Schwarzschild Penrose Diagram with Null Geodesics"
    )

    fig.savefig('examples/output/demo2_penrose.png', dpi=150, bbox_inches='tight')
    print("  Saved: examples/output/demo2_penrose.png")


def demo_causal_classification():
    """
    Demo 3: Classify spacetime points as escaping or trapped.

    This demonstrates the core algorithm for horizon detection:
    launch rays from each point and see if any escape.
    """
    print("\n" + "=" * 60)
    print("Demo 3: Causal Point Classification")
    print("=" * 60)

    spacetime = Schwarzschild(M=1.0)
    classifier = CausalClassifier(
        spacetime,
        n_theta=2,
        n_phi=2,
        escape_radius=50.0
    )

    # Classify points along radial slice
    print("\nClassifying points along radial slice (t=0)...")
    r_values = np.linspace(1.5, 5.0, 8)

    for r in r_values:
        result = classifier.classify_point(np.array([0.0, r, np.pi/2, 0.0]))
        status = "[OK] ESCAPING" if result.classification.value == "escaping" else "[X] TRAPPED"
        print(f"  r = {r:.2f}M: {status} (escape fraction: {result.escape_fraction:.0%})")

    # Find horizon numerically
    print("\nFinding horizon radius by bisection...")
    r_horizon, uncertainty = classifier.find_horizon_radius(
        r_min=1.5,
        r_max=3.0,
        tolerance=0.01
    )
    print(f"  Computed horizon: r = {r_horizon:.4f}M ± {uncertainty:.4f}M")
    print(f"  Exact horizon:    r = {spacetime.r_s:.4f}M")
    print(f"  Error: {abs(r_horizon - spacetime.r_s):.4f}M")


def demo_classification_grid():
    """
    Demo 4: Create a 2D classification grid.

    Visualizes which regions of spacetime can send signals to infinity.
    """
    print("\n" + "=" * 60)
    print("Demo 4: 2D Classification Grid")
    print("=" * 60)

    spacetime = Schwarzschild(M=1.0)
    classifier = CausalClassifier(spacetime, n_theta=2, n_phi=2, escape_radius=50.0)
    plotter = SpacetimePlotter(M=1.0)

    # Create grid (smaller for speed)
    print("\nClassifying 2D grid (this may take a moment)...")
    t_values = np.linspace(-5, 10, 6)
    r_values = np.linspace(1.5, 6.0, 8)

    classification = classifier.classify_grid(t_values, r_values)

    # Count results
    n_escaping = np.sum(classification == 1)
    n_trapped = np.sum(classification == -1)
    n_boundary = np.sum(classification == 0)
    print(f"  Escaping points: {n_escaping}")
    print(f"  Trapped points: {n_trapped}")
    print(f"  Boundary points: {n_boundary}")

    # Visualize
    fig, ax = plotter.plot_classification_grid(
        classification, t_values, r_values,
        title="Causal Classification: J⁻(I⁺) vs Trapped Region"
    )

    fig.savefig('examples/output/demo4_classification.png', dpi=150, bbox_inches='tight')
    print("  Saved: examples/output/demo4_classification.png")


def demo_light_cone():
    """
    Demo 5: Visualize light cones at different positions.

    Shows how light cones tilt as you approach the horizon.
    """
    print("\n" + "=" * 60)
    print("Demo 5: Light Cone Tilting Near Horizon")
    print("=" * 60)

    spacetime = Schwarzschild(M=1.0)
    tracer = NullGeodesicTracer(spacetime)

    # Trace light cones at different radii
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Schwarzschild coordinates
    ax1 = axes[0]
    ax1.set_xlabel(r'$r/M$')
    ax1.set_ylabel(r'$t/M$')
    ax1.set_title('Light Cones in Schwarzschild Coordinates')
    ax1.axvline(x=2.0, color='red', linestyle='--', label='Horizon')

    for r_start in [2.5, 3.0, 4.0, 6.0, 10.0]:
        # Outgoing ray
        out_result = tracer.trace_radial(r_start, outgoing=True, lambda_span=(0, 20))
        ax1.plot(out_result.positions[:, 1], out_result.positions[:, 0],
                 'b-', alpha=0.7)

        # Ingoing ray
        in_result = tracer.trace_radial(r_start, outgoing=False, lambda_span=(0, 20))
        ax1.plot(in_result.positions[:, 1], in_result.positions[:, 0],
                 'r-', alpha=0.7)

        # Mark starting point
        ax1.plot(r_start, 0, 'ko', markersize=5)

    ax1.set_xlim(0, 15)
    ax1.set_ylim(-10, 30)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Penrose diagram with same rays
    diagram = PenroseDiagram(M=1.0)
    ax2 = axes[1]
    ax2.set_xlabel(r"$U'$")
    ax2.set_ylabel(r"$V'$")
    ax2.set_title('Light Cones in Penrose Coordinates')

    # Draw structure
    diagram.shade_trapped_region(ax2)
    diagram.draw_conformal_boundary(ax2)
    diagram.draw_horizon(ax2)
    diagram.draw_singularity(ax2)

    for r_start in [2.5, 3.0, 4.0, 6.0, 10.0]:
        out_result = tracer.trace_radial(r_start, outgoing=True, lambda_span=(0, 20))
        in_result = tracer.trace_radial(r_start, outgoing=False, lambda_span=(0, 20))

        diagram.add_geodesic(ax2, out_result, color='blue', linewidth=1)
        diagram.add_geodesic(ax2, in_result, color='red', linewidth=1)

    ax2.set_xlim(-0.2, np.pi/2 + 0.1)
    ax2.set_ylim(-np.pi/2 - 0.1, np.pi/2 + 0.1)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig('examples/output/demo5_lightcones.png', dpi=150, bbox_inches='tight')
    print("  Saved: examples/output/demo5_lightcones.png")


def main():
    """Run all demos."""
    import os

    # Create output directory
    os.makedirs('examples/output', exist_ok=True)

    print("\n" + "=" * 60)
    print("WHITEHOLE DEMONSTRATION")
    print("Computational Framework for Causal Structure in GR")
    print("=" * 60)

    # Run demos
    demo_basic_geodesics()
    demo_penrose_diagram()
    demo_causal_classification()
    demo_classification_grid()
    demo_light_cone()

    print("\n" + "=" * 60)
    print("All demos complete!")
    print("Output saved to examples/output/")
    print("=" * 60)


if __name__ == "__main__":
    main()
