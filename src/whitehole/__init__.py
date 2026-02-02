"""
Whitehole: Computational Framework for Global Causal Structure in General Relativity

A research-grade tool for computing and visualizing event horizons, causal structure,
and Penrose diagrams from spacetime data.
"""

__version__ = "0.1.0"

from whitehole.causal.classifier import CausalClassifier
from whitehole.causal.geodesic import NullGeodesicTracer
from whitehole.spacetimes.schwarzschild import Schwarzschild
from whitehole.spacetimes.vaidya import Vaidya

__all__ = [
    "Schwarzschild",
    "Vaidya",
    "NullGeodesicTracer",
    "CausalClassifier",
]
