# Development Roadmap

A phased execution plan for building Whitehole from prototype to research platform.

---

## Phase 1: Foundation & Prototype

**Goal:** Working null geodesic tracer with basic Penrose visualization

### Milestone 1.1: Core Infrastructure
- [ ] Set up project structure and dependencies
- [ ] Implement base `SpacetimeMetric` abstract class
- [ ] Create unit testing framework
- [ ] Set up continuous integration

### Milestone 1.2: Schwarzschild Spacetime
- [ ] Implement Schwarzschild metric tensor
- [ ] Implement inverse metric
- [ ] Compute and verify Christoffel symbols
- [ ] Test against known analytic results

### Milestone 1.3: Null Geodesic Engine
- [ ] Implement geodesic equation integrator (RK45)
- [ ] Null constraint enforcement
- [ ] Stopping conditions (singularity, escape)
- [ ] Verify against textbook geodesics

### Milestone 1.4: Basic Visualization
- [ ] Schwarzschild coordinate plotting
- [ ] Kruskal coordinate transformation
- [ ] Penrose compactification
- [ ] Static diagram generation

### Deliverable
A command-line tool that traces null geodesics through Schwarzschild spacetime and plots them in multiple coordinate systems.

---

## Phase 2: Causal Structure & Collapse

**Goal:** Compute J⁻(I⁺), detect horizons, visualize Vaidya collapse

### Milestone 2.1: Causal Classification
- [ ] Point classification algorithm (escaping vs trapped)
- [ ] Grid-based classification
- [ ] Boundary detection (horizon finding)
- [ ] Performance optimization (parallelization)

### Milestone 2.2: Vaidya Spacetime
- [ ] Implement Vaidya metric (ingoing coordinates)
- [ ] Configurable mass function M(v)
- [ ] Christoffel symbols for dynamical metric
- [ ] Test horizon formation

### Milestone 2.3: Retrodictive Horizon Construction
- [ ] Identify asymptotic null generators
- [ ] Backward ray tracing
- [ ] Horizon surface construction
- [ ] Animate horizon formation

### Milestone 2.4: Oppenheimer-Snyder Collapse
- [ ] Interior FLRW dust solution
- [ ] Exterior Schwarzschild matching
- [ ] Junction condition handling
- [ ] Full collapse visualization

### Deliverable
Animated visualization showing event horizon formation during gravitational collapse, demonstrating the global, teleological nature of the event horizon.

---

## Phase 3: Interactive Exploration

**Goal:** User-driven exploration of causal structure

### Milestone 3.1: Interactive Framework
- [ ] Point-and-click event selection
- [ ] Real-time ray launching
- [ ] Parameter adjustment (M, M(v))
- [ ] Responsive UI updates

### Milestone 3.2: Multiple Coordinate Views
- [ ] Side-by-side Schwarzschild/Kruskal/Penrose
- [ ] Synchronized highlighting
- [ ] Coordinate transformation overlays
- [ ] View switching controls

### Milestone 3.3: Causal Structure Tools
- [ ] Light cone visualization at selected point
- [ ] J⁺(p) and J⁻(p) shading
- [ ] Causal diamond construction
- [ ] Path connectivity checker

### Milestone 3.4: Educational Mode
- [ ] Guided tours of key concepts
- [ ] Interactive exercises
- [ ] Comparison with textbook diagrams
- [ ] Export for presentations

### Deliverable
Interactive web application where users can explore causal structure, launch light rays, and see real-time classification of spacetime regions.

---

## Phase 4: Research Extensions

**Goal:** Research-grade capabilities for publication

### Milestone 4.1: Apparent Horizon Detection
- [ ] Spatial slice definition
- [ ] Outgoing null expansion computation
- [ ] MOTS finder algorithm
- [ ] Apparent vs event horizon comparison

### Milestone 4.2: Dynamical Horizons
- [ ] Marginally trapped tube tracking
- [ ] Hayward horizon computation
- [ ] Isolated horizon identification
- [ ] Area evolution plots

### Milestone 4.3: Numerical Relativity Interface
- [ ] HDF5 metric data import
- [ ] Interpolation engine for grid data
- [ ] Time-dependent metric handling
- [ ] Interface with Einstein Toolkit / SpEC

### Milestone 4.4: QFT Overlay (Optional)
- [ ] Scalar field mode propagation
- [ ] Near-horizon redshift visualization
- [ ] Hawking temperature computation
- [ ] Mode scrambling near horizon

### Deliverable
A research platform capable of analyzing numerical relativity data, comparing horizon types, and supporting publication-quality figures.

---

## Phase 5: Advanced Features

**Goal:** Cutting-edge extensions

### Milestone 5.1: Kerr Spacetime
- [ ] Boyer-Lindquist metric implementation
- [ ] Ergosphere visualization
- [ ] Frame dragging geodesics
- [ ] Kerr-Penrose diagram

### Milestone 5.2: Causal Set Approximation
- [ ] Spacetime discretization (sprinkling)
- [ ] Causal relation computation
- [ ] Comparison with continuum
- [ ] Dimension estimators

### Milestone 5.3: GPU Acceleration
- [ ] Parallel geodesic integration
- [ ] Compute shader implementation
- [ ] Real-time classification
- [ ] Large-scale ray tracing

### Deliverable
Extended platform supporting rotating black holes, causal set theory comparisons, and high-performance computation.

---

## Phase 6: Publication & Release

**Goal:** Academic publication and open-source release

### Milestone 6.1: Paper Preparation
- [ ] Method documentation
- [ ] Benchmark comparisons
- [ ] Novel results highlighting
- [ ] Figure preparation

### Milestone 6.2: Open Source Release
- [ ] Code cleanup and documentation
- [ ] Example notebooks
- [ ] Installation guide
- [ ] Community contribution guidelines

### Milestone 6.3: Web Demo
- [ ] Hosted interactive demo
- [ ] Reduced-feature browser version
- [ ] Educational materials
- [ ] API documentation

### Deliverable
Published computational physics paper, open-source repository, and interactive web demo.

---

## Technical Debt & Maintenance

Throughout all phases, maintain:
- [ ] Comprehensive test coverage
- [ ] Performance benchmarks
- [ ] Documentation updates
- [ ] Code review standards

---

## Success Metrics

### Phase 1 Success
- Geodesics match analytic solutions to 6 decimal places
- Penrose diagram visually matches textbooks
- Code runs without crashes on standard inputs

### Phase 2 Success
- Horizon location matches known results (e.g., r = 2M for Schwarzschild)
- Vaidya horizon forms at correct advanced time
- Animation clearly shows global nature of horizon

### Phase 3 Success
- UI is responsive (<100ms latency)
- Users can discover horizon properties through exploration
- Multiple coordinate views are synchronized correctly

### Phase 4 Success
- Results match published numerical relativity papers
- Apparent/event horizon difference is quantitatively correct
- Platform accepts real NR data

### Phase 5 Success
- Kerr geodesics match known analytic/numerical results
- Causal set approximation converges to continuum
- GPU acceleration provides >10x speedup

### Phase 6 Success
- Paper accepted to peer-reviewed journal
- Open-source repository has active users
- Web demo receives traffic from GR community

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| Numerical instability near horizons | Use tortoise coordinates, adaptive stepping |
| Performance issues with ray tracing | GPU acceleration, spatial hashing |
| Coordinate singularities | Multiple coordinate patches, proper handling |
| Scope creep | Strict phase boundaries, MVP focus |
| Mathematical errors | Verification against analytic solutions |

---

## Dependencies & Tools

### Core Dependencies
- **Python 3.10+** or **TypeScript**
- **NumPy/SciPy** — Numerical computation
- **Matplotlib** — Static visualization
- **Three.js / WebGL** — Interactive 3D visualization

### Optional Dependencies
- **h5py** — HDF5 file reading
- **Numba / CUDA** — GPU acceleration
- **SymPy** — Symbolic verification
- **Jupyter** — Notebook exploration

### Development Tools
- **pytest** — Testing
- **mypy** — Type checking
- **black** — Code formatting
- **sphinx** — Documentation

---

## Getting Started Now

1. Clone repository
2. Set up Python environment
3. Implement Schwarzschild metric (Milestone 1.2)
4. Run first geodesic trace
5. Iterate toward Phase 1 deliverable

See [GETTING_STARTED.md](GETTING_STARTED.md) for detailed instructions.
