# Contributing to Whitehole

Thank you for your interest in contributing to Whitehole! This project aims to advance computational visualization of causal structure in general relativity, and we welcome contributions from physicists, mathematicians, and software engineers alike.

## Getting Started

### Prerequisites

- Python 3.10 or higher
- Git

### Development Setup

1. Fork the repository on GitHub

2. Clone your fork:
   ```bash
   git clone https://github.com/YOUR_USERNAME/whitehole.git
   cd whitehole
   ```

3. Create a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

4. Install in development mode with all dependencies:
   ```bash
   pip install -e ".[dev,notebook]"
   ```

5. Run tests to verify setup:
   ```bash
   pytest
   ```

## Development Workflow

### Branching Strategy

- `main` — Stable release branch
- `develop` — Integration branch for features
- `feature/*` — New features
- `fix/*` — Bug fixes
- `docs/*` — Documentation improvements

### Making Changes

1. Create a new branch from `main`:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes, following our code style guidelines

3. Write or update tests as needed

4. Run the test suite:
   ```bash
   pytest
   ```

5. Run type checking:
   ```bash
   mypy src/whitehole
   ```

6. Run the linter:
   ```bash
   ruff check src/whitehole
   ```

7. Commit your changes with a descriptive message:
   ```bash
   git commit -m "Add feature: brief description"
   ```

8. Push to your fork and open a Pull Request

## Code Style

### Python Style

- Follow PEP 8 guidelines
- Use type hints for all function signatures
- Maximum line length: 88 characters (Black default)
- Use descriptive variable names

### Docstrings

Use NumPy-style docstrings for all public functions and classes:

```python
def compute_metric(r: float, theta: float) -> np.ndarray:
    """
    Compute the metric tensor at given coordinates.

    Parameters
    ----------
    r : float
        Radial coordinate (must be positive)
    theta : float
        Angular coordinate in radians

    Returns
    -------
    np.ndarray
        4x4 metric tensor g_μν

    Raises
    ------
    ValueError
        If r <= 0

    Examples
    --------
    >>> metric = compute_metric(10.0, np.pi/2)
    >>> metric.shape
    (4, 4)
    """
```

### Physics Conventions

- Use geometric units (G = c = 1) unless otherwise specified
- Metric signature: (−, +, +, +)
- Index notation: Greek indices (μ, ν) for spacetime (0-3), Latin indices (i, j) for space (1-3)
- Document coordinate systems explicitly

## Types of Contributions

### Bug Reports

- Use the GitHub issue tracker
- Include a minimal reproducible example
- Specify your Python version and OS
- Include the full error traceback

### Feature Requests

- Open an issue describing the feature
- Explain the use case and benefits
- Reference relevant physics literature if applicable

### Code Contributions

We especially welcome:

- New spacetime metrics (analytic or numerical)
- Improved geodesic integration algorithms
- New visualization modes
- Performance optimizations
- Documentation improvements
- Test coverage expansion

### Documentation

- Fix typos or clarify explanations
- Add examples and tutorials
- Improve API documentation
- Add references to physics literature

## Testing

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=whitehole

# Run specific test file
pytest tests/test_schwarzschild.py

# Run tests matching a pattern
pytest -k "geodesic"
```

### Writing Tests

- Place tests in the `tests/` directory
- Mirror the source structure (e.g., `test_schwarzschild.py` for `schwarzschild.py`)
- Test edge cases and boundary conditions
- Include physical validation tests where possible

## Pull Request Guidelines

1. **One feature per PR** — Keep changes focused
2. **Update documentation** — If you change behavior, update docs
3. **Add tests** — New features need tests
4. **Follow commit conventions** — Clear, descriptive commits
5. **Reference issues** — Link related issues in the PR description

### PR Title Format

- `feat: Add Kerr metric implementation`
- `fix: Correct horizon detection near singularity`
- `docs: Improve getting started guide`
- `test: Add integration tests for Vaidya spacetime`
- `refactor: Simplify geodesic tracer interface`

## Community

### Code of Conduct

We are committed to providing a welcoming and inclusive environment. Please be respectful and constructive in all interactions.

### Getting Help

- Open an issue for bugs or feature requests
- Start a discussion for questions or ideas
- Check existing issues before creating new ones

## Recognition

Contributors will be acknowledged in:
- The project README
- Release notes
- Any resulting publications (for significant contributions)

## License

By contributing to Whitehole, you agree that your contributions will be licensed under the MIT License.
