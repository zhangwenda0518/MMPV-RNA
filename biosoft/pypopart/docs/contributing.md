# Contributing to PyPopART

Thank you for your interest in contributing to PyPopART! This document provides guidelines for contributing to the project.

## Ways to Contribute

- Report bugs and issues
- Suggest new features or enhancements
- Improve documentation
- Submit bug fixes
- Add new algorithms or features
- Write tests

## Getting Started

### Development Setup

1. Fork the repository on GitHub
2. Clone your fork locally:

   ```bash
   git clone https://github.com/YOUR_USERNAME/pypopart.git
   cd pypopart
   ```

3. Create a virtual environment:

   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

4. Install in development mode with all dependencies:

   ```bash
   pip install -e ".[dev,test,docs]"
   ```

5. Install pre-commit hooks:
   ```bash
   pre-commit install
   ```

### Making Changes

1. Create a new branch for your changes:

   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes and write tests

3. Run tests to ensure everything works:

   ```bash
   pytest
   ```

4. Run linters and formatters:

   ```bash
   black src tests
   ruff check src tests
   mypy src
   ```

5. Commit your changes with clear messages:

   ```bash
   git commit -m "Add feature: your feature description"
   ```

6. Push to your fork and submit a pull request

## Code Style

- Follow [PEP 8](https://pep8.org/) style guidelines
- Use [Black](https://github.com/psf/black) for code formatting
- Use [Ruff](https://github.com/astral-sh/ruff) for linting
- Add type hints where appropriate
- Write docstrings in NumPy style

## Testing

- Write tests for all new features
- Ensure existing tests pass
- Aim for good test coverage
- Use pytest for testing

## Documentation

The docs are built with [Zensical](https://zensical.org), configured in
`zensical.toml`:

```bash
pip install -e ".[docs]"
zensical serve          # live preview on http://localhost:8000
zensical build --strict # what CI runs; fails on any warning
```

- Update documentation for new features
- Add docstrings to all public functions and classes
  (numpy style; `pre-commit` validates them with numpydoc)
- Include examples in docstrings
- Update the changelog

## Pull Request Process

1. Ensure all tests pass
2. Update documentation as needed
3. Add entry to CHANGELOG.md
4. Submit pull request with clear description
5. Address any review feedback

## Reporting Issues

When reporting issues, please include:

- PyPopART version
- Python version
- Operating system
- Minimal example to reproduce the issue
- Expected vs. actual behavior
- Error messages and stack traces

## Code of Conduct

- Be respectful and inclusive
- Welcome newcomers
- Focus on constructive feedback
- Help create a positive community

## Questions?

If you have questions, feel free to:

- Open an issue on GitHub
- Start a discussion in GitHub Discussions
- Contact the maintainers

Thank you for contributing to PyPopART!
