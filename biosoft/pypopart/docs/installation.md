# Installation

## Requirements

PyPopART requires:

- **Python 3.9 or higher**
- Operating system: Linux, macOS, or Windows

## Installation Methods

### From PyPI (Recommended)

Once released, install PyPopART using pip:

```bash
pip install pypopart
```

### From Source

For the latest development version:

```bash
# Clone the repository
git clone https://github.com/adamtaranto/pypopart.git
cd pypopart

# Install in development mode
pip install -e ".[dev]"
```

### Using conda/mamba

Create a conda environment (optional but recommended):

```bash
# Create environment
conda create -n pypopart python=3.11
conda activate pypopart

# Install pypopart
pip install pypopart
```

## Dependencies

PyPopART automatically installs the following dependencies:

### Core Dependencies

- **biopython**: Sequence I/O and manipulation
- **click**: Command-line interface
- **matplotlib**: Static visualization
- **networkx**: Graph data structures and algorithms
- **numba**: JIT compilation for performance
- **numpy**: Numerical operations
- **pandas**: Data manipulation
- **plotly**: Interactive visualization
- **scipy**: Scientific computing
- **scikit-learn**: Machine learning utilities

### Development Dependencies (optional)

For development, install additional tools:

```bash
pip install pypopart[dev]
```

This includes:

- **pytest**: Testing framework
- **pytest-cov**: Coverage reporting
- **ruff**: Fast linting and formatting
- **mypy**: Static type checking
- **pre-commit**: Git hooks for code quality

## Verify Installation

Check that PyPopART is installed correctly:

```bash
# Check version
pypopart --version

# Show help
pypopart --help

# List available commands
pypopart info --list-algorithms
```

Expected output:

```text
Available Network Construction Algorithms:
  mst - Minimum Spanning Tree
  msn - Minimum Spanning Network
  tcs - Statistical Parsimony (TCS)
  mjn - Median-Joining Network
```

## Testing Your Installation

Run a quick test with sample data:

```bash
# Create a simple FASTA file
cat > test.fasta << 'EOF'
>Seq1
ATCGATCG
>Seq2
ATCGATCG
>Seq3
ATCGATTG
EOF

# Construct a network
pypopart network test.fasta -o test.graphml

# Visualize it
pypopart visualize test.graphml -o test.png
```

If these commands succeed, PyPopART is working correctly!

## Python Version Compatibility

PyPopART is tested on:

- Python 3.9
- Python 3.10
- Python 3.11
- Python 3.12

## Platform Support

PyPopART works on:

- **Linux**: Fully supported
- **macOS**: Fully supported (both Intel and Apple Silicon)
- **Windows**: Fully supported (Windows 10+)

## Troubleshooting

### ImportError: No module named 'pypopart'

The package is not installed. Try:

```bash
pip install pypopart
```

### Command not found: pypopart

The pip install directory is not in your PATH. Try:

```bash
python -m pypopart --help
```

Or add pip's binary directory to your PATH.

### NumPy/SciPy Installation Issues

On some systems, you may need to install NumPy and SciPy separately first:

```bash
pip install numpy scipy
pip install pypopart
```

### Numba Compilation Warnings

Numba may show warnings on first import. These are normal and can be ignored. The code will work correctly even without JIT compilation.

### Visualization Issues on Headless Servers

If running on a server without a display:

```bash
# Set matplotlib backend
export MPLBACKEND=Agg
pypopart visualize network.graphml -o network.png
```

Or use interactive HTML output which doesn't require a display:

```bash
pypopart visualize network.graphml -o network.html --interactive
```

## Updating PyPopART

To update to the latest version:

```bash
pip install --upgrade pypopart
```

## Uninstalling

To remove PyPopART:

```bash
pip uninstall pypopart
```

## Getting Help

If you encounter issues:

1. Check the [FAQ](faq.md)
2. Search [GitHub Issues](https://github.com/adamtaranto/pypopart/issues)
3. Open a new issue with:
   - Your Python version (`python --version`)
   - Your OS and version
   - Complete error message
   - Steps to reproduce the problem

## Next Steps

- **[Quick Start](quickstart.md)**: Create your first haplotype network
- **[CLI Guide](guide/cli.md)**: Learn the command-line interface
- **[Tutorials](tutorials/basic_workflow.md)**: Step-by-step examples
