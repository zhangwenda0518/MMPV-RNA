# PyPopART

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python Version](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)

**PyPopART** is a pure Python implementation of PopART (Population Analysis with Reticulate Trees) for constructing and visualizing haplotype networks from DNA sequence data.

## Features

- **Multiple Network Algorithms**: MST, MSN, TCS (Statistical Parsimony), Median-Joining (MJN), Parsimony Network (PN), and Tight Span Walker (TSW)
- **Distance Metrics**: Hamming, Jukes-Cantor, Kimura 2-parameter, Tamura-Nei
- **Comprehensive Analysis**: Network statistics, topology analysis, population genetics measures
- **Rich Visualization**: Static (matplotlib) and interactive (Dash Cytoscape) network plots
- **Flexible I/O**: Support for FASTA, NEXUS, PHYLIP, GenBank formats
- **Command-Line Interface**: Easy-to-use CLI for all operations
- **Web-based GUI**: Interactive Dash application for network construction and visualization
- **Python API**: Programmatic access for custom workflows

## Installation

### From Source

```bash
git clone https://github.com/adamtaranto/pypopart.git
cd pypopart
pip install -e ".[dev]"
```

### Requirements

- Python 3.11 or higher
- Core dependencies: biopython, click, networkx, numpy
- Optional extras: `pypopart[viz]` (matplotlib, plotly), `pypopart[gui]` (Dash), `pypopart[speed]` (numba), `pypopart[all]` (everything)

## Quick Start

### Entry Points

PyPopART provides two main interfaces:

**1. Command-Line Interface (CLI)** - For scripting and batch processing:

```bash
pypopart --help
```

**2. Web-based GUI** - For interactive analysis:

```bash
pypopart-gui
# Opens web interface at http://localhost:8050
```

### Command-Line Interface

#### 1. Load and Validate Sequence Data

```bash
pypopart load sequences.fasta
```

Output:

```text
Loading sequences from sequences.fasta...
✓ Loaded 50 sequences
  Alignment length: 500 bp

Alignment Statistics:
  Sequences: 50
  Length: 500 bp
  Variable sites: 25
  Parsimony informative: 15
  GC content: 48.5%
```

#### 2. Construct a Haplotype Network

```bash
# Median-Joining Network (default)
pypopart network sequences.fasta -o network.graphml

# Statistical Parsimony (TCS)
pypopart network sequences.fasta -a tcs -o network.graphml

# With custom distance metric
pypopart network sequences.fasta -a mjn -d k2p -o network.graphml
```

#### 3. Analyze Network Statistics

```bash
pypopart analyze network.graphml --stats
```

Output:

```text
Loading network from network.graphml...
✓ Loaded network with 12 nodes

=== Network Statistics ===
nodes: 12
edges: 15
diameter: 4
avg_degree: 2.5
clustering_coefficient: 0.3214
reticulation_index: 0.25
```

#### 4. Visualize the Network

```bash
# Static plot (PNG/PDF/SVG)
pypopart visualize network.graphml -o network.png --layout spring

# Interactive HTML
pypopart visualize network.graphml -o network.html --interactive
```

#### 5. Get Information

```bash
# List available algorithms
pypopart info --list-algorithms

# Output:
# Available Network Construction Algorithms:
#   mst - Minimum Spanning Tree
#   msn - Minimum Spanning Network
#   tcs - Statistical Parsimony (TCS)
#   mjn - Median-Joining Network
#   pn  - Parsimony Network (consensus from multiple trees)
#   tsw - Tight Span Walker (metric-preserving network)

# List distance metrics
pypopart info --list-distances

# List supported formats
pypopart info --list-formats
```

### Web-based GUI

Launch the interactive Dash application:

```bash
# Start GUI on default port 8050
pypopart-gui

# Start on custom port
pypopart-gui --port 8080

# Enable debug mode
pypopart-gui --debug
```

Once started, open your browser to `http://localhost:8050` and follow the workflow:

1. **Upload Data**: Load sequence alignment (FASTA, NEXUS, or PHYLIP) and optional metadata (CSV)
2. **Configure Algorithm**: Choose network algorithm (MST, MSN, TCS, MJN, PN, TSW) and parameters
3. **Compute Network**: Build the haplotype network
4. **Customize Layout**: Adjust node positions, sizes, spacing, and layout algorithms
5. **Export Results**: Download network (GraphML, GML, JSON) or images (PNG, SVG)

**Features:**

- Interactive network visualization with zoom and pan
- Drag-and-drop node repositioning
- Population-based coloring (pie charts for mixed nodes)
- Search and highlight specific haplotypes
- Real-time statistics and haplotype summary
- Multiple layout algorithms (Spring, Hierarchical, Kamada-Kawai, etc.)

### Python API

```python
from pypopart.io import load_alignment
from pypopart.core.haplotype import identify_haplotypes_from_alignment
from pypopart.algorithms import build
from pypopart.visualization import StaticNetworkPlotter

# Load sequences
alignment = load_alignment('sequences.fasta')

# Identify unique haplotypes (informational; algorithms do this internally)
haplotypes = identify_haplotypes_from_alignment(alignment)
print(f'Found {len(haplotypes)} unique haplotypes')

# Construct a Median-Joining Network with K2P distances
mjn = build('mjn', distance_method='k2p', epsilon=0)
network = mjn.build_network(alignment)

# Visualize
plotter = StaticNetworkPlotter(network)
fig, ax = plotter.plot(layout_algorithm='spring')
fig.savefig('network.png', dpi=150, bbox_inches='tight')
```

## Network Construction Algorithms

### Minimum Spanning Tree (MST)

Creates a tree connecting all haplotypes with minimum total distance.

```bash
pypopart network sequences.fasta -a mst -o network.graphml
```

**Use when**: You want the simplest possible network structure without reticulation.

**Properties**: Always produces a tree (no cycles), guaranteed minimum total edge weight.

### Minimum Spanning Network (MSN)

Extends MST by adding alternative connections at equal distance.

```bash
pypopart network sequences.fasta -a msn -o network.graphml
```

**Use when**: You want to show alternative evolutionary pathways at the same genetic distance.

**Properties**: Includes all edges tied for minimum distance, may contain reticulations.

### TCS (Statistical Parsimony)

Connects haplotypes within a parsimony probability limit (default 95%).

```bash
pypopart network sequences.fasta -a tcs -p 0.95 -o network.graphml
```

**Use when**: You want statistically justified connections based on parsimony.

**Properties**: Uses connection limits based on parsimony probability, good for intraspecific data.

### Median-Joining Network (MJN)

Infers ancestral/median sequences and creates a reticulate network.

```bash
pypopart network sequences.fasta -a mjn -e 0 -o network.graphml
```

**Use when**: You want to infer ancestral haplotypes and show complex evolutionary relationships.

**Properties**:

- Infers median vectors (ancestral nodes)
- Epsilon parameter controls complexity (0 = maximum simplification)
- Handles reticulation and homoplasy
- Good for closely related sequences

### Parsimony Network (PN)

Creates a consensus network by sampling edges from multiple random parsimony trees.

```bash
pypopart network sequences.fasta -a pn -o network.graphml
```

**Use when**: You want a consensus approach that captures phylogenetic uncertainty across multiple tree topologies.

**Properties**:

- Samples 100 random parsimony trees by default
- Includes edges that appear frequently across trees
- Can represent reticulation where multiple edges have similar frequencies
- Automatically infers median vertices for multi-mutation edges
- Good for datasets with phylogenetic uncertainty

### Tight Span Walker (TSW)

Constructs networks using the tight span of the distance matrix, preserving all metric properties.

```bash
pypopart network sequences.fasta -a tsw -o network.graphml
```

**Use when**: You need accurate metric-preserving networks for complex evolutionary relationships with reticulation.

**Properties**:

- Preserves all metric properties of the distance matrix
- Automatically infers ancestral/median sequences
- Best for small to medium datasets (n < 100)
- Computationally intensive but highly accurate
- Handles reticulation and complex evolutionary patterns

## Distance Metrics

- **hamming**: Simple count of differences (fastest)
- **jc**: Jukes-Cantor correction for multiple substitutions
- **k2p**: Kimura 2-parameter (transitions vs transversions)
- **tamura_nei**: Accounts for GC content and transition/transversion bias

```bash
pypopart network sequences.fasta -d k2p -o network.graphml
```

## File Format Support

### Input Formats

- **FASTA** (`.fasta`, `.fa`, `.fna`)
- **NEXUS** (`.nexus`, `.nex`) - including traits/metadata
- **PHYLIP** (`.phy`, `.phylip`)
- **GenBank** (`.gb`, `.gbk`)

### Output Formats

- **GraphML** (`.graphml`) - Recommended, preserves all attributes
- **GML** (`.gml`)
- **JSON** (`.json`)
- **NEXUS** (`.nexus`, `.nex`)
- **PNG** (`.png`) - Raster image
- **SVG** (`.svg`) - Vector, web-friendly

## Advanced Usage

### Working with Metadata

```bash
# Load alignment with population metadata
pypopart load sequences.fasta -m metadata.csv

# Visualize the network
pypopart visualize network.graphml -o network.png --show-labels
```

Metadata CSV format:

```csv
id,population,latitude,longitude,color,notes
Hap1,PopA,,,,,
Hap2,PopA,,,,,
Hap3,PopB,,,,,
```

### Network Analysis

```bash
# Comprehensive statistics
pypopart analyze network.graphml --stats --topology --popgen -a sequences.fasta -o results.json
```

### Topology Analysis

```bash
pypopart analyze network.graphml --topology
```

Identifies:

- Connected components
- Star-like patterns
- Central/hub nodes
- Potential ancestral nodes

### Custom Visualization

```bash
# Circular layout with labels
pypopart visualize network.graphml -o network.pdf \
    --layout circular --show-labels --width 1200 --height 1200

# Kamada-Kawai layout, interactive
pypopart visualize network.graphml -o network.html \
    --layout kamada_kawai --interactive
```

## Documentation

Full documentation is available at [https://adamtaranto.github.io/pypopart](https://adamtaranto.github.io/pypopart/)

## Citation

If you use PyPopART in your research, please cite the original [PopART paper](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12410) as well as this repository:

```text
Leigh, J.W., Bryant, D. and Nakagawa, S., 2015. POPART: full-feature software for haplotype network construction. Methods in Ecology & Evolution, 6(9).

Taranto, A. (2025). PyPopART: Pure Python implementation of haplotype network analysis.
GitHub repository: https://github.com/adamtaranto/pypopart
```

### Algorithm References

PyPopART implements algorithms from the following publications:

- **Minimum Spanning Tree/Network**: Excoffier, L. & Smouse, P. E. (1994). Using allele frequencies and geographic subdivision to reconstruct gene trees within a species: molecular variance parsimony. _Genetics_, 136(1), 343-359.

- **TCS (Statistical Parsimony)**: Clement, M., Posada, D., & Crandall, K. A. (2000). TCS: a computer program to estimate gene genealogies. _Molecular Ecology_, 9(10), 1657-1659.

- **Median-Joining Network**: Bandelt, H. J., Forster, P., & Röhl, A. (1999). Median-joining networks for inferring intraspecific phylogenies. _Molecular Biology and Evolution_, 16(1), 37-48.

- **Parsimony Network**: Excoffier, L. & Smouse, P. E. (1994). Using allele frequencies and geographic subdivision to reconstruct gene trees within a species: molecular variance parsimony. _Genetics_, 136(1), 343-359.

- **Tight Span Walker**: Dress, A. W., Huber, K. T., Koolen, J., Moulton, V., & Spillner, A. (2012). _Basic Phylogenetic Combinatorics_. Cambridge University Press.

## License

PyPopART is licensed under the GNU General Public License v3.0 or later. See [LICENSE](LICENSE) for details.

### Development Setup

```bash
git clone https://github.com/adamtaranto/pypopart.git
cd pypopart
pip install -e ".[dev]"
pre-commit install
```

## Related Projects

- [PopART](https://github.com/jessicawleigh/popart-current) - Original PopART software

## Acknowledgments

PyPopART is a python port of the original PopART software developed by Jessica Leigh.
