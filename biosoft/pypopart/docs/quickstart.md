# Quick Start Guide

Get started with PyPopART in minutes! This guide walks you through creating your first haplotype network.

## Prerequisites

- PyPopART installed (see [Installation](installation.md))
- A sequence alignment file (FASTA, NEXUS, PHYLIP, or GenBank format)

## Choosing Your Interface

PyPopART offers two ways to work with haplotype networks:

### Option 1: Web-based GUI (Recommended for Beginners)

Launch the interactive interface:

```bash
pypopart-gui
```

Open your browser to `http://localhost:8050` and follow the visual workflow:

1. Upload sequence file
2. Configure algorithm
3. Compute network
4. Customize layout
5. Export results

**Advantages**: Visual, interactive, no command-line knowledge needed.

### Option 2: Command-Line Interface

For automation and scripting, use the CLI as shown below.

## Your First Network in 3 Steps (CLI)

### Step 1: Prepare Your Data

Create or use an existing sequence alignment. For this example, we'll create a simple FASTA file:

```bash
cat > mysequences.fasta << 'EOF'
>Sample1_PopA
ATCGATCGATCGATCGATCG
>Sample2_PopA
ATCGATCGATCGATCGATCG
>Sample3_PopB
ATCGATCGATCGATTGATCG
>Sample4_PopB
ATCGATCGATCGATTGATCG
>Sample5_PopC
ATCGATCGAGCGATCGATCG
EOF
```

### Step 2: Construct the Network

Use the `network` command to build a haplotype network:

```bash
pypopart network mysequences.fasta -o mynetwork.graphml
```

Output:

```text
Loading sequences from mysequences.fasta...
✓ Loaded 5 sequences (20 bp)
Calculating hamming distances...
✓ Distance matrix computed
Identifying unique haplotypes...
✓ Found 3 unique haplotypes
Building MJN network...
✓ Network constructed

Network Statistics:
  Nodes: 3
  Edges: 2
  Median vectors: 0

✓ Network saved as GRAPHML
```

### Step 3: Visualize the Network

Create a visualization:

```bash
pypopart visualize mynetwork.graphml -o mynetwork.png
```

Or create an interactive HTML visualization:

```bash
pypopart visualize mynetwork.graphml -o mynetwork.html --interactive
```

**That's it!** You've created your first haplotype network.

## Exploring Further

### Try Different Algorithms

PyPopART offers six network construction algorithms:

```bash
# Minimum Spanning Tree (simplest tree)
pypopart network mysequences.fasta -a mst -o network_mst.graphml

# Minimum Spanning Network (shows alternatives)
pypopart network mysequences.fasta -a msn -o network_msn.graphml

# Statistical Parsimony (TCS)
pypopart network mysequences.fasta -a tcs -o network_tcs.graphml

# Median-Joining (default, infers ancestors)
pypopart network mysequences.fasta -a mjn -o network_mjn.graphml

# Parsimony Network (consensus from multiple trees)
pypopart network mysequences.fasta -a pn -o network_pn.graphml

# Tight Span Walker (metric-preserving, for small datasets)
pypopart network mysequences.fasta -a tsw -o network_tsw.graphml
```

**Algorithm Comparison:**

- **MST**: Fastest, simplest tree
- **MSN**: Adds alternative equal-distance connections
- **TCS**: Best for closely related sequences
- **MJN**: Infers ancestral sequences, handles reticulation
- **PN**: Consensus approach, captures uncertainty
- **TSW**: Most accurate metric preservation (slower)

### Use Different Distance Metrics

```bash
# Hamming distance (default, fastest)
pypopart network mysequences.fasta -d hamming -o network.graphml

# Jukes-Cantor correction
pypopart network mysequences.fasta -d jc -o network.graphml

# Kimura 2-parameter
pypopart network mysequences.fasta -d k2p -o network.graphml

# Tamura-Nei (most complex)
pypopart network mysequences.fasta -d tamura_nei -o network.graphml
```

### Analyze Your Network

Get comprehensive statistics:

```bash
pypopart analyze mynetwork.graphml --stats
```

Output:

```text
Loading network from mynetwork.graphml...
✓ Loaded network with 3 nodes

=== Network Statistics ===
nodes: 3
edges: 2
diameter: 2
avg_degree: 1.333
clustering_coefficient: 0.0000
reticulation_index: 0.0000
```

### Customize Visualizations

```bash
# Different layout algorithms
pypopart visualize mynetwork.graphml -o net.png --layout circular
pypopart visualize mynetwork.graphml -o net.png --layout radial

# Show node labels
pypopart visualize mynetwork.graphml -o net.png --show-labels

# Custom size
pypopart visualize mynetwork.graphml -o net.png --width 1200 --height 900

# Save as PDF (publication-ready)
pypopart visualize mynetwork.graphml -o net.pdf
```

## Using the Python API

For more control, use Python directly:

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

## Common Workflows

### Workflow 1: Basic Analysis

```bash
# 1. Load and inspect data
pypopart load mysequences.fasta

# 2. Build network
pypopart network mysequences.fasta -o network.graphml

# 3. Analyze
pypopart analyze network.graphml --stats --topology

# 4. Visualize
pypopart visualize network.graphml -o network.png
```

### Workflow 2: Algorithm Comparison

```bash
# Try all algorithms
for algo in mst msn tcs mjn pn tsw; do
    pypopart network mysequences.fasta -a $algo -o network_${algo}.graphml
    pypopart visualize network_${algo}.graphml -o network_${algo}.png
done
```

### Workflow 3: GUI-based Analysis

For an interactive workflow:

```bash
# Launch GUI
pypopart-gui

# Then use the web interface to:
# 1. Upload sequences and optional metadata
# 2. Try different algorithms with real-time parameter adjustment
# 3. Visualize with multiple layout options
# 4. Drag nodes to reposition
# 5. Export in various formats
```

### Workflow 4: With Metadata

If you have population or location data:

```bash
# Create metadata file (metadata.csv)
cat > metadata.csv << 'EOF'
id,population,location
Sample1,PopA,Site1
Sample2,PopA,Site1
Sample3,PopB,Site2
Sample4,PopB,Site2
Sample5,PopC,Site3
EOF

# Load with metadata
pypopart load mysequences.fasta -m metadata.csv

# Visualize colored by population
pypopart visualize network.graphml -o network.png --color-by population
```

## Next Steps

Now that you've created your first network, learn more:

- **[CLI Guide](guide/cli.md)**: Complete command-line reference
- **[Algorithms](guide/algorithms.md)**: Choosing the right algorithm
- **[Tutorials](tutorials/basic_workflow.md)**: Detailed examples
- **[API Reference](api/core/sequence.md)**: Python API documentation

## Getting Help

- Run `pypopart --help` for command help
- Run `pypopart COMMAND --help` for specific command help
- See [FAQ](faq.md) for common questions
- Check [GitHub Issues](https://github.com/adamtaranto/pypopart/issues) for known problems

## Tips

!!! tip "Performance"
For large datasets (>1000 sequences), consider:

- Using MST or MSN instead of MJN for speed
- Using Hamming distance instead of more complex models
- Increasing epsilon parameter for MJN to reduce complexity

!!! tip "Visualization"

- PNG for presentations and documents
- PDF for publications (vector format, scalable)
- SVG for web use
- HTML for interactive exploration

!!! tip "File Formats"

- GraphML preserves all network attributes (recommended)
- JSON for web applications
- NEXUS for PopART compatibility

Happy network building! 🌳
