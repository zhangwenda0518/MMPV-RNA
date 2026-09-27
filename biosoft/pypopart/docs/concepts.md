# Basic Concepts

Understanding the fundamental concepts behind PyPopART will help you make the most of the software.

## What are Haplotype Networks?

A **haplotype network** is a graphical representation showing the evolutionary relationships between different DNA sequences (haplotypes) in a population. Unlike phylogenetic trees, networks can show:

- **Multiple evolutionary paths** between sequences
- **Reticulation events** like recombination or homoplasy
- **Ancestral sequences** that may still exist in the population
- **Mutation steps** between related sequences

### Key Components

- **Nodes**: Represent observed haplotypes (unique DNA sequences)
- **Edges**: Connect haplotypes differing by mutations
- **Edge length**: Number of mutations between haplotypes
- **Node size**: Often proportional to haplotype frequency
- **Median vectors**: Inferred ancestral or unobserved haplotypes (in some algorithms)

## Haplotypes

A **haplotype** is a unique DNA sequence variant in your dataset. PyPopART:

1. Reads your sequence alignment
2. Identifies unique sequences
3. Counts how many times each appears
4. Uses this information to build the network

### Important Properties

- **Frequency**: How many individuals share this haplotype
- **Sequences**: The actual DNA/protein sequence
- **Metadata**: Associated information (location, population, traits)

## Distance Metrics

PyPopART calculates genetic distances using various evolutionary models:

### Hamming Distance

Simple count of differing positions. Best for closely related sequences.

### Jukes-Cantor

Corrects for multiple mutations at the same site. Assumes equal substitution rates.

### Kimura 2-Parameter (K2P)

Distinguishes between transitions and transversions. More realistic for DNA evolution.

### Tamura-Nei

Most sophisticated, accounts for different base frequencies and transition/transversion ratios.

## Network Algorithms

Different algorithms make different assumptions and are suited for different data types:

### MST (Minimum Spanning Tree)

- Simplest algorithm
- Always produces a tree (no reticulation)
- Connects all haplotypes with minimum total distance
- **Best for**: Initial exploration, small datasets

### MSN (Minimum Spanning Network)

- Extends MST by adding alternative connections
- Shows equally parsimonious paths
- More informative than MST
- **Best for**: Showing alternative evolutionary paths

### TCS (Statistical Parsimony)

- Based on statistical limits of parsimony
- Uses 95% confidence limit for connections
- May produce disconnected networks
- **Best for**: Within-species variation, recent divergence

### MJN (Median-Joining Network)

- Infers ancestral/unobserved haplotypes (median vectors)
- Most comprehensive but complex
- Can show reticulation events
- **Best for**: Complex evolutionary scenarios, larger datasets

### PN (Parsimony Network)

- Consensus approach using multiple MSTs
- Balances between MST simplicity and MJN complexity
- **Best for**: General purpose analysis

### TSW (Tight Span Walker)

- Metric-preserving network construction
- Preserves distance relationships
- **Best for**: When distance accuracy is critical

## Metadata and Traits

PyPopART supports associating metadata with sequences:

- **Populations**: Geographic or demographic groups
- **Sampling dates**: Temporal information
- **Phenotypes**: Traits or characteristics
- **Custom attributes**: Any categorical or numerical data

Metadata can be used for:

- Coloring nodes in visualizations
- Population genetics analyses
- Statistical comparisons
- Pattern identification

## Network Statistics

PyPopART calculates various network properties:

### Topology Metrics

- **Diameter**: Longest shortest path
- **Clustering coefficient**: Network interconnectedness
- **Centrality**: Important nodes in the network

### Population Genetics

- **Diversity indices**: Nucleotide and haplotype diversity
- **Tajima's D**: Test for neutral evolution
- **Fu's Fs**: Another neutrality test
- **FST**: Population differentiation

## Visualization

Networks can be visualized in multiple ways:

### Static Plots

- Publication-quality figures
- PNG, PDF, SVG formats
- Customizable colors, sizes, labels

### Interactive Plots

- HTML-based exploration
- Zoom, pan, hover for details
- Export functionality

### GUI Dashboard

- Real-time parameter adjustment
- Multiple layout algorithms
- Integrated analysis tools

## Workflow Overview

A typical PyPopART analysis:

1. **Load Data**: Import sequence alignment
2. **Calculate Distances**: Choose appropriate metric
3. **Build Network**: Select algorithm
4. **Analyze**: Compute statistics
5. **Visualize**: Create plots
6. **Export**: Save results

## Next Steps

- [Installation](installation.md): Get PyPopART set up
- [Quick Start](quickstart.md): Try your first analysis
- [Tutorials](tutorials/basic_workflow.md): Detailed walkthroughs
- [User Guide](guide/cli.md): Complete documentation
