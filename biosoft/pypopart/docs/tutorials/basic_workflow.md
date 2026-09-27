# Basic Workflow Tutorial

This tutorial walks through a complete PyPopART analysis from start to finish.

## Overview

We'll analyze a sample dataset of mitochondrial DNA sequences to:

1. Load sequence data
2. Build a haplotype network
3. Calculate statistics
4. Create visualizations
5. Export results

## Prerequisites

```bash
pip install pypopart
```

## Step 1: Prepare Your Data

Create a FASTA file with aligned sequences:

```fasta
>Sample1_PopA
ATCGATCGATCGATCGATCGTACGATCG
>Sample2_PopA
ATCGATCGATCGATCGATCGTACGATCG
>Sample3_PopA
ATCGATCGATCGATCGATCGTACGATCG
>Sample4_PopB
ATCGATCGATCGATTGATCGTACGATCG
>Sample5_PopB
ATCGATCGATCGATTGATCGTACGATCG
>Sample6_PopC
ATCGATCGATCGATTGATCGTACGTTCG
```

Save as `sequences.fasta`.

## Step 2: Load Data

```python
from pypopart import Alignment

# Load sequences
alignment = Alignment.from_fasta('sequences.fasta')

# Check alignment
print(f'Sequences: {len(alignment)}')
print(f'Length: {alignment.length}')
print(f'Unique haplotypes: {alignment.n_unique()}')
```

## Step 3: Build Network

```python
from pypopart.algorithms import MSTAlgorithm, MJNAlgorithm

# Start with simple MST
mst_algorithm = MSTAlgorithm()
mst_network = mst_algorithm.build_network(alignment)

print(
    f'MST - Nodes: {mst_network.number_of_nodes()}, Edges: {mst_network.number_of_edges()}'
)

# Try comprehensive MJN
mjn_algorithm = MJNAlgorithm()
mjn_network = mjn_algorithm.build_network(alignment)

print(
    f'MJN - Nodes: {mjn_network.number_of_nodes()}, Edges: {mjn_network.number_of_edges()}'
)
```

## Step 4: Calculate Statistics

```python
from pypopart.stats import NetworkStatistics, PopulationGenetics

# Network statistics
stats = NetworkStatistics(mst_network)
print(f'\nNetwork Statistics:')
print(f'  Diameter: {stats.diameter()}')
print(f'  Avg path length: {stats.average_path_length():.3f}')
print(f'  Clustering: {stats.clustering_coefficient():.3f}')

# Population genetics
popgen = PopulationGenetics(alignment)
print(f'\nPopulation Genetics:')
print(f'  Haplotype diversity: {popgen.haplotype_diversity():.4f}')
print(f'  Nucleotide diversity: {popgen.nucleotide_diversity():.4f}')
print(f"  Tajima's D: {popgen.tajimas_d():.4f}")
```

## Step 5: Visualize

```python
from pypopart.visualization import StaticPlot, InteractivePlot

# Static plot
static = StaticPlot(mst_network, figsize=(10, 10))
static.save('mst_network.png', dpi=300)

# Interactive HTML plot
interactive = InteractivePlot(mjn_network)
interactive.save('mjn_network.html')

print('\nPlots saved!')
```

## Step 6: Export Results

```python
# Save networks
mst_network.save('mst_network.gml')
mjn_network.save('mjn_network.nexus')

# Save statistics report
import json

report = {
    'network': {
        'nodes': stats.number_of_nodes(),
        'edges': stats.number_of_edges(),
        'diameter': stats.diameter(),
    },
    'diversity': {
        'haplotype': popgen.haplotype_diversity(),
        'nucleotide': popgen.nucleotide_diversity(),
    },
}

with open('analysis_report.json', 'w') as f:
    json.dump(report, f, indent=2)

print('Results exported!')
```

## Complete Script

Here's the full workflow in one script:

```python
from pypopart import Alignment
from pypopart.algorithms import MSTAlgorithm, MJNAlgorithm
from pypopart.stats import NetworkStatistics, PopulationGenetics
from pypopart.visualization import StaticPlot, InteractivePlot
import json

# Load data
alignment = Alignment.from_fasta('sequences.fasta')
print(f'Loaded {len(alignment)} sequences')

# Build networks
mst_network = MSTAlgorithm().build_network(alignment)
mjn_network = MJNAlgorithm().build_network(alignment)

# Calculate statistics
stats = NetworkStatistics(mst_network)
popgen = PopulationGenetics(alignment)

# Create visualizations
StaticPlot(mst_network).save('mst_network.png')
InteractivePlot(mjn_network).save('mjn_network.html')

# Export results
mst_network.save('mst_network.gml')
report = {
    'nodes': stats.number_of_nodes(),
    'haplotype_diversity': popgen.haplotype_diversity(),
    'nucleotide_diversity': popgen.nucleotide_diversity(),
}
with open('report.json', 'w') as f:
    json.dump(report, f, indent=2)

print('Analysis complete!')
```

## CLI Equivalent

The same workflow using command-line tools:

```bash
# Build network
pypopart network sequences.fasta -a MST -o mst_network --plot

# Calculate statistics
pypopart stats mst_network.gml --topology --popgen -o stats.txt

# Create visualization
pypopart plot mst_network.gml -o network.png --width 2000 --height 2000
```

## Next Steps

- [Algorithm Comparison](algorithm_comparison.md): Compare different algorithms
- [Visualization Tutorial](visualization.md): Advanced plotting
- [Population Genetics](popgen.md): Detailed population analysis
- [User Guide](../guide/cli.md): Complete documentation
