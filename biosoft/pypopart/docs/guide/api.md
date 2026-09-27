# Python API Guide

Use PyPopART programmatically in your Python scripts and notebooks.

## Basic Workflow

```python
from pypopart import Alignment, Network
from pypopart.algorithms import MSTAlgorithm, MJNAlgorithm
from pypopart.visualization import StaticPlot, InteractivePlot

# 1. Load sequence data
alignment = Alignment.from_fasta('sequences.fasta')

# 2. Build network
algorithm = MSTAlgorithm()
network = algorithm.build_network(alignment)

# 3. Visualize
plot = StaticPlot(network)
plot.save('network.png')
```

## Loading Data

### From Files

```python
from pypopart import Alignment

# FASTA format
alignment = Alignment.from_fasta('sequences.fasta')

# NEXUS format (with metadata)
alignment = Alignment.from_nexus('sequences.nex')

# PHYLIP format
alignment = Alignment.from_phylip('sequences.phy')

# GenBank format
alignment = Alignment.from_genbank('sequences.gb')
```

### From Strings

```python
from pypopart import Alignment

fasta_string = """
>Seq1
ATCGATCGATCG
>Seq2
ATCGATCGATCG
>Seq3
ATCGATTGATCG
"""

alignment = Alignment.from_string(fasta_string, format='fasta')
```

### From BioPython

```python
from Bio import AlignIO
from pypopart import Alignment

# Load with BioPython
bio_alignment = AlignIO.read('sequences.fasta', 'fasta')

# Convert to PyPopART
alignment = Alignment.from_biopython(bio_alignment)
```

## Building Networks

### Algorithm Selection

```python
from pypopart.algorithms import (
    MSTAlgorithm,
    MSNAlgorithm,
    TCSAlgorithm,
    MJNAlgorithm,
    ParsimonyNetAlgorithm,
    TSWAlgorithm,
)

# Minimum Spanning Tree
mst = MSTAlgorithm()
network = mst.build_network(alignment)

# Median-Joining Network
mjn = MJNAlgorithm()
network = mjn.build_network(alignment)

# TCS with custom epsilon
tcs = TCSAlgorithm(epsilon=0.99)
network = tcs.build_network(alignment)
```

### Distance Metrics

```python
from pypopart.core.distance import DistanceCalculator

# Create calculator with specific metric
calc = DistanceCalculator(metric='k2p')  # Kimura 2-parameter
distances = calc.calculate(alignment)

# Available metrics: 'hamming', 'jukes-cantor', 'k2p', 'tamura-nei'

# Use custom distance matrix
algorithm = MSTAlgorithm(distance_matrix=distances)
network = algorithm.build_network(alignment)
```

## Working with Networks

### Network Properties

```python
# Basic properties
print(f'Number of nodes: {network.number_of_nodes()}')
print(f'Number of edges: {network.number_of_edges()}')
print(f'Network density: {network.density()}')

# Get nodes and edges
nodes = list(network.nodes())
edges = list(network.edges())

# Node attributes
for node in network.nodes():
    haplotype = network.nodes[node]['haplotype']
    frequency = network.nodes[node]['frequency']
    print(f'{node}: frequency={frequency}')

# Edge attributes
for u, v in network.edges():
    weight = network[u][v]['weight']
    print(f'{u} -> {v}: distance={weight}')
```

### Network Statistics

```python
from pypopart.stats import NetworkStatistics, TopologyAnalysis

# Calculate basic statistics
stats = NetworkStatistics(network)
print(f'Diameter: {stats.diameter()}')
print(f'Average path length: {stats.average_path_length()}')
print(f'Clustering coefficient: {stats.clustering_coefficient()}')

# Topology analysis
topology = TopologyAnalysis(network)
hubs = topology.identify_hubs()
bridges = topology.identify_bridges()
star_patterns = topology.detect_star_patterns()
```

### Population Genetics

```python
from pypopart.stats import PopulationGenetics

popgen = PopulationGenetics(alignment)

# Diversity measures
print(f'Nucleotide diversity: {popgen.nucleotide_diversity()}')
print(f'Haplotype diversity: {popgen.haplotype_diversity()}')

# Neutrality tests
print(f"Tajima's D: {popgen.tajimas_d()}")
print(f"Fu's Fs: {popgen.fus_fs()}")

# Population differentiation (requires population metadata)
fst = popgen.calculate_fst(population_column='Population')
print(f'FST: {fst}')
```

## Visualization

### Static Plots

```python
from pypopart.visualization import StaticPlot

# Basic plot
plot = StaticPlot(network)
plot.save('network.png')

# Customized plot
plot = StaticPlot(
    network,
    layout='spring',  # or 'circular', 'kamada-kawai'
    node_size=500,
    edge_width=2.0,
    figsize=(12, 10),
    dpi=300,
)

# Color by metadata
plot.color_by_attribute('Population')
plot.save('colored_network.png', format='pdf')
```

### Interactive Plots

```python
from pypopart.visualization import InteractivePlot

# Create interactive HTML plot
plot = InteractivePlot(network)
plot.save('network.html')

# With custom styling
plot = InteractivePlot(
    network,
    layout='spring',
    node_size_by='frequency',
    color_by='Population',
    show_labels=True,
)
plot.save('interactive_network.html')
```

### Layout Algorithms

```python
from pypopart.layout import LayoutAlgorithm

# Use specific layout
layout = LayoutAlgorithm.spring(network, k=1.0, iterations=50)
plot = StaticPlot(network, positions=layout)

# Available layouts
layouts = [
    'spring',  # Force-directed (default)
    'circular',  # Circular arrangement
    'kamada-kawai',  # Energy minimization
    'spectral',  # Eigenvalue-based
    'random',  # Random positions
]
```

## Exporting Results

### Save Networks

```python
# Export in various formats
network.save('network.gml')  # GML format
network.save('network.graphml')  # GraphML format
network.save('network.json')  # JSON format
network.save('network.nexus')  # NEXUS format

# With metadata
network.save('network.nexus', include_traits=True)
```

### Export Distance Matrix

```python
from pypopart.core.distance import DistanceCalculator

calc = DistanceCalculator(metric='k2p')
distances = calc.calculate(alignment)

# Save as CSV
distances.to_csv('distances.csv')

# Save as NumPy array
import numpy as np

np.save('distances.npy', distances.matrix)
```

### Export Statistics

```python
from pypopart.stats import NetworkStatistics
import pandas as pd

stats = NetworkStatistics(network)
results = {
    'diameter': stats.diameter(),
    'avg_path_length': stats.average_path_length(),
    'clustering': stats.clustering_coefficient(),
    'num_nodes': network.number_of_nodes(),
    'num_edges': network.number_of_edges(),
}

# Save as DataFrame
df = pd.DataFrame([results])
df.to_csv('statistics.csv', index=False)
```

## Advanced Usage

### Custom Algorithms

```python
from pypopart.algorithms.base import BaseAlgorithm


class CustomAlgorithm(BaseAlgorithm):
    def build_network(self, alignment):
        # Your implementation
        network = self.create_empty_network()
        # ... add nodes and edges
        return network


# Use your algorithm
algorithm = CustomAlgorithm()
network = algorithm.build_network(alignment)
```

### Batch Processing

```python
from pathlib import Path
from pypopart import Alignment
from pypopart.algorithms import MSTAlgorithm

# Process multiple files
algorithm = MSTAlgorithm()

for fasta_file in Path('data').glob('*.fasta'):
    alignment = Alignment.from_fasta(fasta_file)
    network = algorithm.build_network(alignment)
    network.save(f'networks/{fasta_file.stem}.gml')
```

### Integration with NetworkX

```python
import networkx as nx

# PyPopART networks are NetworkX graphs
# Use any NetworkX function

# Graph algorithms
shortest_paths = nx.shortest_path(network)
betweenness = nx.betweenness_centrality(network)
communities = nx.community.louvain_communities(network)

# Export to other NetworkX formats
nx.write_gexf(network, 'network.gexf')
```

## Error Handling

```python
from pypopart.exceptions import (
    AlignmentError,
    NetworkError,
    DistanceError,
)

try:
    alignment = Alignment.from_fasta('sequences.fasta')
    network = algorithm.build_network(alignment)
except AlignmentError as e:
    print(f'Alignment error: {e}')
except NetworkError as e:
    print(f'Network construction error: {e}')
except FileNotFoundError:
    print('Input file not found')
```

## Next Steps

- [Loading Data Guide](loading_data.md): Detailed data loading instructions
- [Algorithm Guide](algorithms.md): Choose the right algorithm
- [Visualization Guide](visualization.md): Create publication-quality figures
- [API Reference](../api/core/sequence.md): Complete API documentation
