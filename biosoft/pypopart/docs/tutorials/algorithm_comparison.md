# Algorithm Comparison Tutorial

Learn how to compare different network algorithms and choose the best one for your data.

## Overview

We'll compare MST, MSN, TCS, and MJN algorithms on the same dataset to understand their differences.

## Setup

```python
from pypopart import Alignment
from pypopart.algorithms import MSTAlgorithm, MSNAlgorithm, TCSAlgorithm, MJNAlgorithm
from pypopart.visualization import StaticPlot
from pypopart.stats import NetworkStatistics
import matplotlib.pyplot as plt

# Load data
alignment = Alignment.from_fasta('sequences.fasta')
```

## Build Networks with All Algorithms

```python
algorithms = {
    'MST': MSTAlgorithm(),
    'MSN': MSNAlgorithm(),
    'TCS': TCSAlgorithm(epsilon=0.95),
    'MJN': MJNAlgorithm(),
}

networks = {}
for name, algorithm in algorithms.items():
    print(f'Building {name} network...')
    networks[name] = algorithm.build_network(alignment)
    n_nodes = networks[name].number_of_nodes()
    n_edges = networks[name].number_of_edges()
    print(f'  {name}: {n_nodes} nodes, {n_edges} edges')
```

## Compare Network Properties

```python
import pandas as pd

# Collect statistics
results = []
for name, network in networks.items():
    stats = NetworkStatistics(network)
    results.append(
        {
            'Algorithm': name,
            'Nodes': stats.number_of_nodes(),
            'Edges': stats.number_of_edges(),
            'Diameter': stats.diameter(),
            'Avg Path Length': stats.average_path_length(),
            'Clustering': stats.clustering_coefficient(),
        }
    )

# Create comparison table
df = pd.DataFrame(results)
print('\nNetwork Comparison:')
print(df.to_string(index=False))

# Save table
df.to_csv('algorithm_comparison.csv', index=False)
```

## Visualize Side-by-Side

```python
fig, axes = plt.subplots(2, 2, figsize=(16, 16))
axes = axes.flatten()

for idx, (name, network) in enumerate(networks.items()):
    plot = StaticPlot(network, ax=axes[idx], layout='spring')
    axes[idx].set_title(f'{name} Network', fontsize=16, fontweight='bold')

plt.tight_layout()
plt.savefig('algorithm_comparison.png', dpi=300)
print('\nComparison figure saved!')
```

## Analyze Differences

```python
# Compare node sets
mst_nodes = set(networks['MST'].nodes())
mjn_nodes = set(networks['MJN'].nodes())

# Median vectors (inferred nodes) in MJN
inferred = mjn_nodes - mst_nodes
print(f'\nMJN inferred {len(inferred)} median vectors')

# Compare connectivity
for name, network in networks.items():
    density = NetworkStatistics(network).density()
    print(f'{name} density: {density:.3f}')
```

## Decision Guide

```python
def recommend_algorithm(alignment):
    """Suggest best algorithm based on data characteristics."""
    n_seqs = len(alignment)
    diversity = alignment.pairwise_diversity()

    if n_seqs < 20:
        return 'MST - Small dataset, start simple'
    elif diversity < 0.01:
        return 'TCS - Low diversity, within-species'
    elif n_seqs < 100:
        return 'MJN - Medium dataset, comprehensive analysis'
    else:
        return 'MSN - Large dataset, balance speed and information'


recommendation = recommend_algorithm(alignment)
print(f'\nRecommendation: {recommendation}')
```

## Computational Performance

```python
import time

times = {}
for name, algorithm in algorithms.items():
    start = time.time()
    algorithm.build_network(alignment)
    elapsed = time.time() - start
    times[name] = elapsed
    print(f'{name}: {elapsed:.3f} seconds')

# Plot timing
plt.figure(figsize=(8, 6))
plt.bar(times.keys(), times.values())
plt.xlabel('Algorithm')
plt.ylabel('Time (seconds)')
plt.title('Computational Performance')
plt.savefig('algorithm_timing.png', dpi=300)
```

## Complete Comparison Script

```python
from pypopart import Alignment
from pypopart.algorithms import MSTAlgorithm, MSNAlgorithm, TCSAlgorithm, MJNAlgorithm
from pypopart.visualization import StaticPlot
from pypopart.stats import NetworkStatistics
import matplotlib.pyplot as plt
import pandas as pd
import time

# Load data
alignment = Alignment.from_fasta('sequences.fasta')

# Define algorithms
algorithms = {
    'MST': MSTAlgorithm(),
    'MSN': MSNAlgorithm(),
    'TCS': TCSAlgorithm(),
    'MJN': MJNAlgorithm(),
}

# Build and time networks
networks = {}
times = {}
results = []

for name, algorithm in algorithms.items():
    start = time.time()
    networks[name] = algorithm.build_network(alignment)
    times[name] = time.time() - start

    stats = NetworkStatistics(networks[name])
    results.append(
        {
            'Algorithm': name,
            'Nodes': stats.number_of_nodes(),
            'Edges': stats.number_of_edges(),
            'Time (s)': f'{times[name]:.3f}',
        }
    )

# Print comparison
df = pd.DataFrame(results)
print(df.to_string(index=False))

# Create visualizations
fig, axes = plt.subplots(2, 2, figsize=(16, 16))
axes = axes.flatten()

for idx, (name, network) in enumerate(networks.items()):
    StaticPlot(network, ax=axes[idx], layout='spring')
    axes[idx].set_title(f'{name} - {times[name]:.2f}s')

plt.tight_layout()
plt.savefig('algorithm_comparison.png', dpi=300)
print('Comparison complete!')
```

## When to Use Each Algorithm

### MST

- ✅ Quick exploration
- ✅ Small datasets
- ✅ Simple relationships
- ❌ Ignores alternative paths

### MSN

- ✅ Alternative connections
- ✅ Medium datasets
- ✅ Ambiguous relationships
- ❌ Less complete than MJN

### TCS

- ✅ Within-species data
- ✅ Statistical justification
- ✅ Population studies
- ❌ May be disconnected

### MJN

- ✅ Comprehensive analysis
- ✅ Ancestral inference
- ✅ Publication figures
- ❌ Slower computation
- ❌ Complex interpretation

## Next Steps

- [Basic Workflow](basic_workflow.md): Complete analysis tutorial
- [Visualization Tutorial](visualization.md): Advanced plotting
- [Algorithm Guide](../guide/algorithms.md): Detailed documentation
