# Visualization Tutorial

Learn advanced visualization techniques for creating publication-quality figures.

## Overview

This tutorial covers:

- Static plots for publications
- Interactive HTML plots
- Color schemes and styling
- Layouts and customization
- Multi-panel figures

## Basic Visualization

```python
from pypopart import Alignment
from pypopart.algorithms import MJNAlgorithm
from pypopart.visualization import StaticPlot, InteractivePlot

# Load and build network
alignment = Alignment.from_fasta('sequences.fasta')
network = MJNAlgorithm().build_network(alignment)

# Simple static plot
plot = StaticPlot(network)
plot.save('network.png')

# Simple interactive plot
interactive = InteractivePlot(network)
interactive.save('network.html')
```

## Publication-Quality Figure

```python
import matplotlib.pyplot as plt

# Set publication style
plt.rcParams.update(
    {
        'font.size': 12,
        'font.family': 'Arial',
        'axes.linewidth': 1.5,
        'figure.dpi': 300,
    }
)

# Create plot
plot = StaticPlot(
    network,
    figsize=(8, 8),
    layout='spring',
    node_size=400,
    edge_width=1.5,
    show_labels=True,
    font_size=10,
)

# Save in multiple formats
plot.save('figure1.png', dpi=300)  # For presentations
plot.save('figure1.pdf')  # For publication
plot.save('figure1.svg')  # For editing
```

## Coloring by Metadata

```python
# Color nodes by population
plot = StaticPlot(network)
plot.color_by_attribute('Population')
plot.add_legend(title='Population')
plot.save('colored_by_population.png')

# Custom color scheme
color_map = {
    'PopA': '#1f78b4',
    'PopB': '#33a02c',
    'PopC': '#e31a1c',
}
plot.color_by_attribute('Population', color_map=color_map)
plot.save('custom_colors.png')

# Continuous values (e.g., year)
plot.color_by_attribute('Year', cmap='viridis')
plot.add_colorbar(label='Year')
plot.save('temporal.png')
```

## Node Sizing

```python
# Size by frequency
plot = StaticPlot(network)
plot.size_by_frequency(min_size=100, max_size=1000)
plot.save('sized_by_frequency.png')

# Size by custom attribute
plot.size_by_attribute('SampleSize', min_size=50, max_size=800)
plot.save('sized_by_samples.png')
```

## Layout Algorithms

```python
layouts = ['spring', 'circular', 'kamada-kawai', 'spectral']

fig, axes = plt.subplots(2, 2, figsize=(16, 16))
axes = axes.flatten()

for ax, layout in zip(axes, layouts):
    plot = StaticPlot(network, ax=ax, layout=layout)
    ax.set_title(f'{layout.title()} Layout')

plt.tight_layout()
plt.savefig('layout_comparison.png', dpi=300)
```

## Interactive Visualization

```python
# Rich interactive plot
plot = InteractivePlot(
    network,
    layout='spring',
    width=1200,
    height=800,
    title='Interactive Haplotype Network',
    node_size_by='frequency',
    color_by='Population',
    show_edge_labels=True,
    hover_data=['Population', 'Location', 'Frequency'],
)

plot.save('interactive_network.html')
```

## Multi-Panel Figures

```python
from pypopart.algorithms import MSTAlgorithm, MSNAlgorithm, MJNAlgorithm

# Build multiple networks
mst_net = MSTAlgorithm().build_network(alignment)
msn_net = MSNAlgorithm().build_network(alignment)
mjn_net = MJNAlgorithm().build_network(alignment)

# Create multi-panel figure
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

networks = [mst_net, msn_net, mjn_net]
titles = ['MST', 'MSN', 'MJN']

for ax, net, title in zip(axes, networks, titles):
    plot = StaticPlot(net, ax=ax, layout='spring')
    plot.color_by_attribute('Population')
    ax.set_title(title, fontsize=16, fontweight='bold')

plt.tight_layout()
plt.savefig('multi_panel.png', dpi=300)
```

## Highlighting Specific Nodes

```python
# Highlight hub haplotypes
from pypopart.stats import TopologyAnalysis

topology = TopologyAnalysis(network)
hubs = topology.identify_hubs()

plot = StaticPlot(network)
plot.highlight_nodes(
    hubs, color='#ff0000', size=800, edge_color='#ff0000', edge_width=3
)
plot.save('highlighted_hubs.png')
```

## Styling Median Vectors (MJN)

```python
# Different style for observed vs inferred nodes
observed = [n for n in network.nodes() if network.nodes[n].get('type') == 'observed']
median = [n for n in network.nodes() if network.nodes[n].get('type') == 'median']

plot = StaticPlot(network)
plot.style_nodes(observed, shape='o', color='#1f78b4', size=500)
plot.style_nodes(median, shape='s', color='#ff7f00', size=200)
plot.add_legend(
    handles={
        'Observed': {'color': '#1f78b4', 'marker': 'o'},
        'Inferred': {'color': '#ff7f00', 'marker': 's'},
    }
)
plot.save('styled_medians.png')
```

## Adding Annotations

```python
plot = StaticPlot(network, layout='spring')

# Annotate specific nodes
plot.annotate_node('H1', 'Most common', fontsize=12, color='red')
plot.annotate_node('H5', 'Ancestral?', fontsize=12, color='blue')

# Add scale bar
plot.add_scale_bar(mutations_per_unit=1, label='1 mutation')

# Add text box with information
plot.add_textbox(
    'MJN Algorithm\nK2P Distance\nn=50', position=(0.05, 0.95), fontsize=10
)

plot.save('annotated_network.png')
```

## Colorblind-Friendly Palettes

```python
from pypopart.visualization.colors import colorblind_safe

plot = StaticPlot(network)
plot.color_by_attribute('Population', colors=colorblind_safe)
plot.save('colorblind_safe.png')
```

## Export for Other Tools

```python
# For Gephi
network.save('network.gexf')

# For Cytoscape
network.save('network.graphml')

# For R/igraph
network.save('network.gml')

# With metadata
network.save('network_with_traits.nexus', include_traits=True)
```

## Complete Visualization Script

```python
from pypopart import Alignment
from pypopart.algorithms import MJNAlgorithm
from pypopart.visualization import StaticPlot, InteractivePlot
from pypopart.stats import TopologyAnalysis
import matplotlib.pyplot as plt

# Load data and build network
alignment = Alignment.from_fasta('sequences.fasta')
network = MJNAlgorithm().build_network(alignment)

# Identify hubs
hubs = TopologyAnalysis(network).identify_hubs()

# Create static plot
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'

plot = StaticPlot(
    network,
    figsize=(10, 10),
    layout='spring',
    node_size=400,
)

# Style by metadata
plot.color_by_attribute('Population')
plot.size_by_frequency(min_size=200, max_size=1000)
plot.highlight_nodes(hubs, color='#ff0000', size=900)

# Add annotations
plot.add_legend(title='Population')
plot.add_scale_bar(mutations_per_unit=1)

# Save in multiple formats
plot.save('network.png', dpi=300)
plot.save('network.pdf')

# Create interactive version
interactive = InteractivePlot(
    network,
    color_by='Population',
    node_size_by='frequency',
    hover_data=['Population', 'Frequency'],
)
interactive.save('network.html')

print('Visualization complete!')
```

## Tips for Great Figures

1. **Resolution**: Use DPI ≥ 300 for publications
2. **Format**: PDF/SVG for journals, PNG for presentations
3. **Colors**: Use colorblind-safe palettes
4. **Size**: Large enough to be readable when printed
5. **Labels**: Clear and informative
6. **Legend**: Always include for colored nodes
7. **Scale**: Add scale bar for mutation context

## Next Steps

- [Basic Workflow](basic_workflow.md): Complete analysis
- [Population Genetics](popgen.md): Advanced analyses
- [Visualization Guide](../guide/visualization.md): Full documentation
