# Visualization Guide

PyPopART provides multiple visualization options for creating publication-quality figures and interactive explorations.

## Visualization Types

### Static Plots

- **Format**: PNG, PDF, SVG
- **Use**: Publications, presentations, reports
- **Library**: Matplotlib
- **Customization**: High

### Interactive Plots

- **Format**: HTML
- **Use**: Exploration, web embedding
- **Library**: Plotly
- **Features**: Zoom, pan, hover, click

### GUI Dashboard

- **Format**: Web application
- **Use**: Interactive analysis
- **Library**: Dash + Cytoscape
- **Features**: Real-time updates, full control

## Static Visualization

### Basic Plot

```python
from pypopart.visualization import StaticPlot

plot = StaticPlot(network)
plot.save('network.png')
```

CLI:

```bash
pypopart plot network.gml -o figure.png
```

### Customization

```python
plot = StaticPlot(
    network,
    figsize=(12, 10),  # Figure size in inches
    dpi=300,  # Resolution
    layout='spring',  # Layout algorithm
    node_size=500,  # Base node size
    node_color='#1f78b4',  # Default node color
    edge_width=2.0,  # Edge width
    edge_color='#888888',  # Edge color
    font_size=10,  # Label font size
    show_labels=True,  # Show node labels
    title='Haplotype Network',
)
plot.save('custom_network.png')
```

### Output Formats

```python
# PNG (raster)
plot.save('network.png', format='png', dpi=300)

# PDF (vector, best for publications)
plot.save('network.pdf', format='pdf')

# SVG (vector, editable)
plot.save('network.svg', format='svg')

# EPS (vector, legacy journals)
plot.save('network.eps', format='eps')
```

### Layout Algorithms

```python
# Spring layout (force-directed)
plot = StaticPlot(network, layout='spring')

# Circular layout
plot = StaticPlot(network, layout='circular')

# Kamada-Kawai (energy minimization)
plot = StaticPlot(network, layout='kamada-kawai')

# Spectral layout (eigenvalue-based)
plot = StaticPlot(network, layout='spectral')

# Custom positions
positions = {node: (x, y) for node, x, y in ...}
plot = StaticPlot(network, positions=positions)
```

### Coloring by Metadata

```python
# Color by population
plot = StaticPlot(network)
plot.color_by_attribute('Population')
plot.save('colored_network.png')

# Custom color map
color_map = {'PopA': '#ff0000', 'PopB': '#0000ff'}
plot.color_by_attribute('Population', color_map=color_map)

# Continuous values
plot.color_by_attribute('Year', cmap='viridis')
```

### Node Sizing

```python
# Size by frequency
plot = StaticPlot(network, size_by='frequency', node_size=50)

# Size by custom attribute
plot.size_by_attribute('SampleSize', min_size=100, max_size=1000)

# Fixed sizes
plot = StaticPlot(network, node_size=300)
```

### Edge Styling

```python
# Edge width by distance
plot = StaticPlot(network)
plot.edge_width_by_distance(min_width=0.5, max_width=5.0)

# Edge color by type
plot.set_edge_colors({'observed': '#000000', 'inferred': '#888888'})

# Dashed edges for inferred connections
plot.style_inferred_edges(style='dashed', alpha=0.5)
```

## Interactive Visualization

### Basic Interactive Plot

```python
from pypopart.visualization import InteractivePlot

plot = InteractivePlot(network)
plot.save('network.html')
```

Open `network.html` in a browser for interactive exploration.

### Features

- **Hover**: View node/edge details
- **Click**: Select nodes
- **Zoom**: Scroll to zoom
- **Pan**: Drag to pan
- **Export**: Save as PNG from browser

### Customization

```python
plot = InteractivePlot(
    network,
    layout='spring',
    title='Interactive Haplotype Network',
    width=1200,
    height=800,
    node_size_by='frequency',
    color_by='Population',
    show_edge_labels=True,
    hover_data=['Population', 'Location', 'Year'],
)
plot.save('interactive_network.html')
```

### Embedding in Web Pages

```html
<iframe src="network.html" width="100%" height="600px"></iframe>
```

Or use the HTML directly:

```python
html_string = plot.to_html()
# Insert into your web application
```

## GUI Dashboard

### Launch GUI

```bash
pypopart-gui
```

Navigate to `http://localhost:8050`

### Features

1. **File Upload**: Drag-and-drop sequence files
2. **Algorithm Selection**: Choose and configure algorithm
3. **Network Computation**: Real-time building
4. **Layout Adjustment**: Interactive layout control
5. **Visual Customization**: Colors, sizes, labels
6. **Metadata Integration**: Color by traits
7. **Export**: Save networks and figures
8. **Statistics**: View network properties

### Programmatic GUI

```python
from pypopart.gui import NetworkDashboard

dashboard = NetworkDashboard(network=network, alignment=alignment, port=8050)
dashboard.run()
```

## Advanced Visualization

### Multiple Networks Side-by-Side

```python
import matplotlib.pyplot as plt
from pypopart.visualization import StaticPlot

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

algorithms = ['MST', 'MSN', 'MJN']
for ax, name, net in zip(axes, algorithms, networks):
    plot = StaticPlot(net, ax=ax)
    ax.set_title(name)

plt.tight_layout()
plt.savefig('comparison.png', dpi=300)
```

### Subnetwork Visualization

```python
# Extract subnetwork (e.g., one population)
subnetwork = network.subgraph(
    [n for n in network.nodes() if network.nodes[n].get('Population') == 'PopA']
)

plot = StaticPlot(subnetwork)
plot.save('subnetwork.png')
```

### Time Series Animation

```python
from pypopart.visualization import NetworkAnimation

# Create animation showing network growth over time
anim = NetworkAnimation(networks_by_year)
anim.save('network_evolution.gif')
```

### Geographic Overlay

```python
import geopandas as gpd
from pypopart.visualization import GeographicPlot

# Plot network on map (requires coordinates in metadata)
geoplot = GeographicPlot(network, basemap='world')
geoplot.save('geographic_network.png')
```

## Styling Tips

### Publication-Quality Figures

```python
import matplotlib.pyplot as plt

plt.rcParams.update(
    {
        'font.size': 12,
        'font.family': 'Arial',
        'axes.linewidth': 1.5,
        'xtick.major.width': 1.5,
        'ytick.major.width': 1.5,
    }
)

plot = StaticPlot(
    network,
    figsize=(8, 8),
    dpi=300,
    node_size=400,
    edge_width=1.5,
    font_size=11,
)
plot.save('publication_figure.pdf')
```

### Color Schemes

```python
# Colorblind-friendly palettes
from pypopart.visualization.colors import colorblind_safe

plot.color_by_attribute('Population', colors=colorblind_safe)

# Sequential for continuous data
plot.color_by_attribute('Year', cmap='Blues')

# Diverging for deviation from mean
plot.color_by_attribute('Diversity', cmap='RdBu')
```

### Legends and Annotations

```python
# Add legend
plot.add_legend(title='Population', loc='upper right')

# Add scale bar
plot.add_scale_bar(mutations_per_unit=1)

# Annotate specific nodes
plot.annotate_node('Haplotype_1', 'Ancestral', fontsize=12, color='red')

# Add text box
plot.add_textbox('MST Algorithm\nK2P Distance', position=(0.05, 0.95))
```

## Export for External Tools

### Export for Gephi

```python
# Save as GEXF
network.save('network.gexf')
# Open in Gephi for advanced visualization
```

### Export for Cytoscape

```python
# Save as GraphML
network.save('network.graphml')
# Import into Cytoscape
```

### Export for R/igraph

```python
# Save as GML
network.save('network.gml')
```

R code:

```r
library(igraph)
network <- read_graph("network.gml", format="gml")
plot(network)
```

## Performance Optimization

### Large Networks

```python
# Reduce complexity for large networks
plot = StaticPlot(
    network,
    show_labels=False,  # Disable labels
    edge_width=0.5,  # Thin edges
    simplify_edges=True,  # Merge parallel edges
)
```

### Batch Export

```python
# Export multiple visualizations efficiently
layouts = ['spring', 'circular', 'kamada-kawai']

for layout in layouts:
    plot = StaticPlot(network, layout=layout)
    plot.save(f'network_{layout}.png')
```

## Common Customizations

### Highlight Specific Nodes

```python
highlight_nodes = ['H1', 'H2', 'H5']
plot = StaticPlot(network)
plot.highlight_nodes(
    highlight_nodes, color='#ff0000', size=800, edge_color='#ff0000', edge_width=3
)
```

### Two-Color Network

```python
# Example: Mark ingroup vs outgroup
ingroup = [n for n in network.nodes() if network.nodes[n]['group'] == 'ingroup']
outgroup = [n for n in network.nodes() if network.nodes[n]['group'] == 'outgroup']

plot = StaticPlot(network)
plot.color_nodes(ingroup, '#1f78b4')
plot.color_nodes(outgroup, '#ff7f00')
```

### Median Vector Styling

```python
# Different style for median vectors (MJN)
observed = [n for n in network.nodes() if network.nodes[n]['type'] == 'observed']
median = [n for n in network.nodes() if network.nodes[n]['type'] == 'median']

plot = StaticPlot(network)
plot.style_nodes(observed, shape='o', color='#1f78b4')
plot.style_nodes(median, shape='s', color='#ff7f00', size=200)
```

## Troubleshooting

### "Figure too small"

→ Increase `figsize` and `dpi`

### "Overlapping nodes"

→ Try different layout algorithm or increase figure size

### "Labels unreadable"

→ Increase `font_size` or disable labels

### "Colors not distinct"

→ Use fewer categories or colorblind-safe palette

### "File too large"

→ Use vector format (PDF/SVG) or reduce DPI

## Next Steps

- [Analysis Guide](analysis.md): Analyze network properties
- [Tutorials](../tutorials/visualization.md): Detailed examples
- [API Reference](../api/visualization/static.md): Complete API
