# PyPopART GUI Usage Guide

## Overview

PyPopART provides a web-based graphical user interface (GUI) built with Dash for interactive haplotype network analysis. The GUI offers an intuitive workflow for:

1. Uploading sequence alignments
2. Computing haplotype networks
3. Visualizing and analyzing results
4. Exporting networks and figures

## Installation

Ensure you have PyPopART installed with GUI dependencies:

```bash
pip install pypopart
# Or install in development mode:
pip install -e .
```

The GUI requires the following additional dependencies (automatically installed):

- dash
- dash-bootstrap-components

## Launching the GUI

### From Command Line

The easiest way to start the GUI:

```bash
# Start on default port 8050
pypopart-gui

# Start on custom port
pypopart-gui --port 8080

# Enable debug mode
pypopart-gui --debug
```

### From Python

```python
from pypopart.gui.app import main

# Launch with default settings (port 8050, debug=False)
main()

# Or customize:
main(debug=True, port=8080)
```

The GUI will start a local web server. Open your browser to `http://localhost:8050` (or your specified port).

## Using the GUI

### 1. Upload Data

#### Sequence File (Required)

- Click "📁 Select Sequence File"
- Supported formats:
  - FASTA (.fasta, .fa)
  - NEXUS (.nex, .nexus)
  - PHYLIP (.phy, .phylip)
- After upload, you'll see a success message with sequence count and alignment length

#### Metadata File (Optional)

- Click "📊 Select Metadata File"
- Upload a CSV file with population, location, or trait data
- Required columns: `id` (matching sequence IDs)
- Optional columns: `population`, `latitude`, `longitude`, `color`, `notes`
- Population data enables:
  - Colored nodes by population
  - Pie charts for nodes with mixed populations
  - Population-based statistics

### 2. Configure Algorithm

Choose from six network construction algorithms:

#### MST (Minimum Spanning Tree)

- **Parameters**: Distance metric (Hamming, Jukes-Cantor, Kimura 2-parameter)
- **Best for**: Simplest tree-like relationships without reticulation

#### MSN (Minimum Spanning Network)

- **Parameters**: Distance metric
- **Best for**: Showing alternative connections at equal distance

#### TCS (Statistical Parsimony)

- **Parameters**: Connection limit (1-20 mutations)
- **Best for**: Recent divergence, intraspecific analysis
- **Default**: 10 mutations

#### MJN (Median-Joining Network)

- **Parameters**: Epsilon value (0 = automatic)
- **Best for**: Complex reticulate relationships, infers median vectors

#### PN (Parsimony Network)

- **Parameters**: Number of trees to sample (10-500)
- **Best for**: Consensus approach capturing phylogenetic uncertainty
- **Default**: 100 trees

#### TSW (Tight Span Walker)

- **Parameters**: Distance metric
- **Best for**: Metric-preserving networks, small to medium datasets
- **Note**: Computationally intensive

After selecting an algorithm and parameters, click **⚡ Compute Network**.

### 3. Layout Options

Choose how to arrange nodes in the network visualization:

#### Available Layouts

- **Hierarchical (Fast)**: Tree-like hierarchical arrangement, quick computation
- **Spring (Force-directed)**: Physics-based layout, good for general use
- **Spring - Proportional Edge Length**: Spring layout where edge length reflects mutation distance
- **Spectral (Fast, Large networks)**: Eigenvalue-based layout, efficient for large networks
- **Circular**: Nodes arranged in a circle
- **Radial**: Concentric rings from center
- **Kamada-Kawai (High quality, slow)**: Energy-minimization layout, best visual quality
- **Kamada-Kawai - Proportional Edge Length**: KK layout where edge length reflects mutation distance

#### Customization Options

- **Node Spacing**: Adjust spacing between nodes (0.5x - 3.0x)
- **Node Size**: Adjust node size (10 - 100)
- **Edge Width**: Adjust edge thickness (1 - 10)

Click **🎨 Apply Layout** to recompute node positions with new settings.

### 4. Viewing Results

The GUI provides five tabs:

#### Network Tab

- **Interactive visualization** with Dash Cytoscape
- **Node styling**:
  - Size: Proportional to haplotype frequency
  - Color: Single population = solid color, mixed = pie chart
  - Shape: Median vectors shown as gray circles
  - Selection: Click to select (red border), search to highlight
- **Edge styling**:
  - Labels: Show mutation count
  - Thickness: Adjustable via slider
- **Interactions**:
  - Zoom: Scroll wheel
  - Pan: Click and drag background
  - Move nodes: Drag individual nodes to reposition
  - Hover: View haplotype details in tooltip
  - Search: Use dropdown to find and highlight specific haplotypes
- **Legend**: Shows population colors and mixed population indicator

#### Statistics Tab

Displays comprehensive network metrics:

- **Basic Metrics**: Nodes, edges, diameter, clustering coefficient, reticulation index
- **Diversity Metrics**: Haplotype diversity, Shannon index
- **Central Haplotypes**: Degree, betweenness, and closeness centrality measures

#### Haplotype Summary Tab

- **Table view**: H number, type (observed/inferred), frequency, sample IDs, populations
- **Export options**:
  - Download haplotype summary as CSV
  - Download H number label template
  - Upload custom H number labels
- Shows which sequences belong to each haplotype
- Identifies inferred median nodes

#### Metadata Tab

- View uploaded metadata aligned with sequence IDs
- Shows which IDs are in alignment vs. metadata
- Displays population assignments and colors
- Warns about mismatches or duplicates

#### Alignment Tab

- View the uploaded sequence alignment
- Polymorphic sites highlighted in color (A=green, C=blue, G=orange, T=red)
- Shows first 50 sequences with ID and sequence data aligned

### 5. Export

Export your results in various formats:

- **GraphML**: Network format, opens in Cytoscape, Gephi
- **GML**: Graph Markup Language
- **JSON**: JavaScript Object Notation, for web applications
- **PNG**: Raster image, good for presentations
- **SVG**: Vector image, scalable for publications

Click **Download** to save the file.

## Tips and Best Practices

### Data Preparation

- Ensure sequences are aligned before uploading
- All sequences should have the same length
- Include metadata in NEXUS format for population coloring

### Algorithm Selection

- Start with **MSN** for most datasets (good balance)
- Use **TCS** for closely related sequences (e.g., same species)
- Use **MJN** for complex evolutionary scenarios
- **MST** is fastest but simplest

### Layout Optimization

- Try different layout algorithms to find the clearest visualization
- **Spring layout** works well for most networks
- **Radial layout** highlights star-like patterns
- Enable **Snap to Grid** for publication-ready figures

### Performance

- Large datasets (>100 sequences) may take time to compute
- Progress indicators show when computation is running
- Network visualization is interactive even with many nodes

### Troubleshooting

**File Upload Fails**

- Check file format and extension
- Ensure sequences are properly formatted
- Try a smaller test dataset first

**Network Computation Fails**

- Check that all sequences have the same length
- Reduce connection limit for TCS if very large
- Try a different algorithm

**Visualization Issues**

- Refresh the browser page
- Try a different layout algorithm
- Check browser console for errors (F12)

**Export Doesn't Work**

- Ensure network has been computed
- Check browser's download settings
- Try a different export format

## Advanced Features

### Customization

The GUI can be customized by modifying `src/pypopart/gui/app.py`:

- Adjust default parameters
- Change color schemes
- Add custom layout algorithms
- Modify statistics displayed

### Integration

The GUI components can be used programmatically:

```python
from pypopart.gui.app import PyPopARTApp

app = PyPopARTApp(debug=True, port=8050)
# Access app.app for the Dash application object
# Customize before running
app.run()
```

## Support

For issues, questions, or feature requests:

- GitHub Issues: <https://github.com/adamtaranto/pypopart/issues>
- Documentation: <https://github.com/adamtaranto/pypopart>

## Citation

If you use PyPopART in your research, please cite:

[Citation information to be added]
