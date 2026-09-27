# Frequently Asked Questions

## General

### What is PyPopART?

PyPopART is a pure Python implementation of PopART (Population Analysis with Reticulate Trees), a tool for constructing and visualizing haplotype networks from DNA sequence alignments.

### How is PyPopART different from the original PopART?

- **Pure Python**: No Java dependency, easier to install
- **CLI and API**: Both command-line and programmatic access
- **Modern stack**: Uses current Python libraries (NetworkX, Plotly)
- **Open development**: Active development on GitHub
- **Additional features**: More analysis options and statistics

### Is PyPopART compatible with PopART files?

Yes! PyPopART can read NEXUS files with traits/metadata in PopART format, and can export networks in formats compatible with PopART.

## Installation

### Which Python version do I need?

Python 3.9 or higher. We recommend Python 3.11 or 3.12 for best performance.

### Can I use PyPopART on Windows?

Yes! PyPopART works on Windows, macOS, and Linux.

### Installation fails with NumPy/SciPy errors

Try installing NumPy and SciPy first:

```bash
pip install numpy scipy
pip install pypopart
```

On some systems, you may need system dependencies for NumPy/SciPy. See your OS package manager docs.

## Usage

### Which algorithm should I use?

- **MST**: Simplest, always gives a tree
- **MSN**: Shows alternative connections, slightly more complex
- **TCS**: Statistically justified, may give disconnected networks
- **MJN**: Most comprehensive, infers ancestral haplotypes

See the [Algorithm Guide](guide/algorithms.md) for detailed comparison.

### Which distance metric should I use?

- **Hamming**: Simple differences, fast, good for closely related sequences
- **Jukes-Cantor**: Corrects for multiple substitutions
- **K2P**: Accounts for transition/transversion differences
- **Tamura-Nei**: Most complex, accounts for GC content and rate variation

For most cases, Hamming or K2P is sufficient.

### My network is too complex, how can I simplify it?

For MJN networks:

```bash
pypopart network sequences.fasta -a mjn -e 10 -o network.graphml
```

Increase epsilon (e.g., 5, 10, 20) to reduce complexity.

For TCS networks, increase the parsimony limit:

```bash
pypopart network sequences.fasta -a tcs -p 0.99 -o network.graphml
```

Or use a simpler algorithm (MST or MSN).

### My TCS network is disconnected, is this a bug?

No, this is expected behavior! TCS only connects haplotypes within the parsimony limit. Disconnected components may represent distinct lineages or insufficient data.

### How do I add population/metadata information?

Create a CSV file with columns including 'id' and any metadata:

```csv
id,population,location
Seq1,PopA,Site1
Seq2,PopA,Site1
Seq3,PopB,Site2
```

Then:

```bash
pypopart load sequences.fasta -m metadata.csv
pypopart visualize network.graphml --color-by population -o network.png
```

### Can I visualize networks from other software?

If the network is in GraphML, GML, or JSON format, yes! Load it:

```bash
pypopart visualize external_network.graphml -o plot.png
```

## Visualization

### How do I make publication-quality figures?

Use PDF output for vector graphics:

```bash
pypopart visualize network.graphml -o network.pdf --width 1200 --height 1200
```

For highest quality:

- Use larger dimensions (1200+ pixels)
- Use vector formats (PDF, SVG)
- Show labels only if needed
- Use spring layout for complex networks

### Node sizes don't reflect my sample sizes

Node sizes are based on haplotype frequency. Make sure frequency information is preserved in your network file. When constructing from alignment, frequencies are automatically calculated.

### How do I customize colors?

Currently, colors are assigned automatically. For custom coloring:

```python
from pypopart.visualization import StaticVisualizer

viz = StaticVisualizer(network)
custom_colors = {'Hap1': 'red', 'Hap2': 'blue', ...}
viz.plot(node_colors=custom_colors, output_file='network.png')
```

### Can I export interactive plots?

Yes! Use HTML format:

```bash
pypopart visualize network.graphml -o network.html --interactive
```

Open the HTML file in any web browser for interactive exploration.

## Performance

### PyPopART is slow with my large dataset

For datasets >1000 sequences:

1. Use simpler algorithms (MST, MSN instead of MJN)
2. Use Hamming distance instead of complex models
3. For MJN, increase epsilon to reduce complexity
4. Consider subsampling if haplotype diversity is low

### Can I parallelize the analysis?

Currently, PyPopART runs single-threaded for most operations. Parallel processing may be added in future releases.

## Analysis

### What statistics should I report?

Common statistics include:

- Number of haplotypes
- Haplotype diversity
- Nucleotide diversity
- Network diameter
- Number of reticulations

Use:

```bash
pypopart analyze network.graphml --stats --topology
```

### How do I calculate FST between populations?

```python
from pypopart.stats import PopulationGeneticsAnalysis
from pypopart.io import load_alignment

alignment = load_alignment('sequences.fasta')
popgen = PopulationGeneticsAnalysis(alignment)

# Define populations
populations = {'PopA': ['Seq1', 'Seq2'], 'PopB': ['Seq3', 'Seq4']}

# Calculate FST
fst = popgen.calculate_pairwise_fst(populations)
print(f'FST = {fst:.4f}')
```

### What does reticulation index mean?

Reticulation index measures network complexity:

- 0 = tree (no cycles)
- Higher values = more reticulate (more cycles)

It's calculated as: (edges - nodes + 1) / nodes

## Errors

### "Sequence length doesn't match alignment length"

All sequences in an alignment must have the same length. Check your input file for sequences of different lengths.

### "Sequence ID already exists"

Sequence IDs must be unique. Check for duplicate sequence names in your input file.

### "Distance matrix is not square"

This usually indicates a bug. Please report it on GitHub with your data and command.

### "Cannot import name '**version**'"

The package wasn't installed correctly. Try:

```bash
pip uninstall pypopart
pip install pypopart
```

## File Formats

### What input formats are supported?

- FASTA (.fasta, .fa, .fna)
- NEXUS (.nexus, .nex)
- PHYLIP (.phy, .phylip)
- GenBank (.gb, .gbk)

Compressed files (.gz, .zip) are automatically detected and decompressed.

### What output formats are available?

Networks:

- GraphML (.graphml) - Recommended
- GML (.gml)
- JSON (.json)
- NEXUS (.nexus, .nex)

Visualizations:

- PNG (.png) - Raster
- PDF (.pdf) - Vector, publication-ready
- SVG (.svg) - Vector, web-friendly
- HTML (.html) - Interactive

### Can I convert between formats?

Yes:

```bash
# Load in one format, save in another
pypopart network sequences.fasta -o network.graphml
# Convert to JSON
python -c "
from pypopart.io import load_network, save_network
net = load_network('network.graphml')
save_network(net, 'network.json', format='json')
"
```

## Development

### How can I contribute?

See the [Contributing Guide](contributing.md) for guidelines.

### Where should I report bugs?

Open an issue on [GitHub](https://github.com/adamtaranto/pypopart/issues) with:

- Python and PyPopART versions
- Command or code that caused the error
- Full error message
- Sample data if possible

### Can I request features?

Yes! Open an issue on GitHub describing:

- The use case
- Proposed functionality
- Examples from other tools (if applicable)

## Citation

### How do I cite PyPopART?

```text
Taranto, A. (2024). PyPopART: Pure Python implementation of haplotype network analysis.
GitHub repository: https://github.com/adamtaranto/pypopart
```

BibTeX:

```bibtex
@software{pypopart,
  author = {Taranto, Adam},
  title = {PyPopART: Pure Python implementation of haplotype network analysis},
  year = {2024},
  url = {https://github.com/adamtaranto/pypopart}
}
```

## Still Have Questions?

- Check the [User Guide](guide/cli.md)
- Read the [Tutorials](tutorials/basic_workflow.md)
- Search [GitHub Issues](https://github.com/adamtaranto/pypopart/issues)
- Open a new issue with the "question" label
