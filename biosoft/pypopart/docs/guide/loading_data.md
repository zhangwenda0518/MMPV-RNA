# Loading Data

PyPopART supports multiple sequence file formats and provides flexible data loading options.

## Supported Formats

### FASTA Format

The most common sequence format:

```fasta
>Sequence1_PopA
ATCGATCGATCGATCGATCG
>Sequence2_PopA
ATCGATCGATCGATCGATCG
>Sequence3_PopB
ATCGATCGATCGATTGATCG
```

**Load in Python:**

```python
from pypopart import Alignment

alignment = Alignment.from_fasta('sequences.fasta')
```

**CLI:**

```bash
pypopart network sequences.fasta -a MST -o output
```

### NEXUS Format

Supports metadata and traits:

```nexus
#NEXUS
BEGIN TAXA;
    DIMENSIONS NTAX=3;
    TAXLABELS Seq1 Seq2 Seq3;
END;

BEGIN CHARACTERS;
    DIMENSIONS NCHAR=20;
    FORMAT DATATYPE=DNA MISSING=? GAP=-;
    MATRIX
        Seq1 ATCGATCGATCGATCGATCG
        Seq2 ATCGATCGATCGATCGATCG
        Seq3 ATCGATCGATCGATTGATCG
    ;
END;

BEGIN TRAITS;
    DIMENSIONS NTRAITS=2;
    FORMAT LABELS=YES SEPARATOR=,;
    TRAITLABELS Population Location;
    MATRIX
        Seq1 PopA Site1
        Seq2 PopA Site1
        Seq3 PopB Site2
    ;
END;
```

**Load in Python:**

```python
alignment = Alignment.from_nexus('sequences.nex')
metadata = alignment.get_metadata()
```

### PHYLIP Format

Sequential or interleaved:

```phylip
3 20
Seq1      ATCGATCGATCGATCGATCG
Seq2      ATCGATCGATCGATCGATCG
Seq3      ATCGATCGATCGATTGATCG
```

**Load in Python:**

```python
alignment = Alignment.from_phylip('sequences.phy')
```

### GenBank Format

Full GenBank entries:

```python
alignment = Alignment.from_genbank('sequences.gb')
```

## Working with Metadata

### Including Metadata

Metadata can encode population, location, time, or custom traits:

**Method 1: NEXUS Traits Block**

```python
# Automatically loaded from NEXUS file
alignment = Alignment.from_nexus('sequences_with_traits.nex')
print(alignment.metadata)
```

**Method 2: Add Programmatically**

```python
import pandas as pd

# Load sequences
alignment = Alignment.from_fasta('sequences.fasta')

# Add metadata
metadata = pd.DataFrame(
    {
        'sequence_id': ['Seq1', 'Seq2', 'Seq3'],
        'Population': ['PopA', 'PopA', 'PopB'],
        'Location': ['Site1', 'Site1', 'Site2'],
        'Year': [2020, 2020, 2021],
    }
)

alignment.set_metadata(metadata)
```

**Method 3: Parse from Sequence Names**

```python
# If names are like: "Sample1_PopA_Site1"
alignment = Alignment.from_fasta('sequences.fasta')
alignment.parse_names(pattern=r'(?P<sample>\w+)_(?P<population>\w+)_(?P<location>\w+)')
```

### Using Metadata for Analysis

```python
# Color networks by metadata
from pypopart.visualization import StaticPlot

plot = StaticPlot(network)
plot.color_by_attribute('Population')
plot.save('colored_network.png')

# Calculate population statistics
from pypopart.stats import PopulationGenetics

popgen = PopulationGenetics(alignment)
fst = popgen.calculate_fst(population_column='Population')
```

## Data Validation

### Check Alignment Quality

```python
# Verify alignment
print(f'Number of sequences: {len(alignment)}')
print(f'Alignment length: {alignment.length}')
print(f'Valid alignment: {alignment.is_valid()}')

# Check for gaps
if alignment.has_gaps():
    print('Warning: Alignment contains gaps')

# Check for ambiguous bases
if alignment.has_ambiguous():
    print('Warning: Alignment contains ambiguous bases')
```

### Handle Missing Data

```python
# Remove sequences with excessive gaps
alignment = alignment.filter_by_gaps(max_gap_fraction=0.1)

# Remove columns with excessive missing data
alignment = alignment.filter_columns(max_missing=0.2)

# Remove invariant sites
alignment = alignment.remove_invariant_sites()
```

## File Format Detection

PyPopART can auto-detect formats:

```python
# Auto-detect format
alignment = Alignment.from_file('sequences.unknown')
```

Or explicitly specify:

```python
alignment = Alignment.from_file('sequences.txt', format='fasta')
```

## Large Files

### Streaming Data

For very large files:

```python
# Process in chunks
for chunk in Alignment.read_chunks('large_file.fasta', chunk_size=1000):
    # Process each chunk
    network = algorithm.build_network(chunk)
```

### Memory Optimization

```python
# Disable unnecessary features
alignment = Alignment.from_fasta(
    'sequences.fasta', load_metadata=False, compute_stats=False
)
```

## BioPython Integration

Convert between PyPopART and BioPython:

```python
from Bio import AlignIO
from pypopart import Alignment

# BioPython to PyPopART
bio_aln = AlignIO.read('sequences.fasta', 'fasta')
pp_aln = Alignment.from_biopython(bio_aln)

# PyPopART to BioPython
bio_aln = pp_aln.to_biopython()
AlignIO.write(bio_aln, 'output.fasta', 'fasta')
```

## Data Requirements

### Sequence Data

- **Aligned sequences**: All sequences must be the same length
- **DNA/Protein**: PyPopART handles both (specify datatype if needed)
- **No special characters**:除 standard IUPAC codes

### Metadata (Optional)

- **Categorical traits**: Population, location, phenotype
- **Numerical traits**: Year, coordinates, measurements
- **Missing values**: Use 'NA', '?', or leave blank

## Example Datasets

PyPopART includes example datasets:

```python
from pypopart.data import load_example

# Load example data
alignment = load_example('woodmouse')
alignment = load_example('influenza')
alignment = load_example('mtdna')

# Get example file path
path = load_example('woodmouse', return_path=True)
```

## Tips for Data Preparation

1. **Use aligned sequences**: Run alignment tools first (MUSCLE, MAFFT, Clustal)
2. **Include metadata**: Enriches analysis and visualization
3. **Check quality**: Remove poor quality sequences
4. **Remove duplicates**: PyPopART will collapse identical sequences
5. **Consistent naming**: Use systematic sequence names
6. **Check encoding**: Ensure UTF-8 encoding for special characters

## Troubleshooting

### Common Errors

**"Sequences not aligned"**

- All sequences must be same length
- Align sequences before loading

**"Invalid NEXUS format"**

- Check NEXUS syntax
- Ensure all blocks are properly closed

**"Cannot parse metadata"**

- Check trait names match sequence IDs
- Verify separator (comma, tab, space)

**"Memory error"**

- File too large
- Use streaming or chunking
- Filter data first

## Next Steps

- [Distance Metrics Guide](distances.md): Choose appropriate distance
- [Algorithm Guide](algorithms.md): Select network algorithm
- [Visualization Guide](visualization.md): Plot your networks
