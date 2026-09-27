# Network Algorithms

PyPopART implements six algorithms for constructing haplotype networks. Each has different characteristics and is suited for different types of data.

## Algorithm Overview

| Algorithm | Complexity | Reticulation | Median Vectors | Best For               |
| --------- | ---------- | ------------ | -------------- | ---------------------- |
| MST       | Low        | No           | No             | Initial exploration    |
| MSN       | Low-Medium | Yes          | No             | Alternative paths      |
| TCS       | Medium     | Yes          | No             | Within-species         |
| MJN       | High       | Yes          | Yes            | Comprehensive analysis |
| PN        | Medium     | Yes          | No             | General purpose        |
| TSW       | Medium     | Yes          | No             | Distance preservation  |

## Minimum Spanning Tree (MST)

Creates the simplest tree connecting all haplotypes with minimum total distance.

### Characteristics

- **Always produces a tree** (no cycles)
- **Deterministic** (same result each time)
- **Fast** computation
- **No reticulation** (alternative paths)

### Algorithm

1. Start with all haplotypes as separate components
2. Repeatedly add shortest edge that connects components
3. Stop when all haplotypes are connected

### When to Use

- Initial data exploration
- Small datasets (< 50 haplotypes)
- When tree structure is appropriate
- Quick preliminary analysis

### Python

```python
from pypopart.algorithms import MSTAlgorithm

algorithm = MSTAlgorithm()
network = algorithm.build_network(alignment)
```

### CLI

```bash
pypopart network sequences.fasta -a MST -o mst_network
```

### Limitations

- Ignores alternative equally parsimonious connections
- May oversimplify complex relationships
- Single path between haplotypes

## Minimum Spanning Network (MSN)

Extends MST by adding alternative connections of equal parsimony.

### Characteristics

- **Shows alternative paths**
- **Reticulation present**
- **Still relatively simple**
- **Fast computation**

### Algorithm

1. Build MST as foundation
2. Add edges creating cycles if distance equals shortest path
3. Remove redundant edges

### When to Use

- Want to see alternative evolutionary paths
- Data has ambiguous relationships
- Balance between simplicity and completeness
- Following MST analysis

### Python

```python
from pypopart.algorithms import MSNAlgorithm

algorithm = MSNAlgorithm()
network = algorithm.build_network(alignment)
```

### CLI

```bash
pypopart network sequences.fasta -a MSN -o msn_network
```

### Advantages over MST

- Shows multiple equally parsimonious paths
- More informative about uncertainty
- Still computationally efficient

## Statistical Parsimony (TCS)

Connects haplotypes within statistical parsimony limit (95% confidence by default).

### Characteristics

- **Statistically justified connections**
- **May produce disconnected networks**
- **No median vectors**
- **Confidence-based**

### Algorithm

1. Calculate maximum number of differences for connection (epsilon)
2. Connect haplotypes differing by ≤ epsilon mutations
3. Nested clade structure

### Parameters

- `epsilon`: Connection limit (default: 0.95 confidence)
- Calculated from alignment or specified

### When to Use

- Within-species variation
- Recent divergence
- Need statistical justification
- Population-level studies

### Python

```python
from pypopart.algorithms import TCSAlgorithm

# Default epsilon (95% confidence)
algorithm = TCSAlgorithm()
network = algorithm.build_network(alignment)

# Custom epsilon
algorithm = TCSAlgorithm(epsilon=0.99)
network = algorithm.build_network(alignment)
```

### CLI

```bash
# Default
pypopart network sequences.fasta -a TCS -o tcs_network

# Custom epsilon
pypopart network sequences.fasta -a TCS --epsilon 0.99 -o tcs_network
```

### Considerations

- May produce multiple disconnected components
- Sensitive to epsilon parameter
- Assumes neutrality and panmixia

## Median-Joining Network (MJN)

Most comprehensive algorithm inferring ancestral haplotypes.

### Characteristics

- **Infers median vectors** (ancestral/unobserved haplotypes)
- **Shows reticulation**
- **Most complete**
- **Computationally intensive**

### Algorithm

1. Build MSN as starting point
2. Identify and add median vectors (consensus sequences)
3. Remove unnecessary edges
4. Iterate until stable

### When to Use

- Complex evolutionary scenarios
- Large datasets with structure
- Need ancestral sequence inference
- Publication-quality analysis

### Python

```python
from pypopart.algorithms import MJNAlgorithm

algorithm = MJNAlgorithm()
network = algorithm.build_network(alignment)
```

### CLI

```bash
pypopart network sequences.fasta -a MJN -o mjn_network
```

### Advantages

- Most informative
- Shows inferred ancestors
- Handles complex relationships
- Widely used in publications

### Limitations

- Computationally expensive
- Can be complex to interpret
- May overfit small datasets

## Parsimony Network (PN)

Consensus approach using multiple MSTs with different starting points.

### Characteristics

- **Consensus of multiple trees**
- **Shows reticulation**
- **Balanced complexity**
- **Robust to starting conditions**

### Algorithm

1. Build multiple MSTs with different starting haplotypes
2. Combine into consensus network
3. Keep connections appearing in multiple trees

### When to Use

- General purpose analysis
- Want robustness
- Balance between MSN and MJN
- Uncertain about best algorithm

### Python

```python
from pypopart.algorithms import ParsimonyNetAlgorithm

algorithm = ParsimonyNetAlgorithm()
network = algorithm.build_network(alignment)
```

### CLI

```bash
pypopart network sequences.fasta -a PN -o pn_network
```

## Tight Span Walker (TSW)

Metric-preserving network construction.

### Characteristics

- **Preserves distance relationships**
- **Metric space approach**
- **Mathematical rigor**
- **Moderate complexity**

### Algorithm

1. Compute tight span of distance matrix
2. Walk through tight span structure
3. Build network preserving distances

### When to Use

- Distance accuracy critical
- Mathematical rigor needed
- Metric properties important
- Complex distance patterns

### Python

```python
from pypopart.algorithms import TSWAlgorithm

algorithm = TSWAlgorithm()
network = algorithm.build_network(alignment)
```

### CLI

```bash
pypopart network sequences.fasta -a TSW -o tsw_network
```

## Choosing an Algorithm

### By Data Type

**Intraspecific (within species)**

- TCS (recommended)
- MSN
- MJN

**Interspecific (between species)**

- MST
- MSN
- PN

**Population genetics**

- TCS
- MJN

**Phylogeography**

- MJN
- PN

### By Dataset Size

**Small (< 20 haplotypes)**

- Any algorithm
- Start with MST/MSN

**Medium (20-100 haplotypes)**

- TCS
- MSN
- MJN

**Large (> 100 haplotypes)**

- MST (fast exploration)
- TCS
- PN

### By Question

**Exploratory analysis**
→ MST, MSN

**Population structure**
→ TCS, MJN

**Ancestral inference**
→ MJN

**Publication figure**
→ MJN, PN

**Quick overview**
→ MST

**Statistical rigor**
→ TCS

## Comparing Algorithms

Run multiple algorithms on the same data:

```python
from pypopart import Alignment
from pypopart.algorithms import MSTAlgorithm, MSNAlgorithm, TCSAlgorithm, MJNAlgorithm

alignment = Alignment.from_fasta('sequences.fasta')

algorithms = {
    'MST': MSTAlgorithm(),
    'MSN': MSNAlgorithm(),
    'TCS': TCSAlgorithm(),
    'MJN': MJNAlgorithm(),
}

networks = {}
for name, algorithm in algorithms.items():
    networks[name] = algorithm.build_network(alignment)
    print(
        f'{name}: {networks[name].number_of_nodes()} nodes, {networks[name].number_of_edges()} edges'
    )
```

CLI batch comparison:

```bash
for alg in MST MSN TCS MJN; do
    pypopart network sequences.fasta -a $alg -o comparison/${alg}_network
done
```

## Performance Comparison

| Algorithm | Time (50 seqs) | Time (200 seqs) | Memory |
| --------- | -------------- | --------------- | ------ |
| MST       | < 1s           | < 5s            | Low    |
| MSN       | < 2s           | < 10s           | Low    |
| TCS       | < 5s           | < 30s           | Medium |
| MJN       | 10-30s         | 2-5 min         | High   |
| PN        | 5-10s          | 30-60s          | Medium |
| TSW       | 5-15s          | 30-90s          | Medium |

_Times are approximate and depend on divergence level_

## Best Practices

1. **Start simple**: Begin with MST for exploration
2. **Compare algorithms**: Run multiple to understand data
3. **Check assumptions**: Ensure algorithm matches data type
4. **Consider computational cost**: MJN great but slow
5. **Validate results**: Check if network makes biological sense
6. **Report parameters**: Document algorithm and settings used

## Common Issues

### "Network is disconnected"

- Normal for TCS with high divergence
- Increase epsilon or use different algorithm

### "Too many median vectors"

- MJN can generate many inferred nodes
- Normal for large, complex datasets
- Consider MSN or PN for simplicity

### "Computation too slow"

- MJN slow on large datasets
- Try MST/MSN first
- Consider sampling or filtering data

## References

- Kruskal (1956) - MST
- Excoffier & Smouse (1994) - MSN
- Templeton et al. (1992) - TCS
- Bandelt et al. (1999) - Median-Joining
- Dress & Huson (2004) - Tight Span

## Next Steps

- [Visualization Guide](visualization.md): Plot your networks
- [Analysis Guide](analysis.md): Analyze network properties
- [Tutorials](../tutorials/algorithm_comparison.md): Detailed comparison
