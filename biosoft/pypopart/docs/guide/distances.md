# Distance Metrics

Genetic distance measures quantify differences between DNA sequences. Choosing the right metric is crucial for accurate network construction.

## Available Metrics

### Hamming Distance

Simple count of differing positions between sequences.

**Formula:** Number of sites where sequences differ

**Example:**

```
Seq1: ATCGATCG
Seq2: ATCGTTCG
         ^
Distance: 1
```

**When to use:**

- Very similar sequences (< 5% divergence)
- Initial exploration
- No model assumptions needed

**Python:**

```python
from pypopart.core.distance import DistanceCalculator

calc = DistanceCalculator(metric='hamming')
distances = calc.calculate(alignment)
```

**CLI:**

```bash
pypopart distance sequences.fasta -m hamming -o distances.csv
```

### Jukes-Cantor (JC69)

Corrects for multiple mutations at the same site.

**Formula:**

```
d = -3/4 * ln(1 - 4p/3)
```

where p is the proportion of different sites

**Assumptions:**

- Equal base frequencies (25% each)
- Equal substitution rates
- No rate variation

**When to use:**

- Moderate divergence (5-20%)
- Need model-based correction
- Simple correction sufficient

**Python:**

```python
calc = DistanceCalculator(metric='jukes-cantor')
distances = calc.calculate(alignment)
```

**CLI:**

```bash
pypopart distance sequences.fasta -m jukes-cantor -o distances.csv
```

### Kimura 2-Parameter (K2P)

Distinguishes between transitions and transversions.

**Formula:**

```
d = -1/2 * ln((1-2P-Q) * sqrt(1-2Q))
```

where:

- P = transition frequency
- Q = transversion frequency

**Assumptions:**

- Equal base frequencies
- Different rates for transitions vs transversions
- No rate variation

**When to use:**

- Moderate to high divergence (up to 50%)
- Transitions more common than transversions (typical for DNA)
- More realistic than Jukes-Cantor

**Python:**

```python
calc = DistanceCalculator(metric='k2p')
distances = calc.calculate(alignment)
```

**CLI:**

```bash
pypopart distance sequences.fasta -m k2p -o distances.csv
```

### Tamura-Nei

Most sophisticated model with unequal base frequencies.

**Formula:** Accounts for:

- Unequal base frequencies
- Different transition rates
- Transition/transversion ratio
- GC content bias

**Assumptions:**

- Variable base frequencies
- Different rates for transitions vs transversions
- Rate variation across sites

**When to use:**

- High divergence (> 20%)
- Biased base composition
- Most accurate correction needed
- Sufficient data for parameter estimation

**Python:**

```python
calc = DistanceCalculator(metric='tamura-nei')
distances = calc.calculate(alignment)
```

**CLI:**

```bash
pypopart distance sequences.fasta -m tamura-nei -o distances.csv
```

## Choosing a Metric

### Decision Guide

```
Start here
    ↓
Is divergence < 5%?
    YES → Use Hamming
    NO → Continue
         ↓
Is divergence < 20%?
    YES → Use Jukes-Cantor or K2P
    NO → Continue
         ↓
Is there GC bias?
    YES → Use Tamura-Nei
    NO → Use K2P
```

### By Divergence Level

| Divergence | Recommended Metric | Reason                  |
| ---------- | ------------------ | ----------------------- |
| 0-5%       | Hamming            | Simple, no saturation   |
| 5-20%      | Jukes-Cantor / K2P | Model-based correction  |
| 20-50%     | K2P / Tamura-Nei   | Accounts for saturation |
| > 50%      | Tamura-Nei         | Most sophisticated      |

### By Data Type

| Data Type         | Recommended       | Notes                   |
| ----------------- | ----------------- | ----------------------- |
| Mitochondrial DNA | K2P               | High transition bias    |
| Nuclear DNA       | K2P or Tamura-Nei | Variable composition    |
| Coding regions    | K2P               | Consider codon position |
| Non-coding        | Tamura-Nei        | Often GC-biased         |
| Microsatellites   | Hamming           | Step-wise mutation      |

## Working with Distance Matrices

### Calculate and Save

```python
from pypopart.core.distance import DistanceCalculator

calc = DistanceCalculator(metric='k2p')
distances = calc.calculate(alignment)

# Save matrix
distances.to_csv('distances.csv')
distances.to_numpy('distances.npy')
```

### Load Pre-computed Distances

```python
import numpy as np
from pypopart import DistanceMatrix

# Load from file
distances = DistanceMatrix.from_csv('distances.csv')

# Use with algorithm
from pypopart.algorithms import MSTAlgorithm

algorithm = MSTAlgorithm(distance_matrix=distances)
network = algorithm.build_network(alignment)
```

### Inspect Distance Distribution

```python
import matplotlib.pyplot as plt

# Get distance values
values = distances.values()

# Plot histogram
plt.hist(values, bins=30)
plt.xlabel('Genetic Distance')
plt.ylabel('Frequency')
plt.title('Distance Distribution')
plt.show()

# Summary statistics
print(f'Min distance: {min(values):.4f}')
print(f'Max distance: {max(values):.4f}')
print(f'Mean distance: {np.mean(values):.4f}')
print(f'Median distance: {np.median(values):.4f}')
```

## Advanced Options

### Custom Distance Functions

Define your own distance metric:

```python
from pypopart.core.distance import BaseDistance


class CustomDistance(BaseDistance):
    def calculate_pairwise(self, seq1, seq2):
        # Your calculation here
        diff = sum(a != b for a, b in zip(seq1, seq2))
        # Apply custom correction
        return diff * custom_factor


# Use custom metric
calc = DistanceCalculator(metric=CustomDistance())
distances = calc.calculate(alignment)
```

### Gap Handling

Control how gaps are treated:

```python
calc = DistanceCalculator(
    metric='k2p',
    gap_mode='pairwise',  # 'pairwise', 'complete', 'ignore'
)
```

Options:

- `pairwise`: Remove gap positions for each pair
- `complete`: Remove all sites with any gaps
- `ignore`: Treat gaps as fifth character

### Rate Variation

Account for among-site rate variation:

```python
calc = DistanceCalculator(
    metric='k2p',
    gamma=True,  # Use gamma distribution
    alpha=0.5,  # Shape parameter
)
```

## Validation and Diagnostics

### Check for Saturation

```python
# Plot uncorrected vs corrected distances
hamming_calc = DistanceCalculator(metric='hamming')
k2p_calc = DistanceCalculator(metric='k2p')

hamming_dist = hamming_calc.calculate(alignment)
k2p_dist = k2p_calc.calculate(alignment)

# If K2P >> Hamming, saturation is present
ratio = k2p_dist.mean() / hamming_dist.mean()
if ratio > 1.5:
    print('Warning: Significant saturation detected')
```

### Compare Metrics

```python
metrics = ['hamming', 'jukes-cantor', 'k2p', 'tamura-nei']
results = {}

for metric in metrics:
    calc = DistanceCalculator(metric=metric)
    dist = calc.calculate(alignment)
    results[metric] = dist.mean()

print('Mean distances by metric:')
for metric, mean_dist in results.items():
    print(f'  {metric}: {mean_dist:.4f}')
```

## Performance Considerations

### Speed vs Accuracy

| Metric       | Speed    | Accuracy | Use Case                        |
| ------------ | -------- | -------- | ------------------------------- |
| Hamming      | Fastest  | Basic    | Exploration, low divergence     |
| Jukes-Cantor | Fast     | Good     | Moderate divergence             |
| K2P          | Moderate | Better   | General purpose                 |
| Tamura-Nei   | Slower   | Best     | High divergence, final analysis |

### Optimization Tips

```python
# For large datasets
calc = DistanceCalculator(
    metric='k2p',
    parallel=True,  # Use multiprocessing
    n_jobs=4,  # Number of cores
)

# Cache results
calc.calculate(alignment, cache=True)
```

## Common Issues

### "Distance > 1.0"

Some corrected distances can exceed 1.0 for highly divergent sequences:

- Expected for K2P and Tamura-Nei
- Indicates high divergence
- Consider if sequences are too divergent for network analysis

### "Negative distance"

Rare numerical issue with some corrections:

- Check for data quality issues
- Try simpler metric
- May indicate sequences are too divergent

### "NaN values"

Usually caused by:

- All gaps in pairwise comparison
- Identical sequences (distance = 0)
- Numerical underflow in calculations

## Best Practices

1. **Start simple**: Use Hamming for exploration
2. **Match divergence**: Choose metric appropriate for your data
3. **Check assumptions**: Verify metric assumptions match your data
4. **Validate results**: Compare with multiple metrics
5. **Document choice**: Report which metric and why in publications

## References

- Jukes TH, Cantor CR (1969). Evolution of protein molecules
- Kimura M (1980). A simple method for estimating evolutionary rates
- Tamura K, Nei M (1993). Estimation of the number of nucleotide substitutions

## Next Steps

- [Algorithm Guide](algorithms.md): Choose network algorithm
- [Visualization Guide](visualization.md): Plot results
- [API Reference](../api/core/distance.md): Distance API documentation
