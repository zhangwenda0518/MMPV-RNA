# Network Analysis

PyPopART provides comprehensive tools for analyzing haplotype network structure and population genetics.

## Network Statistics

### Basic Properties

```python
from pypopart.stats import NetworkStatistics

stats = NetworkStatistics(network)

# Size metrics
print(f'Number of nodes: {stats.number_of_nodes()}')
print(f'Number of edges: {stats.number_of_edges()}')
print(f'Network density: {stats.density():.3f}')

# Connectivity
print(f'Is connected: {stats.is_connected()}')
print(f'Number of components: {stats.number_of_components()}')

# Path metrics
print(f'Diameter: {stats.diameter()}')
print(f'Average path length: {stats.average_path_length():.3f}')
```

### Centrality Measures

Identify important haplotypes:

```python
# Degree centrality (number of connections)
degree_centrality = stats.degree_centrality()
print('Most connected haplotypes:')
for node, score in sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)[
    :5
]:
    print(f'  {node}: {score:.3f}')

# Betweenness centrality (bridge nodes)
betweenness = stats.betweenness_centrality()

# Closeness centrality (proximity to all nodes)
closeness = stats.closeness_centrality()

# Eigenvector centrality (influence)
eigenvector = stats.eigenvector_centrality()
```

### Clustering

```python
# Clustering coefficient (local interconnectedness)
clustering = stats.clustering_coefficient()
print(f'Average clustering: {clustering:.3f}')

# Per-node clustering
node_clustering = stats.node_clustering()
```

## Topology Analysis

### Hub Detection

Identify central haplotypes:

```python
from pypopart.stats import TopologyAnalysis

topology = TopologyAnalysis(network)

# Find hub haplotypes
hubs = topology.identify_hubs(threshold=3)
print(f'Hub haplotypes: {hubs}')

# Star topology detection
stars = topology.detect_star_patterns()
for center, tips in stars:
    print(f'Star centered on {center} with {len(tips)} tips')
```

### Bridge Analysis

Find edges connecting subpopulations:

```python
# Identify bridge edges
bridges = topology.identify_bridges()
print(f'Bridge connections: {bridges}')

# Articulation points (nodes that disconnect network if removed)
articulation_points = topology.articulation_points()
print(f'Critical nodes: {articulation_points}')
```

### Community Detection

```python
# Detect communities/subgroups
communities = topology.detect_communities(method='louvain')
print(f'Number of communities: {len(communities)}')

for i, community in enumerate(communities):
    print(f'Community {i + 1}: {len(community)} haplotypes')
```

## Population Genetics Statistics

### Diversity Measures

```python
from pypopart.stats import PopulationGenetics

popgen = PopulationGenetics(alignment)

# Haplotype diversity (probability two random samples differ)
h = popgen.haplotype_diversity()
print(f'Haplotype diversity: {h:.4f}')

# Nucleotide diversity (average pairwise differences)
pi = popgen.nucleotide_diversity()
print(f'Nucleotide diversity (π): {pi:.4f}')

# Shannon entropy
shannon = popgen.shannon_entropy()
print(f'Shannon index: {shannon:.4f}')
```

### Neutrality Tests

Test for departure from neutral evolution:

```python
# Tajima's D
tajimas_d = popgen.tajimas_d()
print(f"Tajima's D: {tajimas_d:.4f}")
if tajimas_d < -2:
    print('  → Suggests population expansion or purifying selection')
elif tajimas_d > 2:
    print('  → Suggests balancing selection or population contraction')
else:
    print('  → Consistent with neutral evolution')

# Fu's Fs
fus_fs = popgen.fus_fs()
print(f"Fu's Fs: {fus_fs:.4f}")
if fus_fs < -3:
    print('  → Suggests recent population expansion')
```

### Population Differentiation

Analyze structure between populations:

```python
# FST between populations (requires metadata)
fst = popgen.calculate_fst(population_column='Population')
print(f'FST: {fst:.4f}')

# Pairwise FST
pairwise_fst = popgen.pairwise_fst(population_column='Population')
print('Pairwise FST:')
print(pairwise_fst)

# Gene flow (Nm)
nm = popgen.gene_flow(population_column='Population')
print(f'Gene flow (Nm): {nm:.2f}')
```

### Allele Frequency Spectrum

```python
# Site frequency spectrum
sfs = popgen.site_frequency_spectrum()
print('Site frequency spectrum:', sfs)

# Segregating sites
s = popgen.segregating_sites()
print(f'Number of segregating sites: {s}')
```

## Temporal Analysis

Analyze changes over time:

```python
# Requires temporal metadata
temporal = PopulationGenetics(alignment)

# Diversity over time
diversity_by_year = temporal.diversity_through_time(time_column='Year')

# Expansion/contraction
growth_rate = temporal.estimate_growth_rate(time_column='Year')
print(f'Growth rate: {growth_rate:.4f}')
```

## Spatial Analysis

Analyze geographic structure:

```python
from pypopart.stats import SpatialAnalysis

# Requires geographic metadata
spatial = SpatialAnalysis(alignment, network)

# Isolation by distance
ibd = spatial.isolation_by_distance(lat_column='Latitude', lon_column='Longitude')
print(f'Isolation by distance (R²): {ibd["r_squared"]:.3f}')
print(f'P-value: {ibd["p_value"]:.4f}')

# Geographic structure
amova = spatial.amova(region_column='Region')
print(f'Among-region variance: {amova["among"]:.1f}%')
print(f'Within-region variance: {amova["within"]:.1f}%')
```

## Haplotype Frequency Analysis

```python
# Most common haplotypes
freq = network.haplotype_frequencies()
top_10 = sorted(freq.items(), key=lambda x: x[1], reverse=True)[:10]
print('Top 10 haplotypes by frequency:')
for hap, count in top_10:
    print(f'  {hap}: {count} ({count / sum(freq.values()) * 100:.1f}%)')

# Singleton haplotypes
singletons = [h for h, f in freq.items() if f == 1]
print(f'Number of singleton haplotypes: {len(singletons)}')

# Cumulative frequency
cumulative = popgen.cumulative_frequency_distribution()
```

## Network Comparison

Compare different networks:

```python
from pypopart.stats import NetworkComparison

# Compare two networks
comparison = NetworkComparison(network1, network2)

# Structural similarity
similarity = comparison.structural_similarity()
print(f'Jaccard similarity: {similarity["jaccard"]:.3f}')
print(f'Correlation: {similarity["correlation"]:.3f}')

# Node overlap
node_overlap = comparison.node_overlap()
print(f'Shared nodes: {node_overlap["shared"]}')
print(f'Unique to network1: {node_overlap["unique1"]}')
print(f'Unique to network2: {node_overlap["unique2"]}')

# Topology differences
diff = comparison.topology_differences()
```

## Statistical Testing

### Permutation Tests

```python
from pypopart.stats import PermutationTest

# Test if observed statistic is significant
test = PermutationTest(alignment, n_permutations=1000)

# Test for population differentiation
result = test.test_fst(population_column='Population')
print(f'Observed FST: {result["observed"]:.4f}')
print(f'P-value: {result["p_value"]:.4f}')

# Test for temporal structure
result = test.test_temporal_structure(time_column='Year')
```

### Bootstrap Confidence Intervals

```python
from pypopart.stats import Bootstrap

bootstrap = Bootstrap(alignment, n_replicates=1000)

# Bootstrap diversity estimates
ci = bootstrap.diversity_ci(metric='nucleotide_diversity')
print(f'Nucleotide diversity: {ci["estimate"]:.4f}')
print(f'95% CI: [{ci["lower"]:.4f}, {ci["upper"]:.4f}]')
```

## Export Results

### Summary Report

```python
# Generate comprehensive report
stats = NetworkStatistics(network)
popgen = PopulationGenetics(alignment)

report = {
    'Network': {
        'nodes': stats.number_of_nodes(),
        'edges': stats.number_of_edges(),
        'diameter': stats.diameter(),
        'clustering': stats.clustering_coefficient(),
    },
    'Population': {
        'haplotype_diversity': popgen.haplotype_diversity(),
        'nucleotide_diversity': popgen.nucleotide_diversity(),
        'tajimas_d': popgen.tajimas_d(),
    },
}

# Save as JSON
import json

with open('analysis_report.json', 'w') as f:
    json.dump(report, f, indent=2)

# Save as CSV
import pandas as pd

df = pd.DataFrame(report)
df.to_csv('analysis_report.csv')
```

### Batch Analysis

```python
# Analyze multiple datasets
datasets = ['pop1.fasta', 'pop2.fasta', 'pop3.fasta']
results = []

for dataset in datasets:
    alignment = Alignment.from_fasta(dataset)
    network = algorithm.build_network(alignment)

    stats = NetworkStatistics(network)
    popgen = PopulationGenetics(alignment)

    results.append(
        {
            'dataset': dataset,
            'n_haplotypes': stats.number_of_nodes(),
            'hap_diversity': popgen.haplotype_diversity(),
            'nuc_diversity': popgen.nucleotide_diversity(),
        }
    )

# Save results
df = pd.DataFrame(results)
df.to_csv('batch_analysis.csv', index=False)
```

## Best Practices

1. **Check assumptions**: Ensure statistical tests are appropriate for your data
2. **Multiple testing correction**: Apply when running many tests
3. **Report confidence intervals**: Not just point estimates
4. **Visualize results**: Plots often reveal patterns statistics miss
5. **Biological interpretation**: Statistical significance ≠ biological importance

## Common Analyses

### Publication Checklist

For a typical population genetics paper:

```python
# 1. Basic network properties
stats = NetworkStatistics(network)
print(f'Haplotypes: {stats.number_of_nodes()}')
print(f'Connections: {stats.number_of_edges()}')

# 2. Diversity measures
popgen = PopulationGenetics(alignment)
print(f'Haplotype diversity: {popgen.haplotype_diversity():.4f}')
print(f'Nucleotide diversity: {popgen.nucleotide_diversity():.4f}')

# 3. Neutrality tests
print(f"Tajima's D: {popgen.tajimas_d():.4f}")

# 4. Population structure
print(f'FST: {popgen.calculate_fst(population_column="Population"):.4f}')

# 5. Network topology
topology = TopologyAnalysis(network)
print(f'Hub haplotypes: {topology.identify_hubs()}')
```

## Troubleshooting

### "Insufficient data for test"

→ Ensure adequate sample size (usually n > 20)

### "Test not significant"

→ May indicate true lack of signal or insufficient power

### "Negative diversity"

→ Check for sequencing errors or incorrect alignment

### "FST > 1 or < 0"

→ Calculation error; check population assignments

## Next Steps

- [Visualization Guide](visualization.md): Plot analysis results
- [Tutorials](../tutorials/popgen.md): Population genetics examples
- [API Reference](../api/stats/statistics.md): Complete statistics API
