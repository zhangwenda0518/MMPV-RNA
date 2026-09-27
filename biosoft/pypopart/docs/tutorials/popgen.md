# Population Genetics Tutorial

Comprehensive population genetics analysis with PyPopART.

## Overview

This tutorial demonstrates:

- Diversity measures
- Neutrality tests
- Population structure analysis
- FST calculation
- Demographic inference

## Setup

```python
from pypopart import Alignment
from pypopart.algorithms import TCSAlgorithm
from pypopart.stats import PopulationGenetics, NetworkStatistics
from pypopart.visualization import StaticPlot
import pandas as pd
import matplotlib.pyplot as plt
```

## Load Data with Metadata

```python
# Load sequences
alignment = Alignment.from_nexus('sequences_with_traits.nex')

# Or add metadata programmatically
alignment = Alignment.from_fasta('sequences.fasta')

metadata = pd.DataFrame(
    {
        'sequence_id': ['Seq1', 'Seq2', 'Seq3', 'Seq4', 'Seq5', 'Seq6'],
        'Population': ['PopA', 'PopA', 'PopA', 'PopB', 'PopB', 'PopC'],
        'Location': ['Site1', 'Site1', 'Site1', 'Site2', 'Site2', 'Site3'],
        'Year': [2020, 2020, 2020, 2021, 2021, 2022],
    }
)

alignment.set_metadata(metadata)
```

## Diversity Analysis

```python
popgen = PopulationGenetics(alignment)

# Haplotype diversity (Nei's h)
h = popgen.haplotype_diversity()
print(f'Haplotype diversity: {h:.4f}')

# Nucleotide diversity (π)
pi = popgen.nucleotide_diversity()
print(f'Nucleotide diversity: {pi:.6f}')

# Theta (Watterson's estimator)
theta = popgen.theta_watterson()
print(f'Theta (θw): {theta:.6f}')

# Shannon entropy
shannon = popgen.shannon_entropy()
print(f'Shannon index: {shannon:.4f}')
```

## Neutrality Tests

```python
# Tajima's D
tajimas_d = popgen.tajimas_d()
print(f"\nTajima's D: {tajimas_d:.4f}")

if tajimas_d < -2:
    print('  → Population expansion or purifying selection')
elif tajimas_d > 2:
    print('  → Balancing selection or population contraction')
else:
    print('  → Consistent with neutral evolution')

# Fu's Fs
fus_fs = popgen.fus_fs()
print(f"\nFu's Fs: {fus_fs:.4f}")

if fus_fs < -3:
    print('  → Evidence for population expansion')
else:
    print('  → No strong evidence for expansion')

# Fu and Li's D* and F*
fu_li_d = popgen.fu_li_d_star()
fu_li_f = popgen.fu_li_f_star()
print(f"\nFu and Li's D*: {fu_li_d:.4f}")
print(f"Fu and Li's F*: {fu_li_f:.4f}")
```

## Population Structure

```python
# FST between all populations
fst_total = popgen.calculate_fst(population_column='Population')
print(f'\nOverall FST: {fst_total:.4f}')

if fst_total < 0.05:
    print('  → Little differentiation')
elif fst_total < 0.15:
    print('  → Moderate differentiation')
elif fst_total < 0.25:
    print('  → Great differentiation')
else:
    print('  → Very great differentiation')

# Pairwise FST
pairwise_fst = popgen.pairwise_fst(population_column='Population')
print('\nPairwise FST:')
print(pairwise_fst)

# Gene flow (Nm)
nm = popgen.gene_flow(population_column='Population')
print(f'\nGene flow (Nm): {nm:.2f}')
if nm < 1:
    print('  → Limited gene flow')
else:
    print('  → Substantial gene flow')
```

## Site Frequency Spectrum

```python
# Calculate SFS
sfs = popgen.site_frequency_spectrum()
print(f'\nSite Frequency Spectrum: {sfs}')

# Plot SFS
plt.figure(figsize=(10, 6))
plt.bar(range(1, len(sfs) + 1), sfs)
plt.xlabel('Allele Count')
plt.ylabel('Number of Sites')
plt.title('Site Frequency Spectrum')
plt.savefig('sfs.png', dpi=300)

# Segregating sites
s = popgen.segregating_sites()
print(f'\nSegregating sites: {s}')
```

## Network-Based Analysis

```python
# Build network
network = TCSAlgorithm().build_network(alignment)

# Network statistics
net_stats = NetworkStatistics(network)
print(f'\nNetwork properties:')
print(f'  Haplotypes: {net_stats.number_of_nodes()}')
print(f'  Connections: {net_stats.number_of_edges()}')
print(f'  Diameter: {net_stats.diameter()}')

# Visualize with population colors
plot = StaticPlot(network, figsize=(10, 10))
plot.color_by_attribute('Population')
plot.size_by_frequency(min_size=200, max_size=1000)
plot.add_legend(title='Population')
plot.save('population_network.png', dpi=300)
```

## Per-Population Analysis

```python
# Analyze each population separately
populations = metadata['Population'].unique()

pop_results = []
for pop in populations:
    # Subset alignment
    pop_seqs = metadata[metadata['Population'] == pop]['sequence_id'].tolist()
    pop_aln = alignment.subset(pop_seqs)

    # Calculate diversity
    pop_gen = PopulationGenetics(pop_aln)
    pop_results.append(
        {
            'Population': pop,
            'N': len(pop_aln),
            'Haplotypes': pop_aln.n_unique(),
            'Hap_Diversity': pop_gen.haplotype_diversity(),
            'Nuc_Diversity': pop_gen.nucleotide_diversity(),
            'Tajimas_D': pop_gen.tajimas_d(),
        }
    )

# Create summary table
df = pd.DataFrame(pop_results)
print('\nPer-Population Statistics:')
print(df.to_string(index=False))
df.to_csv('population_statistics.csv', index=False)
```

## AMOVA (Analysis of Molecular Variance)

```python
from pypopart.stats import SpatialAnalysis

spatial = SpatialAnalysis(alignment, network)

# AMOVA by population
amova = spatial.amova(group_column='Population')
print('\nAMOVA Results:')
print(f'  Among populations: {amova["among"]:.2f}%')
print(f'  Within populations: {amova["within"]:.2f}%')
print(f'  FST: {amova["fst"]:.4f}')
print(f'  P-value: {amova["p_value"]:.4f}')
```

## Temporal Analysis

```python
# Diversity through time
years = sorted(metadata['Year'].unique())
temporal_diversity = []

for year in years:
    year_seqs = metadata[metadata['Year'] == year]['sequence_id'].tolist()
    year_aln = alignment.subset(year_seqs)
    year_gen = PopulationGenetics(year_aln)

    temporal_diversity.append(
        {
            'Year': year,
            'Haplotype_Diversity': year_gen.haplotype_diversity(),
            'Nucleotide_Diversity': year_gen.nucleotide_diversity(),
        }
    )

# Plot temporal trends
df_temporal = pd.DataFrame(temporal_diversity)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

ax1.plot(df_temporal['Year'], df_temporal['Haplotype_Diversity'], 'o-')
ax1.set_xlabel('Year')
ax1.set_ylabel('Haplotype Diversity')
ax1.set_title('Haplotype Diversity Over Time')

ax2.plot(df_temporal['Year'], df_temporal['Nucleotide_Diversity'], 'o-')
ax2.set_xlabel('Year')
ax2.set_ylabel('Nucleotide Diversity')
ax2.set_title('Nucleotide Diversity Over Time')

plt.tight_layout()
plt.savefig('temporal_diversity.png', dpi=300)
```

## Complete Analysis Report

```python
# Generate comprehensive report
report = {
    'Dataset': {
        'sequences': len(alignment),
        'length': alignment.length,
        'haplotypes': alignment.n_unique(),
    },
    'Diversity': {
        'haplotype': round(popgen.haplotype_diversity(), 4),
        'nucleotide': round(popgen.nucleotide_diversity(), 6),
        'theta': round(popgen.theta_watterson(), 6),
    },
    'Neutrality': {
        'tajimas_d': round(popgen.tajimas_d(), 4),
        'fus_fs': round(popgen.fus_fs(), 4),
    },
    'Structure': {
        'fst': round(popgen.calculate_fst(population_column='Population'), 4),
        'nm': round(popgen.gene_flow(population_column='Population'), 2),
    },
    'Network': {
        'nodes': net_stats.number_of_nodes(),
        'edges': net_stats.number_of_edges(),
        'diameter': net_stats.diameter(),
    },
}

# Save report
import json

with open('popgen_report.json', 'w') as f:
    json.dump(report, f, indent=2)

# Print summary
print('\n' + '=' * 50)
print('POPULATION GENETICS SUMMARY')
print('=' * 50)
for category, values in report.items():
    print(f'\n{category}:')
    for key, value in values.items():
        print(f'  {key}: {value}')
```

## Interpreting Results

### Tajima's D

- **< -2**: Population expansion or purifying selection
- **-2 to +2**: Neutral evolution
- **> +2**: Balancing selection or bottleneck

### FST

- **0.00-0.05**: Little differentiation
- **0.05-0.15**: Moderate differentiation
- **0.15-0.25**: Great differentiation
- **> 0.25**: Very great differentiation

### Gene Flow (Nm)

- **< 1**: Limited gene flow, drift dominates
- **> 1**: Substantial gene flow
- **> 4**: Panmixia

## Next Steps

- [Basic Workflow](basic_workflow.md): General analysis
- [Algorithm Comparison](algorithm_comparison.md): Network algorithms
- [Analysis Guide](../guide/analysis.md): Complete documentation
