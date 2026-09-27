# Network Layout Algorithms

This guide helps you choose the best layout algorithm for your haplotype network visualization in PyPopART.

## Quick Selection Guide

### By Network Size

| Network Size          | Recommended Algorithm | Alternative    |
| --------------------- | --------------------- | -------------- |
| Small (<50 nodes)     | Kamada-Kawai          | Force-Directed |
| Medium (50-500 nodes) | Force-Directed        | Spectral       |
| Large (>500 nodes)    | Spectral              | Hierarchical   |

### By Purpose

| Purpose                | Algorithm              | Notes                               |
| ---------------------- | ---------------------- | ----------------------------------- |
| General visualization  | Force-Directed         | Good balance of speed and quality   |
| Fast preview           | Hierarchical           | Instant layout, tree structure      |
| Large datasets         | Spectral               | Maintains structure, very fast      |
| Highest quality        | Kamada-Kawai           | Slow but optimal for small networks |
| Highlight central node | Radial                 | Places important node at center     |
| Simple connectivity    | Circular               | Shows connection patterns clearly   |
| Proportional distances | Spring/KK Proportional | Edge length reflects mutation count |

## Algorithm Details

### 1. Hierarchical Layout

**Best for**: Fast previews, tree-like structures

**Speed**: ⚡⚡⚡⚡⚡ Fastest (0.14ms for 100 nodes)

**Quality**: ⭐⭐⭐ Good for hierarchical data

**Description**: Arranges nodes in levels based on distance from a root node. Creates a tree-like structure that's easy to interpret.

**Advantages**:

- Extremely fast, works well for very large networks
- Clear hierarchical relationships
- Handles disconnected components

**Disadvantages**:

- May not show cyclical relationships well
- Layout depends on choice of root node

**Parameters**:

```python
layout = manager.compute_layout(
    'hierarchical',
    root_node='H1',  # Optional: specify root
    vertical=True,  # True = top-down, False = left-right
    width=2.0,  # Layout width
    height=2.0,  # Layout height
)
```

### 2. Spectral Layout

**Best for**: Large networks (100-1000+ nodes)

**Speed**: ⚡⚡⚡⚡ Very fast (6.5ms for 100 nodes)

**Quality**: ⭐⭐⭐⭐ Excellent structure preservation

**Description**: Uses eigenvectors of the graph Laplacian to position nodes. Excellent balance of speed and quality for large networks.

**Advantages**:

- Much faster than force-directed or Kamada-Kawai
- Reveals clustering and community structure
- Good for networks with clear groups

**Disadvantages**:

- May produce less aesthetic layouts than force-directed
- Requires connected graph (handles disconnected with fallback)

**Parameters**:

```python
layout = manager.compute_layout(
    'spectral',
    scale=1.0,  # Scale factor
    center=(0.0, 0.0),  # Center position
)
```

### 3. Force-Directed (Spring) Layout

**Best for**: General purpose, medium networks

**Speed**: ⚡⚡⚡ Moderate (26ms for 100 nodes)

**Quality**: ⭐⭐⭐⭐ Very good aesthetic quality

**Description**: Simulates physical springs between connected nodes. The default choice for most visualizations.

**Advantages**:

- Aesthetically pleasing layouts
- Works well for most network types
- Tunable with iterations parameter

**Disadvantages**:

- Slower than spectral or hierarchical
- Can be slow for large networks (>500 nodes)
- Non-deterministic without seed

**Parameters**:

```python
layout = manager.compute_layout(
    'spring',
    k=None,  # Optimal distance (None = auto)
    iterations=50,  # More = better quality but slower
    seed=42,  # For reproducibility
)
```

**Tuning Tips**:

- Increase `iterations` for better quality (try 100-200)
- Decrease `k` to bring nodes closer together
- Use `seed` for consistent layouts across runs

### 4. Kamada-Kawai Layout

**Best for**: Small networks requiring highest quality

**Speed**: ⚡ Slow (187ms for 100 nodes)

**Quality**: ⭐⭐⭐⭐⭐ Optimal quality

**Description**: Minimizes stress based on graph-theoretic distances. Produces optimal layouts but is computationally expensive.

**Advantages**:

- Highest quality layouts
- Respects graph distances precisely
- Deterministic results

**Disadvantages**:

- Very slow for networks >100 nodes
- O(N³) time complexity
- Not suitable for interactive use with large networks

**Parameters**:

```python
layout = manager.compute_layout(
    'kamada_kawai',
    scale=1.0,  # Scale factor
    center=(0.0, 0.0),  # Center position
)
```

**Warning**: Only use for small networks (<50 nodes) or when layout quality is absolutely critical.

### 5. Circular Layout

**Best for**: Showing connectivity patterns

**Speed**: ⚡⚡⚡⚡⚡ Very fast (0.16ms for 100 nodes)

**Quality**: ⭐⭐ Simple but effective

**Description**: Arranges nodes evenly spaced around a circle.

**Advantages**:

- Extremely fast
- Good for comparing edge densities
- Works well with node coloring

**Disadvantages**:

- Doesn't reflect graph structure
- Can be cluttered for dense networks

**Parameters**:

```python
layout = manager.compute_layout(
    'circular',
    scale=1.0,  # Radius of circle
    center=(0.0, 0.0),  # Center position
)
```

### 6. Radial Layout

**Best for**: Networks with clear central node

**Speed**: ⚡⚡⚡⚡ Fast (0.26ms for 100 nodes)

**Quality**: ⭐⭐⭐ Good for star-like networks

**Description**: Places a central node at origin, others in concentric rings based on distance.

**Advantages**:

- Fast and intuitive
- Highlights central hub
- Shows distance from center clearly

**Disadvantages**:

- Works best for star-like topologies
- May be cluttered if many nodes at same distance

**Parameters**:

```python
layout = manager.compute_layout(
    'radial',
    center_node='H1',  # Node to place at center
    scale=1.0,  # Scale factor
)
```

### 7. Proportional Edge Length Layouts

**Best for**: Showing mutation distances accurately

**Speed**: ⚡⚡⚡ Moderate (similar to base algorithms)

**Quality**: ⭐⭐⭐⭐ Accurate distance representation

**Description**: Variations of spring and Kamada-Kawai layouts where edge lengths are proportional to the number of mutations between haplotypes.

**Advantages**:

- Visually represents genetic distance
- Helps identify distant vs. close relationships
- Useful for understanding evolutionary distances

**Disadvantages**:

- Can produce larger, more spread-out layouts
- May reduce aesthetic appeal for highly variable data

**Available Variations**:

#### Spring Proportional

```python
layout = manager.compute_layout('spring_proportional', iterations=100, seed=42)
```

#### Kamada-Kawai Proportional

```python
layout = manager.compute_layout('kamada_kawai_proportional', scale=1.0)
```

## Performance Benchmarks

All benchmarks run on standard hardware with Python 3.12, NetworkX 3.5.

### Runtime Comparison (100 nodes)

```
Hierarchical:     0.14 ms  ⚡⚡⚡⚡⚡
Circular:         0.16 ms  ⚡⚡⚡⚡⚡
Radial:           0.26 ms  ⚡⚡⚡⚡⚡
Spectral:         6.50 ms  ⚡⚡⚡⚡
Force-Directed:  26.00 ms  ⚡⚡⚡
Kamada-Kawai:   187.00 ms  ⚡
```

### Scalability (200 nodes)

```
Hierarchical:      0.28 ms
Spectral:         23.88 ms
Force-Directed:   87.00 ms
Kamada-Kawai:   1202.00 ms  ⚠️ Not recommended
```

## Best Practices

### 1. Start Fast, Refine Later

For exploratory analysis:

1. Use **Hierarchical** for quick preview
2. Switch to **Force-Directed** for publication
3. Try **Spectral** if force-directed is too slow

### 2. Consider Network Properties

- **Clustered networks**: Use Spectral
- **Hierarchical networks**: Use Hierarchical
- **Star topologies**: Use Radial
- **Dense networks**: Use Circular or Spectral

### 3. Optimization Tips

**For large networks**:

- Reduce force-directed iterations: `iterations=30`
- Use spectral instead of force-directed
- Consider hierarchical for fastest results

**For highest quality**:

- Increase iterations: `iterations=100-200`
- Use Kamada-Kawai for small networks
- Set consistent seed for reproducibility

### 4. Interactive Workflows

```python
# Quick preview
layout = manager.compute_layout('hierarchical')

# If network is small (<100 nodes), upgrade to force-directed
if len(network.nodes) < 100:
    layout = manager.compute_layout('spring', iterations=100, seed=42)

# For large networks, use spectral
if len(network.nodes) > 500:
    layout = manager.compute_layout('spectral')
```

## Common Issues and Solutions

### Layout looks cluttered

**Problem**: Nodes overlap or edges cross excessively

**Solutions**:

- Increase iterations for force-directed: `iterations=100`
- Try spectral layout for better spacing
- Use radial or circular for clearer structure

### Layout is too slow

**Problem**: Layout computation takes too long

**Solutions**:

- Switch to spectral or hierarchical
- Reduce iterations: `iterations=30`
- Use hierarchical for instant preview

### Disconnected components overlap

**Problem**: Separate network components are positioned on top of each other

**Solutions**:

- Use force-directed or spectral (handle components automatically)
- Manually adjust with ManualLayout after initial computation

### Layout changes every time

**Problem**: Non-deterministic layouts make comparison difficult

**Solutions**:

- Set seed for force-directed: `seed=42`
- Use deterministic algorithms: Hierarchical, Kamada-Kawai, Circular, Radial

## Advanced Usage

### Combining Layouts

You can start with one algorithm and refine with another:

```python
# Start with fast hierarchical
initial = manager.compute_layout('hierarchical')

# Refine with force-directed (fewer iterations)
refined = manager.compute_layout('spring', iterations=30, pos=initial)
```

### Manual Adjustments

After automatic layout, use ManualLayout for fine-tuning:

```python
# Compute initial layout
positions = manager.compute_layout('spring')

# Create manual layout for adjustments
manual = ManualLayout(network, initial_positions=positions)

# Move a specific node
manual.set_position('H5', (1.0, 2.0))

# Get final positions
final = manual.compute()
```

### Saving and Loading Layouts

Save time by reusing layouts:

```python
# Compute and save
layout = manager.compute_layout('spring', iterations=200, seed=42)
manager.save_layout(layout, 'my_network_layout.json')

# Load later
layout = manager.load_layout('my_network_layout.json')
```

## References

- **NetworkX Layouts**: https://networkx.org/documentation/stable/reference/drawing.html
- **Fruchterman-Reingold**: "Graph Drawing by Force-directed Placement" (1991)
- **Kamada-Kawai**: "An Algorithm for Drawing General Undirected Graphs" (1989)
- **Spectral Layout**: Uses Laplacian eigendecomposition (Koren 2005)

## Getting Help

If you're unsure which layout to use:

1. Try **Hierarchical** first - it's instant and often good enough
2. If that doesn't work well, try **Force-Directed** with default settings
3. For large networks (>500 nodes), use **Spectral**
4. To show mutation distances accurately, try **Spring Proportional** or **Kamada-Kawai Proportional**

For questions or issues, please see the main PyPopART documentation or open an issue on GitHub.
