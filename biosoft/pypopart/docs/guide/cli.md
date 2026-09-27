# Command-Line Interface

PyPopART installs a `pypopart` command with five subcommands. Every
example below matches the shipped CLI (`pypopart --help` is always the
authoritative reference).

```text
Usage: pypopart [OPTIONS] COMMAND [ARGS]...

Options:
  --version      Show the version and exit.
  -v, --verbose  Increase verbosity (can be repeated: -v, -vv, -vvv)
  -q, --quiet    Suppress all output except errors

Commands:
  analyze    Analyze haplotype network statistics.
  info       Display information about PyPopART capabilities.
  load       Load and validate sequence alignment data.
  network    Construct haplotype network from sequence alignment.
  visualize  Visualize haplotype network.
```

Use `-v` to print full tracebacks on errors, and `-q` to suppress all
informational output (useful in scripts).

## load — inspect an alignment

```bash
pypopart load sequences.fasta
pypopart load sequences.fasta -m metadata.csv       # match metadata rows
pypopart load sequences.nex -o resaved.fasta        # convert formats
```

| Option                                         | Description                                                |
| ---------------------------------------------- | ---------------------------------------------------------- |
| `-f, --format [fasta\|nexus\|phylip\|genbank]` | Input format (auto-detected from the extension if omitted) |
| `-m, --metadata FILE`                          | Metadata CSV; reports how many sequences matched           |
| `-o, --output PATH`                            | Re-save the alignment (format from the output extension)   |

Prints alignment statistics: sequence count, length, variable sites,
parsimony-informative sites, and GC content.

## network — build a haplotype network

```bash
pypopart network sequences.fasta -o network.graphml            # MJN (default)
pypopart network sequences.fasta -a tcs -o network.graphml     # TCS
pypopart network sequences.fasta -a mst -d k2p -o net.graphml  # K2P distances
pypopart network sequences.fasta -a pn --seed 42 -o net.graphml
```

| Option                                              | Description                                                                      |
| --------------------------------------------------- | -------------------------------------------------------------------------------- |
| `-a, --algorithm [mst\|msn\|tcs\|mjn\|pn\|tsw]`     | Construction algorithm (default: `mjn`)                                          |
| `-d, --distance [hamming\|jc\|k2p\|tn\|tamura_nei]` | Distance metric (default: `hamming`; `tamura_nei` is an alias for `tn`)          |
| `-e, --epsilon FLOAT`                               | Epsilon for MSN/MJN (default: 0)                                                 |
| `-p, --parsimony-limit FLOAT`                       | Opt in to a TCS connection limit at this confidence (PopART default is no limit) |
| `--seed INTEGER`                                    | Random seed for stochastic algorithms (`pn`)                                     |
| `-o, --output PATH`                                 | Output network file                                                              |
| `--format [graphml\|gml\|json\|nexus]`              | Output format (default: `graphml`)                                               |

## analyze — network statistics

```bash
pypopart analyze network.graphml                 # statistics (default)
pypopart analyze network.graphml --topology
pypopart analyze network.graphml --popgen -a sequences.fasta
pypopart analyze network.graphml --stats -o results.json
```

| Option                 | Description                                                        |
| ---------------------- | ------------------------------------------------------------------ |
| `--stats`              | Network statistics (the default when no flag is given)             |
| `--topology`           | Topology analysis: components, star patterns, ancestral candidates |
| `--popgen`             | Population genetics (Tajima's D); requires `-a/--alignment`        |
| `-a, --alignment FILE` | The original alignment (needed for `--popgen` and diversity)       |
| `-o, --output PATH`    | Write results as JSON                                              |

## visualize — plot a network

```bash
pypopart visualize network.graphml -o network.png --show-labels
pypopart visualize network.graphml -o network.html            # interactive
pypopart visualize network.graphml -o net.pdf --layout kamada_kawai
```

| Option                                      | Description                                                  |
| ------------------------------------------- | ------------------------------------------------------------ |
| `-o, --output PATH`                         | Output file; `.html` produces an interactive plot (required) |
| `--layout [spring\|circular\|kamada_kawai]` | Layout algorithm (default: `spring`)                         |
| `--width / --height INTEGER`                | Figure size in pixels (defaults: 800x600)                    |
| `--interactive`                             | Force interactive HTML output                                |
| `--show-labels`                             | Show node labels                                             |

Static output (PNG/PDF/SVG) requires the `viz` extra:
`pip install 'pypopart[viz]'`.

## info — capability listings

```bash
pypopart info --list-algorithms
pypopart info --list-distances
pypopart info --list-formats
```

## A complete pipeline

```bash
pypopart load data/sequences.fasta -m data/populations.csv
pypopart network data/sequences.fasta -a mjn -o results/network.graphml
pypopart analyze results/network.graphml --stats --topology -o results/stats.json
pypopart visualize results/network.graphml -o results/network.png --show-labels
```
