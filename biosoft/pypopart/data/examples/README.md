# Example datasets

## `sample.fasta`

Five sequences of 16 bp collapsing to three haplotypes. Small enough to
reason about by hand, which makes it useful for a quick smoke test of the
CLI or the GUI.

## `demo_5pop.fasta` + `demo_5pop_metadata.csv`

Seventy simulated sequences of 100 bp drawn from five populations
(`Coastal`, `Highland`, `Riverine`, `Desert`, `Island`), fourteen each,
collapsing to 43 unique haplotypes.

Each population descends from its own founder, three mutations off a
shared ancestor, and roughly half the individuals in a population carry
private variation on top. The result is five recognisable clusters joined
by a few mutational steps — enough structure to exercise population
colouring, pie-chart nodes, the legend and the figure exports at a size
where layout and label crowding actually show up.

Load the pair together in the GUI: the FASTA as the sequence file and the
CSV as the metadata file.
