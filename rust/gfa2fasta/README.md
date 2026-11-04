# gfa2fasta

A lightweight Rust utility for converting hifiasm GFA assembly output to FASTA format.

The utility extracts sequence-containing segment (`S`) records from the GFA and writes the corresponding contig sequences as FASTA records.

Within the long-read local assembly workflow, `gfa2fasta` converts the phased hifiasm assembly graphs for downstream contig alignment and assembly evaluation.
