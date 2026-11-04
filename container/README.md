# Container

This directory contains the Apptainer/Singularity definition file used to build the software environment for the long-read local assembly workflow.

`nextflow_local_asm.def` uses a multi-stage build to compile and install the required bioinformatics tools, custom Rust `gfa2fasta` utility, and Python analysis environment.

## Included Software

- Samtools 1.22.1
- Hifiasm 0.25.0
- Minimap2 2.30
- Python 3
- `gfa2fasta` 0.1.0
- numpy 2.3.3
- pandas 2.3.2
- pysam 0.23.3
- PyYAML 6.0.2
- Plotly 6.3.0
