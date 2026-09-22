# Nextflow Long-Read Local Assembly

A Nextflow workflow for localized genome assembly from long-read whole-genome sequencing (WGS) data on HPC systems.

## Goal

Genotype structurally complex genomic regions using long-read WGS data.

## Workflow

The workflow performs:

**Coverage QC → Read extraction → Local assembly → GFA-to-FASTA conversion → Reference alignment → Assembly characterization → Optional structural haplotyping**

Core tools include:

- **samtools** — coverage QC and regional read extraction
- **hifiasm** — phased local assembly
- **Rust `gfa2fasta`** — assembly graph conversion
- **minimap2** — assembled-contig alignment
- **Python** — assembly characterization and optional structural haplotyping
- **Nextflow + SLURM** — workflow orchestration and HPC execution

PacBio HiFi (`PB`) and Oxford Nanopore (`ONT`) WGS data are supported.

## Requirements

- Nextflow
- Singularity/Apptainer
- SLURM
- Coordinate-sorted and indexed long-read WGS BAM files

All analysis dependencies are provided in `long_read_local_asm.sif`.

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/nadams2ncsu/nextflow-long-read-local-assembly.git
```
```bash
cd nextflow-long-read-local-assembly
```

### 2. Configure samples

Edit the example configuration `config/samples.tsv`

| Column | Description |
|---|---|
| `sample` | Sample ID |
| `consortium` | Consoritum, dataset, batch, or timepoint label |
| `bam` | Path to aligned and indexed long-read BAM |
| `gene` | Gene or genomic-region identifier |
| `coordinates` | Target coordinates (`chr#:START-END`) |
| `data_type` | `ONT` or `PB` |
| `flanks` | Flanking sequence length(s) in kb |

Supported flank sizes are `50, 100, 200, 300, 400, 500, 1000` kb. Multiple sizes can be provided as a comma-separated list, or use `all` to run all supported sizes.

### 3. Configure coverage QC

Provide a YAML file defining the genomic intervals used for coverage QC. Example configurations are available in `config/qc/`

### 4. Run

Make the submission scripts executable:

```bash
chmod +x submit_local_asm_workflow.sh run_nextflow.sh
```

Run the workflow:

```bash
./submit_local_asm_workflow.sh /submit_local_asm_workflow.sh config/samples.tsv  /path/to/reference/genome.fasta
```

Optionally provide a structural haplotype configuration:

```bash
./submit_local_asm_workflow.sh config/samples.tsv /path/to/reference/genome.fasta config/qc/qc_config.yaml config/haplotypes/haplotypes.yaml
```

Bash script submits the nextflow workflow as an array of jobs by sample to a HPC using the SLURM job scheduler.

## Inputs

| Input | Required | Description |
|---|---|---|
| `Sample configuration` | Yes | Sample, BAM, target-region, sequencing-platform, and flank information |
| `Reference genome` | Yes | FASTA used for assembled-contig alignment |
| `QC configuration` | Yes | YAML defining intervals used for coverage QC |
| `Haplotype configuration` | No | YAML defining known structural haplotypes |

## Outputs

Final outputs are written to `results/`.

| Output | Description |
|---|---|
| `Coverage QC` | Whole-genome and target-region coverage metrics |
| `Assembly FASTA` | Phased locally assembled contigs |
| `Alignment BAM`| Sorted and indexed assembled-contig alignments |
| `Summary TSV` | Combined assembly and variant characterization |
| `Haplotype classification` | Optional structural haplotype assignments |

Nextflow trace, timeline, report, and DAG files are written to `results/pipeline_info/`.

Additional documentation and examples are provided within the relevant subdirectories.

## Citation

If you find this workflow useful in your research or project, please consider citing:

> Adams, Nicole. (2026). *Nextflow Long-Read Local Assembly*. GitHub repository.

https://github.com/nadams2ncsu/nextflow-long-read-local-assembly
