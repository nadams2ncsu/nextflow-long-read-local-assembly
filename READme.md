# Nextflow Long-Read Local Assembly

A containerized Nextflow workflow for localized genome assembly from long-read whole-genome sequencing (WGS) data on HPC systems.

## Goal

Genotype structurally complex genomic regions using long-read WGS data through localized assembly and alignment-based characterization of user-defined genomic regions.

## Workflow Overview
The workflow extracts reads overlapping a target region with **samtools**, performs phased local assembly with **hifiasm**, converts assembly graphs to FASTA using a custom **Rust `gfa2fasta` utility**, aligns assembled contigs to a user-supplied reference with **minimap2**, and summarizes assembly contiguity and alignment-based variation using Python. Optional structural haplotype classification can be performed using a user-supplied YAML configuration.

## Container

The `long_read_local_asm.sif` container includes:

- samtools v1.22.1
- hifiasm v0.25.0
- minimap2 v2.30
- custom Rust `gfa2fasta` utility
- Python 3.12
- NumPy v2.3.3
- pandas v2.3.2
- pysam v0.23.3
- PyYAML v6.0.2

**Nextflow and Singularity/Apptainer must be available on the host system.**

## Requirements

The workflow is designed for HPC execution using **SLURM** and accepts PacBio HiFi (`PB`) and Oxford Nanopore (`ONT`) WGS data.

Input BAM files must be aligned to an appropriate reference genome, coordinate sorted, and indexed. The BAM index must be available alongside the BAM file.

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/nadams2ncsu/nextflow-long-read-local-assembly.git
cd nextflow-long-read-local-assembly
```

### 2. Prepare the sample configuration

An example tab-delimited configuration file is provided at `config/samples.tsv`.

```bash
cp config/samples.tsv samples.tsv
```

| Column | Description |
|---|---|
| `sample` | Sample ID |
| `consortium` | Dataset or consortium label (e.g., `LRSC`, `HPRC`, `HGSVC`, or `NA`) |
| `bam` | Path to the aligned, coordinate-sorted, and indexed long-read BAM |
| `gene` | Gene or genomic-region identifier |
| `coordinates` | Target coordinates in `chr#:START-END` format |
| `data_type` | Sequencing platform (`ONT` or `PB`) |
| `flanks` | Flanking sequence length(s), in kb, added to the target for local assembly |

Example:

```text
sample	consortium	bam	gene	coordinates	data_type	flanks
HG00513	LRSC	/path/to/HG00513.bam	FCGR2_3	chr1:161505457-161678654	ONT	50,100,400
HG00097	HPRC	/path/to/HG00097.bam	FCGR2_3	chr1:161505457-161678654	PB	200,400,500
```

Supported default flank sizes are `50, 100, 200, 300, 400, 500, 1000` kb, where `1000` represents 1 Mb. Multiple flank sizes can be provided as a comma-separated list (e.g., `50,100,200,400`), or `all` can be specified to run all supported flank sizes.

Flanking sequence is used for read extraction and local assembly. Assembly contiguity is evaluated only across the original user-supplied target coordinates.

### 3. Provide the reference genome

Provide the path to the reference genome FASTA used for assembled-contig alignment `/path/to/reference/genome.fasta`

### 4. Run the workflow

Make the submission scripts executable:

```bash
chmod +x submit_local_asm_workflow.sh run_nextflow.sh
```

Run local assembly and alignment analysis:

```bash
./submit_local_asm_workflow.sh  samples.tsv /path/to/reference/genome.fasta
```

The submission script launches the Nextflow controller as a SLURM job. Nextflow then manages individual workflow processes through SLURM using the resources defined in `nextflow.config`.

The haplotype configuration is optional. Without it, the workflow performs local assembly and alignment-based characterization without structural haplotype classification.

## Inputs

| Input | Required | Description |
|---|---|---|
| Sample configuration | Yes | Tab-delimited sample, BAM, target-region, sequencing-platform, and flank information |
| Reference genome | Yes | Reference FASTA used for assembled-contig alignment |
| Haplotype configuration | No | YAML file defining known structural haplotypes |

## Assembly Summary

Each hifiasm phased assembly (`hap1` and `hap2`) is evaluated independently across the original user-supplied target coordinates. The summary includes sample, consortium, target gene/region, flank size, hifiasm haplotype, assembly status, contig information, target coverage, SNVs, small insertions/deletions, and large structural events.

By default, a contig must span at least 90% of the target coordinates to be considered contiguous, and insertions or deletions ≥40 kb are reported as large structural events. Modify `modules/summarize.nf` to change percent of the region covered and large structural variant thresholds.

## Outputs

| Output | Description |
|---|---|
| Regional BAM/FASTQ | Long reads extracted from the target plus selected flanking sequence |
| hifiasm GFA | Phased `hap1` and `hap2` assembly graphs |
| Assembly FASTA | Assembled contigs converted from GFA using `gfa2fasta` |
| Alignment BAM | Assembled contigs aligned to the supplied reference genome |
| Summary TSV | Assembly contiguity, coverage, contig, and variant information |
| Haplotype classification | Optional structural haplotype assignment when a configuration is supplied |

Nextflow execution reports, including the trace, timeline, report, and DAG, are written to `results/pipeline_info/`.

## Custom Utilities

**`gfa2fasta`** is a custom Rust utility used to convert hifiasm GFA sequence records to FASTA. Source code is provided in `rust/gfa2fasta/`, and the compiled executable is included in the container.

**`summarize_alignment.py`** evaluates phased assemblies across the user-supplied target coordinates and reports assembly contiguity, coverage, contig information, and alignment-based variation.

**`classify_haplotype.py`** optionally compares observed structural events with known haplotype definitions supplied in a YAML configuration, allowing locus-specific haplotypes to be defined without hard-coding them into the classification program.

## Notes

Successful workflow completion does not necessarily indicate that the target locus was assembled as a single contiguous sequence. Assembly success can depend on sequencing coverage, read alignment, sequence similarity, structural complexity, and the amount of flanking sequence used. 

Assembly contiguity should therefore be evaluated using the generated summary results.

## Citation

If you find this workflow useful in your research or project, please consider citing:

> Adams, Nicole. (2026). *Nextflow Long-Read Local Assembly*. GitHub repository.

https://github.com/nadams2ncsu/nextflow-long-read-local-assembly
