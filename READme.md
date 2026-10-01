# Nextflow Long-Read Local Assembly

A Nextflow DSL2 workflow for local, haplotype-resolved assembly of long-read whole-genome sequencing (WGS) data on SLURM-based HPC systems.

The workflow performs regional read extraction, phased assembly, reference-based assembly evaluation, *optional* structural haplotyping, sequence reconstruction, cross-dataset comparison, and automated reporting. PacBio HiFi (`PB`) and Oxford Nanopore (`ONT`) WGS data are supported.

---

## Goal

Provide a reproducible and configurable workflow for local assembly and characterization of structurally complex genomic regions from long-read WGS data. The human **FCGR2/3** locus is included as a test case.

---

## Workflow

![Nextflow Long-Read Local Assembly workflow](images/metro_map.png)

Only assembly runs in which **both hap1 and hap2 are contiguous** are retained for final target-region sequence reconstruction. When the same sample is represented by multiple datasets, the resulting assemblies can also be compared across groups.

---

## Implementation

The workflow is implemented in **Nextflow DSL2** for execution with **SLURM** and **Apptainer/Singularity**. Analysis steps use **samtools**, **hifiasm**, **minimap2**, **Python**, and a custom **Rust `gfa2fasta`** utility. Analysis dependencies are packaged in `long_read_local_asm.sif`.

### Requirements

- **Nextflow v25.10.3**
- **Apptainer/Singularity** (tested with Apptainer v1.4.2-1)
- **SLURM**
- coordinate-sorted and indexed long-read WGS BAM files
- indexed reference genome FASTA

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/nadams2ncsu/nextflow-long-read-local-assembly.git
cd nextflow-long-read-local-assembly
```

### 2. Configure samples

Edit the example sample configuration at `config/samples.tsv`.

The sample table contains one row for each input dataset and genomic region.

| Column | Description |
|---|---|
| `sample` | Sample ID |
| `group` | Dataset, consortium, sequencing source, batch, or comparison-group label |
| `bam` | Path to coordinate-sorted and indexed long-read BAM |
| `gene` | Gene or genomic-region identifier |
| `coordinates` | Target coordinates in the coordinate system of the input BAM (`chr#:START-END`) |
| `data_type` | Sequencing data type: `ONT` or `PB` |
| `flanks` | Flanking-region size(s), in kb |
| `qc_config` | Path to locus-specific coverage QC YAML |
| `haplotype_config` | Path to optional structural haplotype YAML, or `NA` |

Supported flank sizes are `50, 100, 200, 300, 400, 500, 1000` kb.

Multiple flank sizes can be supplied as a comma-separated list such as `100,200,400`, or use `all` to evaluate all supported flank sizes.

### 3. Configure locus-specific QC

Coverage QC configurations are provided in `config/qc/`.

Each YAML defines the genomic intervals used to evaluate coverage for a target locus. Different coordinate systems can be associated with different input groups when the original BAMs were aligned to different references.

### 4. Optional structural haplotyping

Known structural haplotypes can be defined in `config/haplotypes/`.

These configurations describe the target evaluation region and expected structural events used for locus-specific haplotype classification.

For regions where structural haplotyping is not required, set `haplotype_config` to `NA` in `config/samples.tsv`.

### 5. Run

Make submission scripts as well as process python scripts and modules executable:

```bash
chmod +x submit_local_asm_workflow.sh run_nextflow.sh
chmod +x bin/*
chmod +x modules/*
```

Submit the workflow:

```bash
./submit_local_asm_workflow.sh config/samples.tsv /path/to/reference/genome.fasta
```

The submission script launches the Nextflow controller as a SLURM job. Nextflow manages sample × flank expansion and submits individual workflow processes to SLURM.

---

## Outputs

Final workflow outputs are written to `results/`.

| Output | Description |
|---|---|
| **Coverage QC** | Coverage measurements across configured locus-specific intervals |
| **Read statistics** | Regional long-read statistics for each input dataset |
| **Assembly summary** | Assembly status, sequence characteristics, and optional structural haplotype assignments |
| **Reconstructed FASTA** | Final hap1 and hap2 target-region sequences from paired-contiguous assemblies |
| **Assembly-reference differences** | Sequence differences between paired-contiguous assemblies and the workflow reference |
| **Across-group comparisons** | Sequence and assembly-reference difference concordance between independently generated assemblies |
| **HTML report** | Interactive summary of assembly outcomes, QC, sequence comparisons, and assembly-reference difference concordance |

Representative workflow outputs are provided in `example/results/`. These files demonstrate the final output structure without requiring execution of the complete workflow.

Intermediate files generated during execution are maintained by Nextflow in the `work/` directory and are not included with the example results.

Nextflow execution metadata, including the DAG, trace, timeline, and execution report, are written to `results/pipeline_info/`.

---

## Test Dataset

A reduced test dataset is provided in `test/`.

The test dataset contains reduced versions of the same representative samples used to generate the outputs in `example/results/`. It is intended to exercise the complete workflow with smaller input files.

### Test inputs

Downsampled long-read BAM files are provided in `test/data/`. These BAMs retain reads surrounding the target locus while removing reads from the remainder of the genome, reducing input size while preserving the reads required for local assembly.

The reference required for the test workflow is provided in `test/reference/`.

Test-specific sample, QC, and haplotype configurations are provided in `test/config/`.

### Run the test dataset

From the repository root:

```bash
./submit_local_asm_workflow.sh test/config/samples.tsv test/reference/GRCh38_chr1.fasta
```

The test dataset exercises the same assembly, evaluation, reconstruction, comparison, and reporting steps as the full workflow.

Because the test BAMs contain reads only around the FCGR2/3 region, genome-wide coverage metrics are expected to be low and may trigger coverage QC warnings. These warnings are expected for the reduced test dataset.

The resulting outputs can be compared with the representative results in `example/results/`.

---

## Configuration Notes

### Input coordinates

The `coordinates` column in `samples.tsv` describes the target region in the coordinate system of the **original input BAM**.

These coordinates are used for regional read extraction and input read statistics. Assembled contigs are subsequently aligned to the reference genome supplied when the workflow is launched.

For loci with a structural haplotype configuration, `evaluation_coordinates` in the haplotype YAML define the target interval used for downstream assembly evaluation and sequence reconstruction.

This allows input datasets aligned to different reference assemblies to be processed while downstream assemblies are evaluated against a common workflow reference.

### Paired-contiguous assemblies

Final sequence reconstruction requires both phased assemblies for a `sample + group + gene + flank` combination to be classified as **contiguous**.

If either haplotype fails this requirement, neither haplotype from that assembly run is retained as a final reconstructed assembly.

Samples represented by only one group can still produce final reconstructed sequences. Multiple groups are required only for across-group comparison.

---

## Citation

If you use this workflow in your research or project, please consider citing:

> Adams, Nicole. (2026). *Nextflow Long-Read Local Assembly*. GitHub repository.

https://github.com/nadams2ncsu/nextflow-long-read-local-assembly
