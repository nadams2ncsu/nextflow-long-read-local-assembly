# Nextflow Modules

This directory contains the Nextflow DSL2 process modules that define the individual stages of the long-read local assembly workflow.

Each module encapsulates a specific workflow task, container execution, inputs and outputs, and publication of workflow results. Analysis logic is primarily implemented by command-line tools or the Python utilities in `bin/`, while these modules define how those tasks are executed and connected by Nextflow.

## Modules

| Module | Description |
|---|---|
| `extract_reads.nf` | Extracts reads from the target genomic region and user-specified flanking sequence from each input BAM. |
| `bam2fastq.nf` | Converts extracted regional BAM reads to FASTQ for local assembly. |
| `hifiasm.nf` | Performs phased local assembly with hifiasm and generates haplotype-resolved assembly graphs. |
| `gfa2fasta.nf` | Converts phased hifiasm GFA output to FASTA using the custom Rust `gfa2fasta` utility. |
| `align_contigs.nf` | Aligns assembled contigs to the workflow reference using minimap2. |
| `sort_index.nf` | Sorts and indexes contig-to-reference alignments for downstream assembly evaluation. |
| `coverage_qc.nf` | Evaluates coverage across locus-specific QC intervals in each original input BAM. |
| `combine_qc.nf` | Combines per-dataset coverage QC results into a workflow-level summary. |
| `read_stats.nf` | Calculates regional long-read statistics from each original input BAM. |
| `combine_read_stats.nf` | Combines regional read statistics across input datasets. |
| `summarize.nf` | Evaluates aligned assemblies for target-region coverage, contiguity, sequence characteristics, and assembly-reference differences. |
| `classify_haplotypes.nf` | Optionally assigns known structural haplotypes using a locus-specific haplotype configuration. |
| `combine_results.nf` | Combines individual assembly evaluations into the workflow-level assembly summary. |
| `select_contiguous_pairs.nf` | Identifies sample, group, locus, and flank combinations for which both hap1 and hap2 assemblies are contiguous. |
| `reconstruct_sequence.nf` | Reconstructs final target-region sequences for paired-contiguous assemblies. |
| `note_assembly_differences.nf` | Identifies sequence differences between each reconstructed assembly and the workflow reference. |
| `merge_assembly_differences.nf` | Combines hap1 and hap2 assembly-reference differences for each retained assembly run. |
| `across_groups_comparison.nf` | Compares independently generated assemblies when the same sample and locus are represented by multiple input groups. |
| `compare_group_differences.nf` | Evaluates concordance of assembly-reference differences across independently generated assemblies. |
| `generate_report.nf` | Generates the final interactive HTML report from workflow QC, assembly, and comparison summaries. |

## Workflow Organization

The modules form three related analysis paths.

### Local assembly and evaluation

Input BAMs are expanded across the requested flank sizes and processed through the primary assembly path:

`EXTRACT_READS` → `BAM2FASTQ` → `HIFIASM` → `GFA2FASTA` → `ALIGN_CONTIGS` → `SORT_INDEX` → `SUMMARIZE`

When a structural haplotype configuration is supplied, assembly summaries are additionally processed by `CLASSIFY_HAPLOTYPES`.

### Input quality control

Coverage and read-level QC are calculated directly from the original input BAMs:

`COVERAGE_QC` → `COMBINE_QC`

`READ_STATS` → `COMBINE_READ_STATS`

These analyses are independent of the sample × flank assembly expansion.

### Final assembly characterization

Assembly runs are retained for final reconstruction only when both phased assemblies satisfy the contiguity requirement:

`COMBINE_RESULTS` → `SELECT_CONTIGUOUS_PAIRS` → `RECONSTRUCT_SEQUENCE`

Retained assemblies are then characterized relative to the workflow reference and, when multiple datasets are available for the same sample, compared across groups:

`RECONSTRUCT_SEQUENCE` → `NOTE_ASSEMBLY_DIFFERENCES` → `MERGE_ASSEMBLY_DIFFERENCES`

`SELECT_CONTIGUOUS_PAIRS` → `ACROSS_GROUPS_COMPARISON`

The resulting QC, assembly, and comparison summaries are combined by `GENERATE_REPORT` to produce the final HTML report.

## Module Design

Modules are designed to keep workflow orchestration separate from analysis implementation:

- **Nextflow modules** define process inputs and outputs, resource requests, container execution, and data flow.
- **Python utilities in `bin/`** perform QC, assembly evaluation, sequence reconstruction, comparison, and report generation.
- **`gfa2fasta` in `rust/`** provides the custom GFA-to-FASTA conversion used after hifiasm assembly.
- **External command-line tools** including samtools, hifiasm, and minimap2 perform read processing, assembly, and alignment.

The modules are imported and connected in `main.nf`, which controls sample × flank expansion, optional analysis branches, paired-contiguous filtering, and downstream aggregation.

## Execution

The modules are components of the complete workflow and are not intended to be launched individually.

Run the workflow from the repository root with:

```bash
./submit_local_asm_workflow.sh config/samples.tsv /path/to/reference/genome.fasta
```

The submission script launches the Nextflow controller through SLURM. Nextflow then schedules individual module processes according to the resource definitions in `nextflow.config`.

## Software Environment

All workflow modules execute using the Apptainer/Singularity image `long_read_local_asm.sif`, providing a consistent software environment across compute nodes.

Process-specific CPU and memory requirements are defined in `nextflow.config`.
