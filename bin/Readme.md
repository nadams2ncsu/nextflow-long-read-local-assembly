# Workflow Utilities

This directory contains Python utilities used by the Nextflow processes in `modules/` for assembly evaluation, quality control, sequence reconstruction, cross-dataset comparison, and report generation.

These scripts implement analysis logic that is kept separate from the Nextflow workflow definition. Nextflow manages data flow and execution, while the utilities in this directory perform the underlying calculations and file transformations.

## Utilities

| Script | Description |
|---|---|
| `coverage_qc.py` | Calculates coverage across locus-specific intervals defined in the QC configuration and evaluates regional coverage thresholds. |
| `combine_qc.py` | Combines per-sample coverage QC results into the workflow-level coverage QC summary. |
| `read_stats.py` | Calculates long-read statistics for the target region of each input BAM. |
| `combine_read_stats.py` | Combines regional read statistics across input datasets. |
| `summarize_alignment.py` | Evaluates assembled contigs after reference alignment, including target-region coverage, assembly contiguity, sequence length, SNVs, and small and large indels. |
| `classify_haplotypes.py` | Optionally assigns known structural haplotypes using expected structural events defined in a locus-specific haplotype YAML. |
| `combine_results.py` | Combines per-assembly evaluation results into the workflow-level assembly summary. |
| `select_contiguous_pairs.py` | Identifies assembly runs in which both hap1 and hap2 satisfy the contiguity requirement for downstream reconstruction. |
| `reconstruct_sequence.py` | Reconstructs the target-region assembly sequence from the contig-to-reference alignment while retaining assembly insertions within the target interval. |
| `note_assembly_differences.py` | Identifies SNVs, insertions, and deletions between a reconstructed assembly and the workflow reference within the target region. |
| `merge_assembly_differences.py` | Merges hap1 and hap2 assembly-reference differences for each paired-contiguous assembly. |
| `across_groups_comparison.py` | Compares independently generated assemblies for samples represented in multiple input groups and evaluates haplotype correspondence and sequence similarity. |
| `compare_group_differences.py` | Compares assembly-reference differences across groups and summarizes concordant and group-specific differences. |
| `generate_report.py` | Generates the final interactive HTML report from workflow QC, assembly, and cross-group comparison summaries. |

## Configuration

Several utilities use information supplied through the workflow configuration files:

- `config/samples.tsv` defines input datasets, target regions, sequencing data types, flank sizes, and QC and haplotype configuration files.
- `config/qc/` contains locus-specific intervals used for coverage QC.
- `config/haplotypes/` optionally defines evaluation coordinates and expected structural events for structural haplotype classification.

The `coordinates` field in `samples.tsv` describes the target region in the coordinate system of the original input BAM and is used for regional read extraction and input read statistics.

When a haplotype configuration is supplied, its `evaluation_coordinates` define the interval used for downstream assembly evaluation, sequence reconstruction, and assembly-reference difference analysis against the workflow reference.

## Workflow Integration

The utilities in this directory are called by the corresponding Nextflow processes in `modules/`. They are not intended to replace the main workflow entry point.

For normal analysis, launch the complete workflow from the repository root using:

```bash
./submit_local_asm_workflow.sh config/samples.tsv /path/to/reference/genome.fasta
```

Nextflow manages dependencies between utilities, parallel execution, intermediate files, and publication of final results.

## Dependencies

Python dependencies required by these utilities are included in the workflow Apptainer/Singularity image. Major Python packages include:

- `pandas v2.3.2`
- `numpy v2.3.3`
- `pysam v0.23.3`
- `PyYAML v6.0.2`
- `plotly v6.3.0`

See the repository container definition for the complete software environment and version information.
