# Test Dataset

This directory contains a downsampled dataset for testing the long-read local assembly workflow using the human **FCGR2/3** region.

The test data are the same long-read WGS samples used during workflow development. Input BAMs contain reads from a limited genomic region surrounding FCGR2/3, allowing the workflow to be tested without requiring complete WGS datasets.

## Directory Structure
- `config/` contains the sample configuration used for the test workflow.
- `data/` contains downsampled long-read BAM files and their indexes.
- `reference/` contains chromosome 1 from the GRCh38 reference and its index.
- `logs/` contains example workflow execution logs.
- `test_results/` contains representative outputs from a completed test run.

## Running the Test

From the repository root:

```bash
./submit_local_asm_workflow.sh test/config/samples.tsv test/reference/GRCh38_chr1.fasta
```

The workflow uses the same Nextflow processes, container, QC procedures, assembly evaluation, and reporting steps as a full analysis.

*Direct paths to all input files, QC and haplotype yaml files need to be updated for each user in the `samples.tsv`.*

## Test Data

The downsampled BAM files retain reads surrounding the FCGR2/3 region rather than complete genome-wide sequencing data. They are intended to provide a smaller dataset for evaluating workflow execution and representative outputs.

Because these BAMs contain reads only around the target region, genome-wide coverage metrics are expected to be low and may trigger coverage QC warnings. These warnings are expected for the reduced test dataset.

## Expected Results

Representative outputs from the test workflow are provided in `test_results/`.

These results can be used to inspect the expected output structure and confirm successful execution of the major workflow stages, including local assembly, assembly evaluation, sequence reconstruction, assembly-reference comparison, cross-group comparison, and HTML report generation.
