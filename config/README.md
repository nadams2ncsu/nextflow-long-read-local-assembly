# Workflow Configuration

This directory contains the input and locus-specific configuration files used by the long-read local assembly workflow.

## Sample Configuration

`samples.tsv` defines the input datasets and regions to analyze. It contains one row for each input dataset and genomic region.

| Column | Description |
|---|---|
| `sample` | Sample identifier |
| `group` | Dataset, consortium, sequencing source, batch, or other comparison-group label |
| `bam` | Path to the coordinate-sorted and indexed long-read WGS BAM |
| `gene` | Gene or genomic-region identifier |
| `coordinates` | Target coordinates in the coordinate system of the original input BAM (`chr#:START-END`) |
| `data_type` | Sequencing data type: `ONT` or `PB` |
| `flanks` | One or more flanking-region sizes used for local assembly, in kb |
| `qc_config` | Path to the locus-specific coverage QC YAML |
| `haplotype_config` | Path to an optional structural haplotype YAML, or `NA` |

Supported flank sizes are `50, 100, 200, 300, 400, 500, 1000`

Multiple flank sizes can be supplied as a comma-separated list, such as `100,200,400`.

Use `all` to evaluate all supported flank sizes.

### Input Coordinates

The `coordinates` field describes the target region in the coordinate system of the **original input BAM**.

These coordinates are used for regional read extraction, regional read statistics, and calculation of the assembly region after adding the requested flanking sequence.

## Coverage QC Configuration

Files in `qc/` define the genomic intervals used to evaluate coverage for a locus.

A QC configuration can contain coordinates for multiple reference assemblies. Input groups are associated with the appropriate reference coordinate system so that coverage is evaluated against the same reference used to generate the original BAM.

For example, the FCGR2/3 configuration contains locus-specific intervals for datasets aligned to GRCh38 and T2T/CHM13.

Conceptually, the configuration is organized as:

```yaml
region: FCGR2_3

datasets:
  LRS:
    reference: GRCh38
  HGSVC:
    reference: T2T

references:
  GRCh38:
    FCGR2A:
      coordinates: [...]
    FCGR3A:
      coordinates: [...]

  T2T:
    FCGR2A:
      coordinates: [...]
    FCGR3A:
      coordinates: [...]
```

The workflow uses the `group` value from `samples.tsv` to select the appropriate coordinate system.

Some features may contain multiple intervals when more than one corresponding region is present in a reference assembly.

## Structural Haplotype Configuration

Files in `haplotypes/` optionally define known structural haplotypes for locus-specific classification.

If structural haplotyping is not required for a target region, set `NA` in the `haplotype_config` column of `samples.tsv`.

A haplotype configuration defines:

- the locus or region,
- downstream evaluation coordinates,
- assembly evaluation parameters, and
- structural events associated with known haplotypes.

The structural events describe expected insertions or deletions used to classify an assembled haplotype.

## Coordinate Handling

The workflow may accept input BAMs aligned to different reference assemblies.

The `coordinates` field in `samples.tsv` specifies the target region on the reference used by each **input BAM**. These coordinates are used to extract the reads needed for local assembly and to calculate input read statistics.

After assembly, all assembled contigs are aligned to the **workflow reference** supplied when the pipeline is launched. When a `haplotype_config` is provided, its `evaluation_coordinates` define the corresponding target region on this reference for:

- assembly contiguity evaluation,
- target-region sequence reconstruction,
- structural haplotype classification, and
- assembly-reference difference analysis.

When `haplotype_config` is `NA`, the workflow uses the `coordinates` value from `samples.tsv` for downstream evaluation.

This allows datasets aligned to different references, such as GRCh38 and T2T/CHM13, to be assembled independently and then evaluated against the same workflow reference.

## Adding a New Dataset

To add another dataset for an existing locus:

1. Add a row to `samples.tsv`.
2. Specify the appropriate input-BAM coordinates.
3. Assign a `group` label.
4. Set `data_type` to `ONT` or `PB`.
5. Select the desired flank size or sizes.
6. Reference the appropriate QC configuration.
7. Provide a structural haplotype configuration or set `haplotype_config` to `NA`.

If the new group uses a different reference coordinate system for its input BAM, add the corresponding intervals to the locus-specific QC YAML and associate the group with that reference.

## Running with a Configuration

The sample configuration is supplied when the workflow is launched:

```bash
./submit_local_asm_workflow.sh config/samples.tsv /path/to/reference/genome.fasta
```

Paths specified within `samples.tsv` should be accessible from the compute environment where the workflow is executed.
