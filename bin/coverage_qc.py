#!/usr/bin/env python3

import argparse
import subprocess
import pandas as pd
import yaml

#######################################333##############################################################
# Calculate coverage fold enrichment of target genes/region compared to genome average
#
#	- Calculates avergage coverage of gene & genome-wide average coverage
#	- Calculates coverage fold enrichment for each input coordinate
#	- Flags samples with genome-wide and/or genes/regions with < 20x coverage 
#
# Inputs: samples.tsv file && WGS BAM
# Output: Sample && region specific depth of coverage values

#################################
# Helper functions

## generic coverage calculation func
def run_coverage(bam, region=None):
    cmd = [
        "samtools",
        "coverage"
    ]

    if region:
        cmd.extend([
            "-r",
            region
        ])

    cmd.append(bam)

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True
    )

    rows = []

    for line in result.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        rows.append({
            "chrom": fields[0],
            "start": int(fields[1]),
            "end": int(fields[2]),
            "mean_depth": float(fields[6])
        })

    return rows

## Calculate whole-genome average coverage
def calculate_genome_depth(bam):
    rows = run_coverage(bam)

    if not rows:
        raise RuntimeError(
            "samtools coverage returned no whole-genome coverage data."
        )

    total_bases = 0
    weighted_depth = 0.0

    for row in rows:
        length = row["end"] - row["start"] + 1
        total_bases += length
        weighted_depth += (
            row["mean_depth"] * length
        )

    if total_bases == 0:
        raise RuntimeError(
            "Genome length calculated from samtools coverage was zero."
        )

    return weighted_depth / total_bases


## calculate gene/region coverage from coordinates in config/qc/YAML
def calculate_interval_depth(bam, region):
    rows = run_coverage(
        bam,
        region
    )

    if not rows:
        return 0.0

    total_bases = 0
    weighted_depth = 0.0

    for row in rows:
        length = row["end"] - row["start"] + 1
        total_bases += length
        weighted_depth += (
            row["mean_depth"] * length
        )

    if total_bases == 0:
        return 0.0
    return weighted_depth / total_bases


#####################################################################
# Calculate depth of coverage values for whole-genome & region

def main():

    ## print statement of proceess
    parser = argparse.ArgumentParser(
        description=(
            "Calculate whole-genome and gene-level "
            "coverage QC from an aligned BAM."
        )
    )

    ## required arguments
    parser.add_argument("--bam",required=True)
    parser.add_argument(
        "--config",
        required=True,
        help="QC YAML containing dataset/reference and gene coordinates."
    )

    parser.add_argument("--sample",required=True)
    parser.add_argument("--group",required=True)
    parser.add_argument(
        "--gene",
        required=True,
        help="Region/locus name, for example FCGR2_3."
    )

    parser.add_argument("--output",required=True)
    args = parser.parse_args()

    ### Read QC configuration
    with open(args.config) as handle:
        config = yaml.safe_load(handle)

    ### Determine BAM reference genome
    datasets = config.get(
        "datasets",
        {}
    )

    if args.group not in datasets:
        raise ValueError(
            f"No QC reference mapping configured for "
            f"group '{args.group}'."
        )

    reference_build = datasets[
        args.group
    ]["reference"]

    references = config.get(
        "references",
        {}
    )

    if reference_build not in references:
        raise ValueError(
            f"Reference '{reference_build}' is not "
            f"defined in the QC configuration."
        )

    gene_config = references[
        reference_build
    ]

    ### Whole-genome coverage calculation
    genome_depth = calculate_genome_depth(
        args.bam
    )

    ### if depth is <= 20 at whole genome level issue warning
    if genome_depth <= 20:
        genome_warning = (
            "WARNING: whole-genome coverage <=20x; "
            "input coverage may be too low for reliable assembly."
        )

    else:
        genome_warning = "PASS"


    ### Gene-level coverage at user-supplied region from samples.tsv
    output_rows = []

    for gene_name, settings in gene_config.items():
        coordinates = settings.get(
            "coordinates",
            []
        )

        if not coordinates:
            continue

        interval_depths = []

        for region in coordinates:
            depth = calculate_interval_depth(
                args.bam,
                region
            )

            interval_depths.append(
                depth
            )

        # Multiple configured intervals are deliberately
        # summed before normalization.
        # intended for extended haplotypes with the same gene name
        summed_gene_depth = sum(
            interval_depths
        )

        if genome_depth > 0:
            normalized_depth = (
                summed_gene_depth /
                genome_depth
            )

        else:
            normalized_depth = float("nan")

        #### Gene coverage warning; same as the whole-genome level
        low_intervals = [
            coordinates[i]
            for i, depth in enumerate(interval_depths)
            if depth <= 20
        ]

        if low_intervals:
            gene_warning = (
                "WARNING: gene interval coverage <=20x; "
                "low coverage may reflect a deletion or "
                "insufficient input data for reliable assembly."
            )

        else:
            gene_warning = "PASS"


        #### Store output
        output_rows.append({

            "sample":
                args.sample,

            "group":
                args.group,

            "region":
                args.gene,

            "reference_build":
                reference_build,

            "genome_depth":
                round(genome_depth, 3),

            "genome_qc":
                genome_warning,

            "gene":
                gene_name,

            "interval_count":
                len(coordinates),

            "coordinates":
                ";".join(coordinates),

            "interval_depths":
                ";".join(
                    f"{depth:.3f}"
                    for depth in interval_depths
                ),

            "summed_gene_depth":
                round(summed_gene_depth, 3),

            "normalized_depth":
                round(normalized_depth, 3),

            "gene_qc":
                gene_warning
        })


    #################################
    # Write output

    df = pd.DataFrame(
        output_rows
    )

    df.to_csv(
        args.output,
        sep="\t",
        index=False
    )


if __name__ == "__main__":

    main()
