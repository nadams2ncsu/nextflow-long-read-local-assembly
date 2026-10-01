#!/usr/bin/env python3

import argparse
import statistics
import pysam
import pandas as pd

################################################################################
# Calculate regional read-length statistics from an aligned BAM.
#
# 	- All reads that overlap the region are counted towards the statistic
#
# - inputs: samples.tsv && WGS BAM
# - outputs: read number of lengths (min, median, max) summary table


###########################
# Helper functions

def parse_coordinates(coordinates):
    chrom, coords = coordinates.split(":")
    start, end = map(
        int,
        coords.replace(",", "").split("-")
    )

    # Convert to 0-based half-open coordinates for pysam
    return chrom, start - 1, end

##################################################
# Parse user-supplied BAM file coordinate region
#	for read statistics

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate read-length statistics for reads overlapping "
            "a user-supplied genomic interval."
        )
    )

    ## required arguments (eg, supplied in samples.tsv)
    parser.add_argument("--bam",required=True)
    parser.add_argument("--coordinates",required=True)
    parser.add_argument("--sample",required=True)
    parser.add_argument("--group",required=True)
    parser.add_argument("--gene",required=True)

    ## minimum length, issues warning if below but analysis continues
    parser.add_argument(
        "--min-median-read-length",
        type=int,
        default=5000,
        help=(
            "Minimum median read length before a WARNING is reported. "
            "Default: 5000 bp."
        )
    )

    parser.add_argument("--output",required=True)
    args = parser.parse_args()

    ## Parse target coordinates
    chrom, start, end = parse_coordinates(
        args.coordinates
    )


    ################################
    # Collect read lengths 

    read_lengths = []

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(
            chrom,
            start,
            end
        ):

            # Ignore non-primary alignments
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
            ):
                continue

            read_length = read.query_length

            if read_length is not None:
                read_lengths.append(read_length)


    ####################################
    # Calculate statistics

    if read_lengths:

        read_count = len(read_lengths)
        min_length = min(read_lengths)
        median_length = statistics.median(read_lengths)
        max_length = max(read_lengths)

        qc = (
            "WARNING"
            if median_length <= args.min_median_read_length
            else "PASS"
        )

    else:

        read_count = 0
        min_length = "NA"
        median_length = "NA"
        max_length = "NA"
        qc = "WARNING"

    ############
    # Output

    result = {
        "sample": args.sample,
        "group": args.group,
        "gene": args.gene,
        "coordinates": args.coordinates,
        "read_count": read_count,
        "min_read_length": min_length,
        "median_read_length": median_length,
        "max_read_length": max_length,
        "read_length_qc": qc
    }

    pd.DataFrame(
        [result]
    ).to_csv(
        args.output,
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    main()
