#!/usr/bin/env python3

import argparse
import pandas as pd

################################################################################
# Combine local assembly results
# 
# Inputs: sample alignment summary tsv files
# Outputs:
#   1. Complete summary containing all assembly attempts.
#   2. Contiguous summary containing only sample/group/gene/flank combinations
#      where both hap1 and hap2 are contiguous.
#

###############################################################################
# Combine all sample-gene-flank specific alignment summary tsv files into 1

def main():
    parser = argparse.ArgumentParser(
        description="Combine per-assembly haplotype TSV files."
    )

    ## required arguments
    parser.add_argument(
        "--input",
        nargs="+",
        required=True
    )

    parser.add_argument("--output",required=True)
    parser.add_argument("--contiguous-output",required=True)
    args = parser.parse_args()

    ### Read and combine results
    dfs = [
        pd.read_csv(
            path,
            sep="\t",
            keep_default_na=False
        )
        for path in args.input
    ]

    combined = pd.concat(
        dfs,
        ignore_index=True
    )


    ### Normalize structural haplotype column
    if "structural_haplotype" not in combined.columns:
        combined["structural_haplotype"] = "NA"
    else:
        combined["structural_haplotype"] = (
            combined["structural_haplotype"]
            .fillna("NA")
            .replace("", "NA")
        )


    ### Write complete summary
    combined.to_csv(
        args.output,
        sep="\t",
        index=False
    )


    #########################
    # Identify paired contiguous assemblies

    group_columns = [
        "sample",
        "group",
        "gene",
        "flank_kb"
    ]

    contiguous = combined[
        combined["assembly_status"] == "contiguous"
    ]

    successful_groups = (
        contiguous
        .groupby(group_columns)["assembly_haplotype"]
        .agg(lambda x: set(x))
        .reset_index()
    )

    successful_groups = successful_groups[
        successful_groups["assembly_haplotype"].apply(
            lambda x: {"hap1", "hap2"}.issubset(x)
        )
    ][group_columns]


    #########################
    # Keep both haplotypes from successful groups

    contiguous_summary = combined.merge(
        successful_groups,
        on=group_columns,
        how="inner"
    )

    contiguous_summary = contiguous_summary[
        (contiguous_summary["assembly_status"] == "contiguous")
        & (
            contiguous_summary["assembly_haplotype"].isin(
                ["hap1", "hap2"]
            )
        )
    ]


    #########################
    # Write contiguous summary

    contiguous_summary.to_csv(
        args.contiguous_output,
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    main()
