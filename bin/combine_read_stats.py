#!/usr/bin/env python3

import argparse
import pandas as pd

########################################################################3
# Combine all sample read stats into a summary text file
# 
# - inputs: read stat qc files
#  -outputs: summary of all samples read stats


################################
#  combine files

def main():

    ## process print statement
    parser = argparse.ArgumentParser(
        description="Combine regional read-statistic TSV files."
    )

    ## required arguments
    parser.add_argument(
        "--input",
        nargs="+",
        required=True
    )

    parser.add_argument("--output",required=True)
    args = parser.parse_args()

    ## store results in summary df
    dataframes = [
        pd.read_csv(
            file,
            sep="\t"
        )
        for file in args.input
    ]

    if dataframes:
        combined = pd.concat(
            dataframes,
            ignore_index=True
        )

        combined = combined.sort_values(
            by=[
                "sample",
                "group",
                "gene"
            ]
        )

    else:
        combined = pd.DataFrame()

    combined.to_csv(
        args.output,
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    main()
