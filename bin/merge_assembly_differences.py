#!/usr/bin/env python3

import argparse
import pandas as pd

###################################################################################
# Combine the variants observed in haplotype assemblies together for each sample
#
#	- this is done for assemblies with the same sampleID, group, gene, & flank size
#
# Inputs: differences (eg, variants) bed file for each sample
# Outputs: merges differences across haplotypes for the same sample into 1 file

################
# Input

## expected columns
COLUMNS = [
    "chrom",
    "start",
    "end",
    "difference_type",
    "haplotype",
    "notes",
    "reference_sequence",
    "assembly_sequence"
]

###########################
# Helper Functions

# 1. read in bed file with variants
def read_bed(path):
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        names=COLUMNS,
        dtype=str,
        keep_default_na=False
    )

    if df.empty:
        return df

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    return df

####################
# Merge assemblies 

def main():
    parser = argparse.ArgumentParser(
        description="Merge hap1 and hap2 assembly-reference differences."
    )

    parser.add_argument("--hap1", required=True)
    parser.add_argument("--hap2", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    hap1 = read_bed(args.hap1)
    hap2 = read_bed(args.hap2)

    key_columns = [
        "chrom",
        "start",
        "end",
        "difference_type",
        "reference_sequence",
        "assembly_sequence"
    ]

    position_columns = [
        "chrom",
        "start",
        "end",
        "difference_type"
    ]

    hap1_records = {
        tuple(row[col] for col in key_columns): row
        for _, row in hap1.iterrows()
    }

    hap2_records = {
        tuple(row[col] for col in key_columns): row
        for _, row in hap2.iterrows()
    }

    hap1_positions = {
        tuple(row[col] for col in position_columns)
        for _, row in hap1.iterrows()
    }

    hap2_positions = {
        tuple(row[col] for col in position_columns)
        for _, row in hap2.iterrows()
    }

    all_keys = set(hap1_records) | set(hap2_records)

    output_rows = []

    for key in sorted(
        all_keys,
        key=lambda x: (x[0], int(x[1]), int(x[2]), x[3], x[4], x[5])
    ):

        in_hap1 = key in hap1_records
        in_hap2 = key in hap2_records

        row = (
            hap1_records[key]
            if in_hap1
            else hap2_records[key]
        )

        position_key = tuple(
            row[col]
            for col in position_columns
        )

        if in_hap1 and in_hap2:

            haplotype = "hap1,hap2"
            notes = "Shared identical difference in hap1 and hap2"

        elif (
            position_key in hap1_positions
            and position_key in hap2_positions
        ):

            haplotype = "hap1" if in_hap1 else "hap2"
            notes = (
                "Same genomic position in hap1 and hap2 "
                "but different sequence"
            )

        else:

            haplotype = "hap1" if in_hap1 else "hap2"
            notes = "Haplotype-specific difference"

        output_rows.append({
            "chrom": row["chrom"],
            "start": row["start"],
            "end": row["end"],
            "difference_type": row["difference_type"],
            "haplotype": haplotype,
            "notes": notes,
            "reference_sequence": row["reference_sequence"],
            "assembly_sequence": row["assembly_sequence"]
        })

    output_df = pd.DataFrame(
        output_rows,
        columns=COLUMNS
    )

    with open(args.output, "w") as handle:

        handle.write(
            "#chrom\tstart\tend\tdifference_type\t"
            "haplotype\tnotes\treference_sequence\t"
            "assembly_sequence\n"
        )

        if not output_df.empty:

            output_df.to_csv(
                handle,
                sep="\t",
                index=False,
                header=False
            )


if __name__ == "__main__":
    main()
