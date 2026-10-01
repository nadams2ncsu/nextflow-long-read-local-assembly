#!/usr/bin/env python3

import argparse
import pandas as pd

#############################################################################################################
# Determines which samples have contiguous haplotype assemblies
#
#	- those with contiguous assemblies undergo:
#		i) trimming to the user-supplied gene coordinate
#		ii) comparison across groups (if applicable; eg, same sample sequenced twice)
#
#	- non-contiguous assemblies stop here (eg, saved in work/ but not copied over the final results/ directory
#
# Input: Per-haplotype summary TSV files
# Output: Manifest of paired-contiguous assemblies
#


#####################################3
# Contiguous assembly selector

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Identify sample/group/gene/flank assemblies where "
            "both hap1 and hap2 are contiguous."
        )
    )

    parser.add_argument(
        "--input",
        nargs="+",
        required=True,
        help="Per-haplotype SUMMARIZE TSV files"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output paired-contiguous manifest"
    )

    args = parser.parse_args()

    dfs = [
        pd.read_csv(path, sep="\t")
        for path in args.input
    ]

    df = pd.concat(
        dfs,
        ignore_index=True
    )

    required = {
        "sample",
        "group",
        "gene",
        "flank_kb",
        "assembly_haplotype",
        "assembly_status"
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(sorted(missing))
        )

    key_columns = [
        "sample",
        "group",
        "gene",
        "flank_kb"
    ]

    valid_rows = []

    for key, group_df in df.groupby(
        key_columns,
        dropna=False
    ):

        contiguous = group_df[
            group_df["assembly_status"] == "contiguous"
        ]

        haplotypes = set(
            contiguous["assembly_haplotype"]
        )

        # Require exactly one hap1 and exactly one hap2.
        hap1_count = (
            contiguous["assembly_haplotype"]
            == "hap1"
        ).sum()

        hap2_count = (
            contiguous["assembly_haplotype"]
            == "hap2"
        ).sum()

        if (
            haplotypes == {"hap1", "hap2"}
            and hap1_count == 1
            and hap2_count == 1
        ):

            sample, group, gene, flank_kb = key

            valid_rows.append({
                "sample": sample,
                "group": group,
                "gene": gene,
                "flank_kb": flank_kb
            })

    output_df = pd.DataFrame(
        valid_rows,
        columns=key_columns
    )

    output_df.to_csv(
        args.output,
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    main()
