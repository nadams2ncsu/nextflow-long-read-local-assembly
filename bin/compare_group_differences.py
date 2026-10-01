#!/usr/bin/env python3

import argparse
import itertools
import pandas as pd

################################################################################
#
# Compare merged assembly-reference difference BED files across groups.
#
# For each sample + gene:
#
#   - compare every pair of groups
#   - compare every available flank combination between those groups
#   - identify shared and discordant assembly-reference differences


####################
# I/O

## expected input structure
DIFFERENCE_COLUMNS = [
    "chrom",
    "start",
    "end",
    "difference_type",
    "haplotype",
    "notes",
    "reference_sequence",
    "assembly_sequence"
]

####################################33
# Helper functions

## read in sample combinations
def read_manifest(path):
    df = pd.read_csv(
        path,
        sep="\t",
        keep_default_na=False
    )

    required = {
        "sample",
        "group",
        "gene",
        "flank_kb",
        "bed"
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            "Manifest missing required columns: "
            + ", ".join(sorted(missing))
        )

    return df

## read in sample differences (eg, variants) between assembly && reference genome
def read_differences(path):
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        names=DIFFERENCE_COLUMNS,
        dtype=str,
        keep_default_na=False
    )

    return df

## setting up which samples to compare
## same sample require, but group & flank can vary
def difference_key(row):
    return (
        str(row["chrom"]),
        int(row["start"]),
        int(row["end"]),
        str(row["difference_type"]),
        str(row["reference_sequence"]),
        str(row["assembly_sequence"])
    )


## reading in differences
def load_difference_set(path):
    df = read_differences(path)

    differences = {}

    for _, row in df.iterrows():
        key = difference_key(row)

        if key in differences:
            raise ValueError(
                f"Duplicate difference within merged BED {path}: {key}"
            )

        differences[key] = {
            "chrom": str(row["chrom"]),
            "start": int(row["start"]),
            "end": int(row["end"]),
            "difference_type": str(row["difference_type"]),
            "haplotype": str(row["haplotype"]),
            "reference_sequence": str(
                row["reference_sequence"]
            ),
            "assembly_sequence": str(
                row["assembly_sequence"]
            )
        }

    return differences

###################################################################3########
# Compare sample differences with reference genome across groups
# (eg, compare variants of the same individual sequenced twice)

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare merged assembly-reference differences "
            "across groups."
        )
    )

    parser.add_argument("--manifest",required=True)
    parser.add_argument("--summary-output",required=True)
    parser.add_argument("--discordance-output",required=True)

    args = parser.parse_args()

    manifest = read_manifest(args.manifest)

    ### Load each merged BED
    assemblies = {}

    for _, row in manifest.iterrows():
        key = (
            str(row["sample"]),
            str(row["group"]),
            str(row["gene"]),
            str(row["flank_kb"])
        )

        if key in assemblies:
            raise ValueError(
                f"Duplicate merged difference BED: {key}"
            )

        assemblies[key] = load_difference_set(
            row["bed"]
        )

    summary_results = []
    discordant_results = []

    ### Compare independently for each sample + gene
    sample_genes = (
        manifest[
            ["sample", "gene"]
        ]
        .astype(str)
        .drop_duplicates()
    )

    for _, sample_gene in sample_genes.iterrows():
        sample = sample_gene["sample"]
        gene = sample_gene["gene"]
        sample_rows = manifest[
            (
                manifest["sample"].astype(str)
                == sample
            )
            &
            (
                manifest["gene"].astype(str)
                == gene
            )
        ]

        groups = sorted(
            sample_rows["group"]
            .astype(str)
            .unique()
        )

        # Nothing to compare across groups.
        if len(groups) < 2:
            continue

        ### Compare every pair of groups
        for group1, group2 in itertools.combinations(
            groups,
            2
        ):

            group1_rows = sample_rows[
                sample_rows["group"].astype(str)
                == group1
            ]

            group2_rows = sample_rows[
                sample_rows["group"].astype(str)
                == group2
            ]

            group1_flanks = sorted(
                group1_rows["flank_kb"]
                .astype(str)
                .unique(),
                key=lambda x: int(x)
            )

            group2_flanks = sorted(
                group2_rows["flank_kb"]
                .astype(str)
                .unique(),
                key=lambda x: int(x)
            )

            ### Compare every cross-group flank combination
            for flank1, flank2 in itertools.product(
                group1_flanks,
                group2_flanks
            ):

                key1 = (
                    sample,
                    group1,
                    gene,
                    flank1
                )

                key2 = (
                    sample,
                    group2,
                    gene,
                    flank2
                )

                if (
                    key1 not in assemblies
                    or key2 not in assemblies
                ):
                    continue

                differences1 = assemblies[key1]
                differences2 = assemblies[key2]

                set1 = set(
                    differences1.keys()
                )

                set2 = set(
                    differences2.keys()
                )

                shared = set1 & set2
                group1_only = set1 - set2
                group2_only = set2 - set1
                union = set1 | set2

                ### Jaccard concordance
                if union:
                    concordance = (
                        len(shared)
                        / len(union)
                    ) * 100.0

                else:

                    # Both assemblies have zero differences from reference.
                    concordance = 100.0

                summary_results.append({
                    "sample":
                        sample,
                    "gene":
                        gene,
                    "group1":
                        group1,
                    "group1_flank_kb":
                        flank1,
                    "group2":
                        group2,
                    "group2_flank_kb":
                        flank2,
                    "group1_difference_count":
                        len(set1),
                    "group2_difference_count":
                        len(set2),
                    "shared_difference_count":
                        len(shared),
                    "group1_only_count":
                        len(group1_only),
                    "group2_only_count":
                        len(group2_only),
                    "union_difference_count":
                        len(union),
                    "difference_concordance":
                        round(concordance, 4)
                })

                ### Store group1-only differences
                for key in sorted(group1_only):
                    difference = differences1[key]
                    discordant_results.append({
                        "chrom":
                            difference["chrom"],

                        "start":
                            difference["start"],
                        "end":
                            difference["end"],
                        "difference_type":
                            difference["difference_type"],
                        "haplotype":
                            difference["haplotype"],
                        "notes":
                            f"Present only in {group1} comparison",
                        "reference_sequence":
                            difference["reference_sequence"],
                        "assembly_sequence":
                            difference["assembly_sequence"],
                        "sample":
                            sample,
                        "gene":
                            gene,
                        "group1":
                            group1,
                        "group1_flank_kb":
                            flank1,
                        "group2":
                            group2,
                        "group2_flank_kb":
                            flank2,
                        "discordance":
                            "group1_only"
                    })

                ### Store group2-only differences
                for key in sorted(group2_only):
                    difference = differences2[key]
                    discordant_results.append({
                        "chrom":
                            difference["chrom"],
                        "start":
                            difference["start"],
                        "end":
                            difference["end"],
                        "difference_type":
                            difference["difference_type"],
                        "haplotype":
                            difference["haplotype"],
                        "notes":
                            f"Present only in {group2} comparison",
                        "reference_sequence":
                            difference["reference_sequence"],
                        "assembly_sequence":
                            difference["assembly_sequence"],
                        "sample":
                            sample,
                        "gene":
                            gene,
                        "group1":
                            group1,
                        "group1_flank_kb":
                            flank1,
                        "group2":
                            group2,
                        "group2_flank_kb":
                            flank2,
                        "discordance":
                            "group2_only"
                    })

    ######################
    ### Summary output

    summary_columns = [
        "sample",
        "gene",
        "group1",
        "group1_flank_kb",
        "group2",
        "group2_flank_kb",
        "group1_difference_count",
        "group2_difference_count",
        "shared_difference_count",
        "group1_only_count",
        "group2_only_count",
        "union_difference_count",
        "difference_concordance"
    ]

    summary_df = pd.DataFrame(
        summary_results,
        columns=summary_columns
    )

    summary_df.to_csv(
        args.summary_output,
        sep="\t",
        index=False
    )

    ### Detailed discordance BED
    discordance_columns = [
        "chrom",
        "start",
        "end",
        "difference_type",
        "haplotype",
        "notes",
        "reference_sequence",
        "assembly_sequence",
        "sample",
        "gene",
        "group1",
        "group1_flank_kb",
        "group2",
        "group2_flank_kb",
        "discordance"
    ]

    discordance_df = pd.DataFrame(
        discordant_results,
        columns=discordance_columns
    )

    discordance_df.to_csv(
        args.discordance_output,
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    main()
