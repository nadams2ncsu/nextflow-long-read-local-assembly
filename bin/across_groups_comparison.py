#!/usr/bin/env python3

import argparse
import itertools
import pandas as pd
import pysam

###################################################################################################################
#
# Compare paired-contiguous haplotype assemblies across groups.
#
# 	- Samples represented in two or more groups for the same gene are compared.
#	- All valid paired-contiguous flank combinations are compared across groups.
# 		The assignment with the highest mean sequence similarity is retained.
#
# Inputs: Contig alignment files
# Outputs: Summary of sequence similarity of haplotype assemblies across groups for the same individual
#


################################
# Helper functions

## load contig alignment file to reference genome
def load_contig(bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        alignments = [
            aln
            for aln in bam.fetch(until_eof=True)
            if not (
                aln.is_unmapped
                or aln.is_secondary
                or aln.is_supplementary
            )
        ]

    if not alignments:
        raise ValueError(
            f"No primary alignment found in {bam_path}"
        )

    return max(
        alignments,
        key=lambda x: x.query_alignment_length
    )


## only compare sequence in user-supplied coordinate region in samples.tsv
def reference_anchored_sequence(aln):
    query_sequence = aln.query_sequence

    if query_sequence is None:
        raise ValueError(
            f"No query sequence stored for contig {aln.query_name}"
        )

    positions = {}
    query_pos = 0
    ref_pos = aln.reference_start

    for operation, length in aln.cigartuples or []:
        # M, =, X
        if operation in (0, 7, 8):
            for i in range(length):
                positions[ref_pos + i] = (
                    query_sequence[query_pos + i]
                )

            query_pos += length
            ref_pos += length

        # Insertion relative to reference
        elif operation == 1:
            query_pos += length

        # Deletion / skipped reference region
        elif operation in (2, 3):

            for i in range(length):
                positions[ref_pos + i] = "-"

            ref_pos += length

        # Soft clipping
        elif operation == 4:

            query_pos += length

        # Hard clipping / padding
        elif operation in (5, 6):
            continue
    return positions

## calculate per base similarity
def calculate_similarity(contig1, contig2):
    seq1 = contig1["sequence"]
    seq2 = contig2["sequence"]
    shared_positions = (
        set(seq1.keys())
        & set(seq2.keys())
    )
    if not shared_positions:
        return None

    matches = sum(
        seq1[pos].upper() == seq2[pos].upper()
        for pos in shared_positions
    )

    return (
        matches / len(shared_positions)
    ) * 100.0


def assignment_mean(similarity1, similarity2):
    # Both haplotypes must have a valid comparison.
    if similarity1 is None or similarity2 is None:
        return -1.0

    return (
        similarity1 + similarity2
    ) / 2.0

## read in which samples to compare
## sample && region will be the same; flank, data_type will vary
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
        "assembly_haplotype",
        "bam"
    }

    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            "Manifest missing required columns: "
            + ", ".join(sorted(missing))
        )

    return df

## read in aligned contg
def load_assembly(row):
    aln = load_contig(
        row["bam"]
    )

    return {
        "sample": str(row["sample"]),
        "group": str(row["group"]),
        "gene": str(row["gene"]),
        "flank_kb": str(row["flank_kb"]),
        "haplotype": str(row["assembly_haplotype"]),
        "contig_id": aln.query_name,
        "contig_length": aln.query_length,
        "sequence": reference_anchored_sequence(aln)
    }

#################################################################
# compare assemblies across groups

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare paired-contiguous haplotype assemblies "
            "across groups for the same sample and gene."
        )
    )

    ## required argumemts
    parser.add_argument("--manifest",required=True)
    parser.add_argument(
        "--similarity-threshold",
        type=float,
        default=95.0
    )

    parser.add_argument("--output",required=True)
    args = parser.parse_args()
    manifest = read_manifest(
        args.manifest
    )

    ### Load assemblies
    assemblies = {}

    for _, row in manifest.iterrows():
        key = (
            str(row["sample"]),
            str(row["group"]),
            str(row["gene"]),
            str(row["flank_kb"]),
            str(row["assembly_haplotype"])
        )

        if key in assemblies:
            raise ValueError(
                f"Duplicate assembly in manifest: {key}"
            )

        assemblies[key] = load_assembly(
            row
        )

    results = []

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

        # Samples present in only one group are skipped.
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

            ### Compare every valid cross-group flank combination
            for flank1, flank2 in itertools.product(
                group1_flanks,
                group2_flanks
            ):

                key_1_hap1 = (
                    sample,
                    group1,
                    gene,
                    flank1,
                    "hap1"
                )

                key_1_hap2 = (
                    sample,
                    group1,
                    gene,
                    flank1,
                    "hap2"
                )

                key_2_hap1 = (
                    sample,
                    group2,
                    gene,
                    flank2,
                    "hap1"
                )

                key_2_hap2 = (
                    sample,
                    group2,
                    gene,
                    flank2,
                    "hap2"
                )

                required_keys = [
                    key_1_hap1,
                    key_1_hap2,
                    key_2_hap1,
                    key_2_hap2
                ]

                if not all(
                    key in assemblies
                    for key in required_keys
                ):
                    continue

                g1_h1 = assemblies[key_1_hap1]
                g1_h2 = assemblies[key_1_hap2]

                g2_h1 = assemblies[key_2_hap1]
                g2_h2 = assemblies[key_2_hap2]

                ### Same assignment (hifiasm labels were the same)
                same_1 = calculate_similarity(
                    g1_h1,
                    g2_h1
                )

                same_2 = calculate_similarity(
                    g1_h2,
                    g2_h2
                )

                same_mean = assignment_mean(
                    same_1,
                    same_2
                )

                ### Swapped assignment (hifiasm labels were swapped)
                swap_1 = calculate_similarity(
                    g1_h1,
                    g2_h2
                )

                swap_2 = calculate_similarity(
                    g1_h2,
                    g2_h1
                )

                swap_mean = assignment_mean(
                    swap_1,
                    swap_2
                )

                ### Retain best haplotype assignment
                if same_mean >= swap_mean:
                    assignment = "same"
                    comparisons = [
                        (
                            g1_h1,
                            g2_h1,
                            same_1
                        ),
                        (
                            g1_h2,
                            g2_h2,
                            same_2
                        )
                    ]

                else:

                    assignment = "swapped"

                    comparisons = [
                        (
                            g1_h1,
                            g2_h2,
                            swap_1
                        ),
                        (
                            g1_h2,
                            g2_h1,
                            swap_2
                        )
                    ]

                ### Store selected comparisons
                for contig1, contig2, similarity in comparisons:
                    if similarity is None:
                        similarity_output = None
                        qc = "NO_SHARED_ALIGNMENT"

                    else:

                        similarity_output = round(
                            similarity,
                            4
                        )

                        qc = (
                            "WARNING"
                            if similarity
                            < args.similarity_threshold
                            else "PASS"
                        )

                    results.append({

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

                        "haplotype_assignment":
                            assignment,

                        "group1_assembly_haplotype":
                            contig1["haplotype"],

                        "group2_assembly_haplotype":
                            contig2["haplotype"],

                        "contig_id_1":
                            contig1["contig_id"],

                        "contig_id_2":
                            contig2["contig_id"],

                        "contig_length_1":
                            contig1["contig_length"],

                        "contig_length_2":
                            contig2["contig_length"],

                        "sequence_similarity":
                            similarity_output,

                        "comparison_qc":
                            qc

                    })
    ###########
    # Output

    columns = [
        "sample",
        "gene",
        "group1",
        "group1_flank_kb",
        "group2",
        "group2_flank_kb",
        "haplotype_assignment",
        "group1_assembly_haplotype",
        "group2_assembly_haplotype",
        "contig_id_1",
        "contig_id_2",
        "contig_length_1",
        "contig_length_2",
        "sequence_similarity",
        "comparison_qc",
        "sample_gene_all_comparisons_pass"
    ]

    if results:
        result_df = pd.DataFrame(
            results
        )

        #### Overall comparison QC for each sample + gene
        flags = (
            result_df
            .groupby(
                ["sample", "gene"],
                sort=False
            )["comparison_qc"]
            .apply(
                lambda x: bool(
                    (x == "PASS").all()
                )
            )
            .rename(
                "sample_gene_all_comparisons_pass"
            )
            .reset_index()
        )

        result_df = result_df.merge(
            flags,
            on=["sample", "gene"],
            how="left"
        )

        result_df = result_df[
            columns
        ]

    else:

        result_df = pd.DataFrame(
            columns=columns
        )

    result_df.to_csv(
        args.output,
        sep="\t",
        index=False,
        na_rep="NA"
    )


if __name__ == "__main__":
    main()
