#!/usr/bin/env python3

import argparse
import pysam
import yaml

########################################################################
# Reconstruct haplotype sequences from aligned contigs
#   - Only occurs if both haplotype sequences are contiguous
#   - Trims haplotype sequences to gene region
#   - Preserves large SVs
#
# input: aligned contig BAM files
# output: haplotype-specific fasta files

#########################
# Helper Functions

## parse user-supplied coordinates
def parse_coordinates(coordinates):
    chrom, positions = coordinates.split(":")
    start, end = positions.split("-")
    return chrom, int(start) - 1, int(end)


# Get locus/gene coordinates used to reconstruct the assembly
def get_evaluation_coordinates(args):

    # Default to coordinates supplied in samples.tsv
    evaluation_coordinates = args.coordinates

    # If a haplotype configuration is supplied, use its
    # workflow-reference evaluation coordinates instead
    if args.haplotype_config:

        with open(args.haplotype_config) as handle:
            config = yaml.safe_load(handle)

        if not isinstance(config, dict):
            raise ValueError(
                "Haplotype configuration must be a YAML mapping."
            )

        evaluation_coordinates = config.get(
            "evaluation_coordinates"
        )

        if not evaluation_coordinates:
            raise ValueError(
                f"Haplotype config {args.haplotype_config} does not contain "
                "'evaluation_coordinates'. Add the coordinates of the complete "
                "target locus in the workflow-reference coordinate system."
            )

    if not evaluation_coordinates:
        raise ValueError(
            "No evaluation coordinates were available."
        )

    return str(evaluation_coordinates)


###################################################
# Read alignment CIGAR string for reconstruction
#   - possible with the "--eqx" minimap2 parameter

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct an assembled contig sequence corresponding "
            "to a reference-coordinate interval."
        )
    )

    ## required arguments
    parser.add_argument("--bam",required=True)
    parser.add_argument("--coordinates",required=True)
    parser.add_argument("--haplotype-config",default=None)
    parser.add_argument("--output",required=True)
    args = parser.parse_args()


    ### Determine evaluation coordinates
    evaluation_coordinates = get_evaluation_coordinates(
        args
    )

    chrom, target_start, target_end = parse_coordinates(
        evaluation_coordinates
    )


    ### Find primary alignment overlapping target region
    with pysam.AlignmentFile(args.bam, "rb") as bam:

        alignments = [
            aln
            for aln in bam.fetch(
                chrom,
                target_start,
                target_end
            )
            if not (
                aln.is_unmapped
                or aln.is_secondary
                or aln.is_supplementary
            )
        ]


    if not alignments:
        raise ValueError(
            f"No primary alignment overlaps "
            f"{evaluation_coordinates} in {args.bam}"
        )


    ###################################################
    # Use primary alignment with greatest overlap
    # with target interval

    aln = max(
        alignments,
        key=lambda x: (
            min(x.reference_end, target_end)
            - max(x.reference_start, target_start)
        )
    )


    if aln.query_sequence is None:
        raise ValueError(
            f"No query sequence stored for contig "
            f"{aln.query_name}"
        )


    ###################################################
    # Reconstruct query sequence corresponding to
    # target reference interval

    query_sequence = aln.query_sequence
    query_pos = 0
    ref_pos = aln.reference_start

    sequence_parts = []

    for operation, length in aln.cigartuples or []:

        # M, =, X
        if operation in (0, 7, 8):

            for i in range(length):
                current_ref = ref_pos + i

                if (
                    target_start
                    <= current_ref
                    < target_end
                ):

                    sequence_parts.append(
                        query_sequence[
                            query_pos + i
                        ]
                    )

            query_pos += length
            ref_pos += length


        # Insertion
        elif operation == 1:

            # Retain insertions occurring inside
            # the target interval
            if (
                target_start
                <= ref_pos
                < target_end
            ):

                sequence_parts.append(
                    query_sequence[
                        query_pos:
                        query_pos + length
                    ]
                )

            query_pos += length


        # Deletion / reference skip
        elif operation in (2, 3):
            ref_pos += length


        # Soft clipping
        elif operation == 4:
            query_pos += length

        # Hard clipping / padding
        elif operation in (5, 6):
            continue


    ###################################################
    # Join reconstructed sequence

    reconstructed_sequence = "".join(
        sequence_parts
    )

    if not reconstructed_sequence:
        raise ValueError(
            f"No sequence could be reconstructed for "
            f"{evaluation_coordinates} from {args.bam}"
        )


    ###################################################
    # Write FASTA

    with open(args.output, "w") as out:
        header = (
            f">{aln.query_name}"
            f"|region={evaluation_coordinates}"
        )

        out.write(
            header + "\n"
        )

        for i in range(
            0,
            len(reconstructed_sequence),
            80
        ):

            out.write(
                reconstructed_sequence[
                    i:i + 80
                ]
                + "\n"
            )


if __name__ == "__main__":
    main()
