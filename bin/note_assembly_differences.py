#!/usr/bin/env python3

import argparse
import pysam
import yaml

###############################################################
# Note differences between assemblies and reference genome
#   - Reads CIGAR to note assembly-reference differences
#     including SNVs, INDELs, and SVs
#
# Input: aligned contig BAM files
# Output: BED files for each haplotype assembly
#

###########################
# Helper Functions

# 1. Parse BAM file coordinates
def parse_coordinates(value):

    chrom, positions = value.split(":")
    start, end = positions.replace(",", "").split("-")

    # Input coordinates: 1-based inclusive.
    # Internal coordinates: 0-based half-open.
    return chrom, int(start) - 1, int(end)


# 2. Get coordinates used to evaluate assembly-reference differences
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


# 3. Get primary alignment for contig
def get_primary_alignment(
    bam_path,
    chrom,
    start,
    end
):

    with pysam.AlignmentFile(bam_path, "rb") as bam:

        alignments = [
            aln
            for aln in bam.fetch(
                chrom,
                start,
                end
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
            f"{chrom}:{start + 1}-{end}"
        )

    return max(
        alignments,
        key=lambda aln: (
            min(aln.reference_end, end)
            - max(aln.reference_start, start)
        )
    )


#########################################################################################
# Parse BAM file CIGAR to note differences between assemblies and reference genome

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Record sequence differences between a haplotype assembly "
            "and the reference genome."
        )
    )

    ### required arguments
    parser.add_argument("--bam",required=True)
    parser.add_argument("--reference",required=True)
    parser.add_argument("--haplotype",required=True)
    parser.add_argument("--coordinates",required=True)
    parser.add_argument("--haplotype-config",default=None)
    parser.add_argument(
        "--large-indel-threshold",
        type=int,
        default=40000
    )

    parser.add_argument("--output",required=True)
    args = parser.parse_args()


    ### Determine workflow-reference evaluation coordinates
    evaluation_coordinates = get_evaluation_coordinates(
        args
    )

    chrom, target_start, target_end = parse_coordinates(
        evaluation_coordinates
    )


    ### Open reference
    reference = pysam.FastaFile(
        args.reference
    )


    ### Get primary alignment overlapping evaluation region
    aln = get_primary_alignment(
        args.bam,
        chrom,
        target_start,
        target_end
    )

    query_sequence = aln.query_sequence

    if query_sequence is None:
        raise ValueError(
            f"No query sequence stored for {aln.query_name}"
        )


    ### Identify assembly-reference differences
    rows = []
    query_pos = 0
    ref_pos = aln.reference_start

    for operation, length in aln.cigartuples or []:

        ##########################
        # Match (=)

        if operation == 7:

            query_pos += length
            ref_pos += length


        ##########################
        # Mismatch (X)

        elif operation == 8:

            for i in range(length):

                position = ref_pos + i

                if not (
                    target_start
                    <= position
                    < target_end
                ):
                    continue

                reference_base = reference.fetch(
                    chrom,
                    position,
                    position + 1
                ).upper()

                assembly_base = query_sequence[
                    query_pos + i
                ].upper()

                rows.append([
                    chrom,
                    position,
                    position + 1,
                    "SNV",
                    args.haplotype,
                    "Haplotype-specific assembly-reference difference",
                    reference_base,
                    assembly_base
                ])

            query_pos += length
            ref_pos += length


        ##########################
        # Generic M

        elif operation == 0:

            query_pos += length
            ref_pos += length


        ##################################
        # Insertion relative to reference

        elif operation == 1:

            inserted_sequence = query_sequence[
                query_pos:
                query_pos + length
            ].upper()

            if (
                target_start
                <= ref_pos
                < target_end
            ):

                difference_type = (
                    "large_insertion"
                    if length
                    >= args.large_indel_threshold
                    else "insertion"
                )

                rows.append([
                    chrom,
                    ref_pos,
                    ref_pos + 1,
                    difference_type,
                    args.haplotype,
                    "Haplotype-specific assembly-reference difference",
                    "-",
                    inserted_sequence
                ])

            query_pos += length


        ##################################
        # Deletion relative to reference

        elif operation == 2:

            deletion_start = ref_pos
            deletion_end = ref_pos + length

            # Report the complete deletion if it overlaps
            # the target evaluation region.
            if (
                deletion_start < target_end
                and deletion_end > target_start
            ):

                reference_sequence = reference.fetch(
                    chrom,
                    deletion_start,
                    deletion_end
                ).upper()

                difference_type = (
                    "large_deletion"
                    if length
                    >= args.large_indel_threshold
                    else "deletion"
                )

                rows.append([
                    chrom,
                    deletion_start,
                    deletion_end,
                    difference_type,
                    args.haplotype,
                    "Haplotype-specific assembly-reference difference",
                    reference_sequence,
                    "-"
                ])

            ref_pos += length


        ##########################
        # Reference skip

        elif operation == 3:
            ref_pos += length


        ##########################
        # Soft clipping

        elif operation == 4:
            query_pos += length


        ##########################
        # Hard clipping / padding

        elif operation in (5, 6):
            continue

    reference.close()


    ###############################################################
    # Write BED-like output

    with open(args.output, "w") as out:

        out.write(
            "#chrom\tstart\tend\tdifference_type\t"
            "haplotype\tnotes\treference_sequence\t"
            "assembly_sequence\n"
        )

        for row in rows:

            out.write(
                "\t".join(
                    map(str, row)
                )
                + "\n"
            )


if __name__ == "__main__":
    main()
