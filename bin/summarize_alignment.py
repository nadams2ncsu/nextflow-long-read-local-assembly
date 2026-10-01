#!/usr/bin/env python3

import argparse
import pandas as pd
import pysam
import yaml


################################################################################
# Summarizes local assembly/alignment status
#
# Reports:
#   - contiguity status: contiguous, fragmented, discontiguous, unphased
#   - assembly region length
#   - SNVs
#   - small insertions/deletions
#   - large insertions/deletions
# 
# Inputs: samples.tsv, aligned contig BAM files
# Outputs: sample summary on contiguity, variants, and haplotypes


############################################
# Helper functions

# Read genomic coordinates
def parse_coordinates(coordinates):
    if coordinates is None:
        raise ValueError(
            "Evaluation coordinates are None. "
            "Check --coordinates and --haplotype-config."
        )

    chrom, coords = str(coordinates).split(":")
    start, end = map(
        int,
        coords.replace(",", "").split("-")
    )

    if start < 1 or end < start:
        raise ValueError(
            f"Invalid coordinates: {coordinates}"
        )

    return chrom, start, end

# Determine locus/gene coordinates used to evaluate the assembly
def get_evaluation_coordinates(args):

    # Default:
    # use coordinates supplied through samples.tsv
    evaluation_coordinates = args.coordinates

    # If a haplotype YAML is supplied, use its workflow-reference
    # evaluation coordinates instead.
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
            "No evaluation coordinates were available. "
            "Provide --coordinates or a haplotype YAML containing "
            "'evaluation_coordinates'."
        )

    return str(evaluation_coordinates)


# Calculate collective coverage if assembly is fragmented
def calculate_union_coverage(
    alignments,
    target_start,
    target_end
):

    intervals = []

    for aln, _, _ in alignments:

        overlap_start = max(
            aln.reference_start,
            target_start
        )

        overlap_end = min(
            aln.reference_end,
            target_end
        )

        if overlap_end > overlap_start:
            intervals.append(
                (overlap_start, overlap_end)
            )

    if not intervals:
        return 0

    intervals.sort()

    merged = []

    for start, end in intervals:

        if not merged or start > merged[-1][1]:
            merged.append(
                [start, end]
            )

        else:
            merged[-1][1] = max(
                merged[-1][1],
                end
            )

    return sum(
        end - start
        for start, end in merged
    )


# Calculate query-sequence length corresponding to target region
def calculate_assembly_region_length(
    aln,
    target_start,
    target_end
):

    assembly_length = 0
    ref_pos = aln.reference_start

    for operation, length in aln.cigartuples or []:

        # M, =, X
        if operation in (0, 7, 8):

            op_start = ref_pos
            op_end = ref_pos + length

            overlap_start = max(
                op_start,
                target_start
            )

            overlap_end = min(
                op_end,
                target_end
            )

            if overlap_end > overlap_start:
                assembly_length += (
                    overlap_end - overlap_start
                )

            ref_pos = op_end

        # I
        elif operation == 1:

            if target_start <= ref_pos < target_end:
                assembly_length += length

        # D, N
        elif operation in (2, 3):
            ref_pos += length

        # S, H, P
        elif operation in (4, 5, 6):
            continue

    return assembly_length


# Generate result for incomplete/unphased assembly
def empty_result(
    sample,
    group,
    gene,
    flank,
    coordinates,
    haplotype,
    status,
    region_length,
    covered_bases="NA",
    coverage_fraction="NA",
    contig_name="NA",
    assembly_region_length="NA"
):

    return {
        "sample": sample,
        "group": group,
        "gene": gene,
        "flank_kb": flank,
        "coordinates": coordinates,
        "assembly_haplotype": haplotype,
        "assembly_status": status,
        "contig_name": contig_name,
        "assembly_region_length": assembly_region_length,
        "region_length": region_length,
        "region_covered_bases": covered_bases,
        "region_coverage_fraction": coverage_fraction,
        "snv_count": "NA",
        "small_insertion_count": "NA",
        "small_deletion_count": "NA",
        "large_insertion_count": "NA",
        "large_insertion_positions": "NA",
        "large_deletion_count": "NA",
        "large_deletion_positions": "NA"
    }


# Write output
def write_result(
    result,
    output
):

    pd.DataFrame(
        [result]
    ).to_csv(
        output,
        sep="\t",
        index=False
    )


############################################
# Generate sample summary for local assembly

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize assembled-contig alignment across "
            "workflow-reference evaluation coordinates."
        )
    )

    parser.add_argument("--bam",required=True)
    parser.add_argument(
        "--coordinates",
        required=True,
        help=(
            "Fallback evaluation coordinates. "
            "When --haplotype-config is provided, "
            "evaluation_coordinates from the YAML are used instead."
        )
    )

    parser.add_argument(
        "--haplotype-config",
        default=None,
        help=(
            "Optional haplotype YAML containing evaluation_coordinates "
            "in the workflow-reference coordinate system."
        )
    )

    parser.add_argument("--gene",required=True)
    parser.add_argument("--sample",required=True)
    parser.add_argument("--group",required=True)
    parser.add_argument(
        "--flank",
        type=int,
        required=True
    )

    parser.add_argument("--assembly-haplotype",required=True)
    parser.add_argument(
        "--large-indel-threshold",
        type=int,
        default=40000
    )

    parser.add_argument(
        "--min-region-coverage",
        type=float,
        default=0.90
    )

    parser.add_argument("--output",required=True)
    args = parser.parse_args()


    ### Validate arguments
    if not 0 < args.min_region_coverage <= 1:
        raise ValueError(
            "--min-region-coverage must be > 0 and <= 1"
        )

    if args.large_indel_threshold < 1:
        raise ValueError(
            "--large-indel-threshold must be >= 1"
        )


    ### Determine evaluation coordinates
    evaluation_coordinates = get_evaluation_coordinates(
        args
    )

    chrom, region_start, region_end = parse_coordinates(
        evaluation_coordinates
    )

    # Input coordinates are 1-based inclusive.
    # pysam uses 0-based half-open coordinates.
    target_start = region_start - 1
    target_end = region_end

    region_length = (
        target_end - target_start
    )


    ### Read primary alignments overlapping target
    bam = pysam.AlignmentFile(
        args.bam,
        "rb"
    )

    if chrom not in bam.references:
        bam.close()
        raise ValueError(
            f"Chromosome {chrom} from {evaluation_coordinates} "
            "is not present in the assembled-contig BAM."
        )

    alignments = []

    for aln in bam.fetch(
        chrom,
        target_start,
        target_end
    ):

        if (
            aln.is_unmapped
            or aln.is_secondary
            or aln.is_supplementary
        ):
            continue

        overlap_start = max(
            aln.reference_start,
            target_start
        )

        overlap_end = min(
            aln.reference_end,
            target_end
        )

        covered_bases = max(
            0,
            overlap_end - overlap_start
        )

        if covered_bases == 0:
            continue

        coverage_fraction = (
            covered_bases /
            region_length
        )

        alignments.append(
            (
                aln,
                covered_bases,
                coverage_fraction
            )
        )


    ### No primary alignment
    if not alignments:
        bam.close()
        result = empty_result(
            args.sample,
            args.group,
            args.gene,
            args.flank,
            evaluation_coordinates,
            args.assembly_haplotype,
            "unphased",
            region_length
        )

        write_result(
            result,
            args.output
        )

        return


    ### Select best single contig
    aln, covered_bases, coverage_fraction = max(
        alignments,
        key=lambda x: x[1]
    )

    contig_name = aln.query_name

    ### Determine assembly status
    if coverage_fraction >= args.min_region_coverage:
        assembly_status = "contiguous"
        final_covered_bases = covered_bases
        final_coverage_fraction = coverage_fraction

    else:

        union_covered_bases = calculate_union_coverage(
            alignments,
            target_start,
            target_end
        )

        union_coverage_fraction = (
            union_covered_bases /
            region_length
        )

        distinct_contigs = {
            item[0].query_name
            for item in alignments
        }

        if (
            len(distinct_contigs) >= 2
            and union_coverage_fraction
            >= args.min_region_coverage
        ):

            assembly_status = "fragmented"

        else:

            assembly_status = "discontiguous"

        final_covered_bases = union_covered_bases
        final_coverage_fraction = union_coverage_fraction


    ### Incomplete assembly
    if assembly_status != "contiguous":
        bam.close()
        result = empty_result(
            args.sample,
            args.group,
            args.gene,
            args.flank,
            evaluation_coordinates,
            args.assembly_haplotype,
            assembly_status,
            region_length,
            final_covered_bases,
            round(
                final_coverage_fraction,
                4
            ),
            contig_name,
            "NA"
        )

        write_result(
            result,
            args.output
        )

        return

    ### Assembly sequence length across target
    assembly_region_length = calculate_assembly_region_length(
        aln,
        target_start,
        target_end
    )


    ### Sequence differences
    snvs = 0

    small_insertions = 0
    small_deletions = 0

    large_insertions = []
    large_deletions = []

    ref_pos = aln.reference_start

    for operation, length in aln.cigartuples or []:

        # M
        if operation == 0:
            ref_pos += length

        # I
        elif operation == 1:

            if target_start <= ref_pos < target_end:

                if length >= args.large_indel_threshold:

                    anchor = ref_pos + 1

                    large_insertions.append(
                        (
                            anchor,
                            anchor,
                            length
                        )
                    )

                else:

                    small_insertions += 1

        # D
        elif operation == 2:

            deletion_start = ref_pos
            deletion_end = ref_pos + length

            overlaps_target = (
                deletion_start < target_end
                and deletion_end > target_start
            )

            if overlaps_target:

                if length >= args.large_indel_threshold:

                    deletion_start_1based = (
                        deletion_start + 1
                    )

                    deletion_end_1based = (
                        deletion_end
                    )

                    large_deletions.append(
                        (
                            deletion_start_1based,
                            deletion_end_1based,
                            length
                        )
                    )

                else:

                    small_deletions += 1

            ref_pos += length

        # N
        elif operation == 3:
            ref_pos += length

        # =
        elif operation == 7:
            ref_pos += length

        # X
        elif operation == 8:

            mismatch_start = ref_pos
            mismatch_end = ref_pos + length

            overlap_start = max(
                mismatch_start,
                target_start
            )

            overlap_end = min(
                mismatch_end,
                target_end
            )

            if overlap_end > overlap_start:

                snvs += (
                    overlap_end - overlap_start
                )

            ref_pos += length

        # S, H, P
        elif operation in (4, 5, 6):
            continue


    bam.close()


    ### Sanity checks
    if not 0 <= coverage_fraction <= 1:
        raise RuntimeError(
            f"Invalid coverage fraction: {coverage_fraction}"
        )

    if covered_bases > region_length:
        raise RuntimeError(
            "Covered bases exceed evaluation-region length."
        )

    if assembly_region_length <= 0:
        raise RuntimeError(
            "Contiguous assembly has assembly_region_length <= 0."
        )


    ### Final result
    result = {
        "sample": args.sample,
        "group": args.group,
        "gene": args.gene,
        "flank_kb": args.flank,

        # These are the coordinates actually used to evaluate
        # the assembled-contig BAM.
        "coordinates": evaluation_coordinates,

        "assembly_haplotype": args.assembly_haplotype,
        "assembly_status": "contiguous",
        "contig_name": contig_name,
        "assembly_region_length": assembly_region_length,
        "region_length": region_length,
        "region_covered_bases": covered_bases,
        "region_coverage_fraction": round(
            coverage_fraction,
            4
        ),
        "snv_count": snvs,
        "small_insertion_count": small_insertions,
        "small_deletion_count": small_deletions,
        "large_insertion_count": len(
            large_insertions
        ),
        "large_insertion_positions": ";".join(
            f"{chrom}:{start}:{end}:{length}"
            for start, end, length
            in large_insertions
        ) or "NA",
        "large_deletion_count": len(
            large_deletions
        ),
        "large_deletion_positions": ";".join(
            f"{chrom}:{start}:{end}:{length}"
            for start, end, length
            in large_deletions
        ) or "NA"
    }

    write_result(
        result,
        args.output
    )


if __name__ == "__main__":
    main()
