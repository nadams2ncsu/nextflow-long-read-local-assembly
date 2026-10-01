#!/usr/bin/env python3

import argparse
import pandas as pd
import yaml


################################################################################
# Classify structural haplotypes from large structural events identified by
# summarize_alignment.py.
#
# Observed structural-event format:
#
#     chromosome:start:end:length
#
# Known structural events are matched using:
#
#     - event type
#     - chromosome
#     - approximate event length
#     - overlap with the configured structural-event region
#
# When more than one configured haplotype is compatible, the haplotype with
# the strongest interval overlap is selected.

###########################################
# Helper Functions

## parse user supplied coordinate regions
def parse_region(region):
    chrom, coords = region.split(":")
    start, end = map(
        int,
        coords.replace(",", "").split("-")
    )
    return chrom, start, end

## parse known haplotype events
def parse_events(value):
    if (
        pd.isna(value)
        or str(value).strip() in ["", "NA"]
    ):
        return []

    events = []

    for event in str(value).split(";"):
        fields = event.split(":")
        if len(fields) != 4:
            raise ValueError(
                f"Invalid structural-event format: {event}. "
                "Expected chrom:start:end:length"
            )

        chrom, start, end, length = fields
        events.append(
            (
                chrom,
                int(start),
                int(end),
                int(length)
            )
        )

    return events

## haplotype overlap && aligned contig overlap
def interval_overlap(
    start1,
    end1,
    start2,
    end2
):

    return max(
        0,
        min(end1, end2)
        - max(start1, start2)
        + 1
    )

## classifying SV haplotypes with some flexibility
def score_event(
    observed_event,
    expected_event,
    length_tolerance,
    breakpoint_tolerance
):
    """
    Return a matching score or None if the event is incompatible.

    Deletions:
        Require compatible chromosome and approximate length.
        Rank using interval overlap with the configured deletion region.

    Insertions:
        Require compatible chromosome and approximate length.
        Rank using distance of insertion anchor from configured interval.
    """

    (
        observed_chrom,
        observed_start,
        observed_end,
        observed_length
    ) = observed_event

    (
        expected_chrom,
        expected_start,
        expected_end
    ) = parse_region(
        expected_event["region"]
    )

    if observed_chrom != expected_chrom:
        return None

    expected_length = expected_event.get(
        "length"
    )

    if expected_length is not None:

        length_difference = abs(
            observed_length
            - int(expected_length)
        )

        if length_difference > length_tolerance:
            return None

    else:
        length_difference = 0

    event_type = expected_event["type"]

    #########
    ### Deletion
    if event_type == "deletion":

        overlap = interval_overlap(
            observed_start,
            observed_end,
            expected_start,
            expected_end
        )

        # Allow breakpoint tolerance to rescue a near-boundary event,
        # but unrelated intervals must not match.
        if overlap == 0:

            expanded_start = (
                expected_start
                - breakpoint_tolerance
            )

            expanded_end = (
                expected_end
                + breakpoint_tolerance
            )

            overlap = interval_overlap(
                observed_start,
                observed_end,
                expanded_start,
                expanded_end
            )

            if overlap == 0:
                return None

        expected_span = (
            expected_end
            - expected_start
            + 1
        )

        observed_span = (
            observed_end
            - observed_start
            + 1
        )

        expected_overlap_fraction = (
            overlap /
            expected_span
        )

        observed_overlap_fraction = (
            overlap /
            observed_span
        )

        midpoint_observed = (
            observed_start
            + observed_end
        ) / 2

        midpoint_expected = (
            expected_start
            + expected_end
        ) / 2

        midpoint_distance = abs(
            midpoint_observed
            - midpoint_expected
        )

        # Higher tuple is better.
        return (
            expected_overlap_fraction,
            observed_overlap_fraction,
            -length_difference,
            -midpoint_distance
        )


   ####################
    # Insertion

    elif event_type == "insertion":
        anchor = observed_start
        expanded_start = (
            expected_start
            - breakpoint_tolerance
        )

        expanded_end = (
            expected_end
            + breakpoint_tolerance
        )

        if not (
            expanded_start
            <= anchor
            <= expanded_end
        ):
            return None

        if expected_start <= anchor <= expected_end:
            anchor_distance = 0

        else:
            anchor_distance = min(
                abs(anchor - expected_start),
                abs(anchor - expected_end)
            )

        return (
            1.0,
            1.0,
            -length_difference,
            -anchor_distance
        )

    else:

        raise ValueError(
            f"Unsupported structural event type: "
            f"{event_type}"
        )

## best SV assignment by alignment score 
def best_rule_score(
    observed_events,
    expected_event,
    length_tolerance,
    breakpoint_tolerance
):

    expected_count = int(
        expected_event.get(
            "count",
            1
        )
    )

    scored_events = []

    for observed_event in observed_events:
        score = score_event(
            observed_event,
            expected_event,
            length_tolerance,
            breakpoint_tolerance
        )

        if score is not None:
            scored_events.append(score)

    if len(scored_events) < expected_count:
        return None

    scored_events.sort(
        reverse=True
    )

    selected = scored_events[
        :expected_count
    ]

    # Current FCGR configurations use count = 1.
    # For count > 1, use the weakest selected event so all required
    # events must be compatible.
    return min(selected)

### haplotype classifier
def classify(row, config):

    ### Only classify complete assemblies
    if row["assembly_status"] != "contiguous":
        return row["assembly_status"]

    ### Observed structural events
    insertions = parse_events(
        row["large_insertion_positions"]
    )

    deletions = parse_events(
        row["large_deletion_positions"]
    )

    ### Configuration
    defaults = config.get(
        "defaults",
        {}
    )

    breakpoint_tolerance = int(
        defaults.get(
            "breakpoint_tolerance",
            0
        )
    )

    length_tolerance = int(
        defaults.get(
            "length_tolerance",
            0
        )
    )


    ### No large structural events = reference-like haplotype
    if not insertions and not deletions:
        for haplotype, rules in config["haplotypes"].items():
            structural_events = rules.get(
                "structural_events",
                []
            )

            if not structural_events:
                return haplotype
        return "novel"

    ### Score non-reference structural haplotypes
    candidates = []

    for haplotype, rules in config["haplotypes"].items():
        structural_events = rules.get(
            "structural_events",
            []
        )

        if not structural_events:
            continue

        if len(structural_events) != 1:
            continue

        expected_event = structural_events[0]
        event_type = expected_event["type"]

        if event_type == "insertion":

            # Simple insertion haplotypes should not also contain
            # a large deletion.
            if deletions:
                continue

            score = best_rule_score(
                insertions,
                expected_event,
                length_tolerance,
                breakpoint_tolerance
            )

        elif event_type == "deletion":

            # Simple deletion haplotypes should not also contain
            # a large insertion.
            if insertions:
                continue

            score = best_rule_score(
                deletions,
                expected_event,
                length_tolerance,
                breakpoint_tolerance
            )

        else:
            raise ValueError(
                f"Unsupported structural event type "
                f"'{event_type}' for haplotype "
                f"'{haplotype}'."
            )

        if score is not None:
            candidates.append(
                (
                    score,
                    haplotype
                )
            )


   ### No known structural haplotype matched
    if not candidates:
        return "novel"


    ### Choose strongest structural match
    candidates.sort(
        key=lambda x: x[0],
        reverse=True
    )

    return candidates[0][1]

##########################################################################
# Classify contigs as haplotypes based on alignment to reference genome

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Classify assembled haplotypes using "
            "user-defined structural configurations."
        )
    )

    parser.add_argument("--summary",required=True)
    parser.add_argument("--config",required=True)
    parser.add_argument("--output",required=True)
    args = parser.parse_args()

    ### Load configuration
    with open(args.config) as handle:
        config = yaml.safe_load(
            handle
        )

    if not isinstance(config, dict):
        raise ValueError(
            "Haplotype configuration must be a YAML mapping."
        )

    if "haplotypes" not in config:
        raise ValueError(
            "Haplotype configuration is missing "
            "the 'haplotypes' section."
        )

    ### Read alignment summary
    df = pd.read_csv(
        args.summary,
        sep="\t",
        keep_default_na=False
    )

    required_columns = [
        "assembly_status",
        "large_insertion_positions",
        "large_deletion_positions"
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in df.columns
    ]

    if missing_columns:

        raise ValueError(
            "Summary file is missing required column(s): "
            + ", ".join(missing_columns)
        )


    # Classify assemblies using known haplotype coordinate from yaml
    df["structural_haplotype"] = df.apply(
        lambda row: classify(
            row,
            config
        ),
        axis=1
    )


    ### Write
    df.to_csv(
        args.output,
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    main()
