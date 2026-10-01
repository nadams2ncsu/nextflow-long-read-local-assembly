#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import plotly.express as px
import plotly.io as pio
from plotly.offline import get_plotlyjs


################################################################################
# Generate HTML analysis report
#
# Summarizes:
#   - assembly genotype / outcome
#   - coverage QC
#   - regional read statistics
#   - representative across-group sequence similarity
#   - representative across-group assembly-reference difference concordance
#
# Input: all summary tsv files
# Outputs: HTML report
#

#######################
# helper functions

## read summary tsv files
def read_tsv(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()

    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

## html plotter
def plot_html(fig):
    fig.update_layout(
        template="plotly_white",
        margin=dict(l=60, r=30, t=60, b=60)
    )

    return pio.to_html(
        fig,
        full_html=False,
        include_plotlyjs=False
    )


## html tabler
def table_html(df, max_rows=100):
    if df.empty:
        return "<p>No data available.</p>"

    return df.head(max_rows).to_html(
        index=False,
        classes="data-table",
        border=0
    )

### formatting functions
def card(label, value):
    return f"""
    <div class="card">
        <div class="card-value">{value}</div>
        <div class="card-label">{label}</div>
    </div>
    """


def format_flank(value):
    try:
        value = int(float(value))

    except (TypeError, ValueError):
        return str(value)
    if value == 1000:
        return "1 Mb"
    return f"{value} kb"

## sort genotype combinations
def genotype_sort_key(haplotype):
    text = str(haplotype)

    if text.startswith("H"):
        try:
            return int(text[1:])

        except ValueError:
            pass

    return 999


### Build one genotype / outcome for each sample + group + gene + flank
def build_assembly_results(df):
    required = {
        "sample",
        "group",
        "gene",
        "flank_kb",
        "assembly_haplotype",
        "assembly_status",
        "structural_haplotype"
    }

    if df.empty or not required.issubset(df.columns):
        return pd.DataFrame()

    rows = []

    group_columns = [
        "sample",
        "group",
        "gene",
        "flank_kb"
    ]

    for keys, group_df in df.groupby(
        group_columns,
        dropna=False
    ):

        sample, group, gene, flank_kb = keys

        hap1 = group_df[
            group_df["assembly_haplotype"] == "hap1"
        ]

        hap2 = group_df[
            group_df["assembly_haplotype"] == "hap2"
        ]

        genotype = None

        # Both haplotypes must exist and be contiguous
        if len(hap1) == 1 and len(hap2) == 1:

            status1 = str(
                hap1.iloc[0]["assembly_status"]
            )

            status2 = str(
                hap2.iloc[0]["assembly_status"]
            )

            if (
                status1 == "contiguous"
                and status2 == "contiguous"
            ):

                structural1 = hap1.iloc[0][
                    "structural_haplotype"
                ]

                structural2 = hap2.iloc[0][
                    "structural_haplotype"
                ]

                if (
                    pd.notna(structural1)
                    and pd.notna(structural2)
                    and str(structural1)
                    not in {"NA", "unphased"}
                    and str(structural2)
                    not in {"NA", "unphased"}
                ):

                    pair = sorted(
                        [
                            str(structural1),
                            str(structural2)
                        ],
                        key=genotype_sort_key
                    )

                    genotype = "/".join(pair)

                else:
                    genotype = "unphased"

        # If a paired contiguous genotype was not obtained,
        # summarize the assembly outcome.
        if genotype is None:

            statuses = set(
                group_df[
                    "assembly_status"
                ]
                .dropna()
                .astype(str)
            )

            if "fragmented" in statuses:
                genotype = "fragmented"

            elif "discontiguous" in statuses:
                genotype = "discontiguous"

            elif "unphased" in statuses:
                genotype = "unphased"

            elif statuses:
                genotype = "/".join(
                    sorted(statuses)
                )

            else:
                genotype = "unphased"

        rows.append(
            {
                "sample": sample,
                "group": group,
                "gene": gene,
                "flank_kb": flank_kb,
                "genotype_or_outcome": genotype
            }
        )

    return pd.DataFrame(rows)


### Collapse identical genotype calls across flank sizes
def collapse_assembly_results(df):
    if df.empty:
        return df

    rows = []

    for keys, group_df in df.groupby(
        [
            "sample",
            "group",
            "gene",
            "genotype_or_outcome"
        ],
        dropna=False
    ):

        sample, group, gene, genotype = keys

        flank_values = sorted(
            group_df[
                "flank_kb"
            ].dropna().unique(),
            key=lambda x: float(x)
        )

        flank_labels = ", ".join(
            format_flank(x)
            for x in flank_values
        )

        rows.append(
            {
                "sample": sample,
                "group": group,
                "gene": gene,
                "genotype_or_outcome": genotype,
                "flanks": flank_labels
            }
        )

    result = pd.DataFrame(rows)

    return result.sort_values(
        [
            "sample",
            "group",
            "gene",
            "genotype_or_outcome"
        ]
    ).reset_index(drop=True)


### Comparison label
def add_comparison_label(df):
    required = {
        "sample",
        "gene",
        "group1",
        "group1_flank_kb",
        "group2",
        "group2_flank_kb"
    }

    if df.empty or not required.issubset(df.columns):
        return df

    df = df.copy()

    df["comparison"] = (
        df["sample"].astype(str)
        + " | "
        + df["gene"].astype(str)
        + " | "
        + df["group1"].astype(str)
        + " "
        + df["group1_flank_kb"].apply(
            format_flank
        )
        + " vs "
        + df["group2"].astype(str)
        + " "
        + df["group2_flank_kb"].apply(
            format_flank
        )
    )

    return df


# Select one representative comparison per sample + gene + group pair
#
# Priority:
#   1. Same flank size in both groups
#   2. Smallest difference between flank sizes
#   3. Smallest maximum flank size
#   4. Smallest group1 flank
#   5. Smallest group2 flank

def select_representative_comparisons(df):
    required = {
        "sample",
        "gene",
        "group1",
        "group1_flank_kb",
        "group2",
        "group2_flank_kb"
    }

    if df.empty or not required.issubset(df.columns):
        return df

    data = df.copy()

    data["group1_flank_kb"] = pd.to_numeric(
        data["group1_flank_kb"],
        errors="coerce"
    )

    data["group2_flank_kb"] = pd.to_numeric(
        data["group2_flank_kb"],
        errors="coerce"
    )

    data["same_flank"] = (
        data["group1_flank_kb"]
        == data["group2_flank_kb"]
    )

    data["flank_difference"] = (
        data["group1_flank_kb"]
        - data["group2_flank_kb"]
    ).abs()

    data["largest_flank"] = data[
        [
            "group1_flank_kb",
            "group2_flank_kb"
        ]
    ].max(axis=1)

    selected = []

    group_columns = [
        "sample",
        "gene",
        "group1",
        "group2"
    ]

    for _, group_df in data.groupby(
        group_columns,
        dropna=False
    ):

        comparison_options = (
            group_df[
                [
                    "group1_flank_kb",
                    "group2_flank_kb",
                    "same_flank",
                    "flank_difference",
                    "largest_flank"
                ]
            ]
            .drop_duplicates()
            .sort_values(
                [
                    "same_flank",
                    "flank_difference",
                    "largest_flank",
                    "group1_flank_kb",
                    "group2_flank_kb"
                ],
                ascending=[
                    False,
                    True,
                    True,
                    True,
                    True
                ]
            )
        )

        best = comparison_options.iloc[0]

        subset = group_df[
            (
                group_df["group1_flank_kb"]
                == best["group1_flank_kb"]
            )
            &
            (
                group_df["group2_flank_kb"]
                == best["group2_flank_kb"]
            )
        ].copy()

        selected.append(
            subset
        )

    if not selected:
        return pd.DataFrame()

    result = pd.concat(
        selected,
        ignore_index=True
    )

    return result.drop(
        columns=[
            "same_flank",
            "flank_difference",
            "largest_flank"
        ]
    )


################################################################################
# Generate HTML report for local assembly workflow

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate HTML report for localized "
            "long-read assembly."
        )
    )

    parser.add_argument(
        "--assembly-summary",
        required=True
    )

    parser.add_argument(
        "--coverage-qc",
        required=True
    )

    parser.add_argument(
        "--read-stats",
        required=True
    )

    parser.add_argument(
        "--comparison-summary",
        required=True
    )

    parser.add_argument(
        "--difference-comparison",
        required=True
    )

    parser.add_argument(
        "--output",
        required=True
    )

    args = parser.parse_args()


    ### Read inputs

    assembly = read_tsv(
        args.assembly_summary
    )

    coverage = read_tsv(
        args.coverage_qc
    )

    reads = read_tsv(
        args.read_stats
    )

    comparisons = read_tsv(
        args.comparison_summary
    )

    difference_comparisons = read_tsv(
        args.difference_comparison
    )


    ### Prepare assembly genotype / outcome summary

    assembly_by_flank = build_assembly_results(
        assembly
    )

    assembly_results = collapse_assembly_results(
        assembly_by_flank
    )

    sections = []


    ### Run overview

    samples = (
        assembly["sample"].nunique()
        if "sample" in assembly.columns
        else 0
    )

    groups = (
        assembly["group"].nunique()
        if "group" in assembly.columns
        else 0
    )

    genes = (
        assembly["gene"].nunique()
        if "gene" in assembly.columns
        else 0
    )

    if (
        not assembly.empty
        and {"sample", "group"}.issubset(
            assembly.columns
        )
    ):

        sample_groups = (
            assembly[
                [
                    "sample",
                    "group"
                ]
            ]
            .drop_duplicates()
            .shape[0]
        )

    else:
        sample_groups = 0


    cards = [
        card(
            "Samples",
            samples
        ),

        card(
            "Sample-group datasets",
            sample_groups
        ),

        card(
            "Groups",
            groups
        ),

        card(
            "Genes",
            genes
        )
    ]


    sections.append(
        """
        <section>

            <h2>Run Overview</h2>

            <div class="cards">
                {}
            </div>

        </section>
        """.format(
            "".join(cards)
        )
    )


    ### Assembly results

    sections.append(
        """
        <section>

            <h2>Assembly Results</h2>

            <p>
            Structural haplotypes from paired-contiguous hap1 and hap2
            assemblies are combined into genotypes. Identical genotype or
            assembly-outcome calls across flank sizes are collapsed into one
            sample-group result.
            </p>

            {}

        </section>
        """.format(
            table_html(
                assembly_results
            )
        )
    )


    ### Genotype / outcome distribution

    if not assembly_results.empty:

        genotype_counts = (
            assembly_results[
                "genotype_or_outcome"
            ]
            .value_counts()
            .rename_axis(
                "genotype_or_outcome"
            )
            .reset_index(
                name="count"
            )
        )

        fig = px.bar(
            genotype_counts,
            x="genotype_or_outcome",
            y="count",
            title=(
                "Structural Genotype / "
                "Assembly Outcome"
            ),
            labels={
                "genotype_or_outcome":
                    "Genotype / assembly outcome",

                "count":
                    "Number of sample-group results"
            }
        )

        sections.append(
            "<section>"
            "<h2>Structural Genotype Summary</h2>"
            + plot_html(fig)
            + "</section>"
        )


    ### Coverage QC

    coverage_display = coverage.copy()

    if (
        not coverage_display.empty
        and {"sample", "group"}.issubset(
            coverage_display.columns
        )
    ):

        coverage_display = (
            coverage_display.sort_values(
                [
                    "sample",
                    "group"
                ]
            )
        )

    sections.append(
        """
        <section>

            <h2>Coverage QC</h2>

            <p>
            Coverage QC is reported separately for each sample-group dataset.
            </p>

            {}

        </section>
        """.format(
            table_html(
                coverage_display
            )
        )
    )


    ### Read statistics QC

    reads_display = reads.copy()

    if (
        not reads_display.empty
        and {"sample", "group"}.issubset(
            reads_display.columns
        )
    ):

        reads_display = (
            reads_display.sort_values(
                [
                    "sample",
                    "group"
                ]
            )
        )

    sections.append(
        """
        <section>

            <h2>Read Statistics QC</h2>

            <p>
            Regional read statistics are reported separately for each
            sample-group dataset.
            </p>

            {}

        </section>
        """.format(
            table_html(
                reads_display
            )
        )
    )


    ############################################################################
    # Across-group sequence comparison
    #
    # One representative flank comparison is retained for each sample,
    # gene, and group pair. The individual haplotype rows belonging to
    # that comparison are retained in the table.

    if (
        not comparisons.empty
        and "sequence_similarity"
        in comparisons.columns
    ):

        comparison_data = (
            select_representative_comparisons(
                comparisons
            )
        )

        comparison_data = add_comparison_label(
            comparison_data
        )

        sections.append(
            """
            <section>

                <h2>Across-Group Sequence Comparison</h2>

                <p>
                One representative assembly comparison is shown for each
                sample, gene, and pair of groups. Comparisons using the same
                flank size are preferred when available.
                </p>

                {}

            </section>
            """.format(
                table_html(
                    comparison_data
                )
            )
        )

    else:

        sections.append(
            """
            <section>

                <h2>Across-Group Sequence Comparison</h2>

                <p>
                No samples had paired-contiguous assemblies available
                in two or more groups. No across-group sequence
                comparisons were performed.
                </p>

            </section>
            """
        )


    ############################################################################
    # Across-group assembly-reference difference comparison

    required_difference_columns = {
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
    }

    if (
        not difference_comparisons.empty
        and required_difference_columns.issubset(
            difference_comparisons.columns
        )
    ):

        ########################################################################
        # Select one representative comparison per sample + gene + group pair

        difference_data = (
            select_representative_comparisons(
                difference_comparisons
            )
        )

        difference_data = add_comparison_label(
            difference_data
        )


        ########################################################################
        # Difference concordance

        fig = px.bar(
            difference_data,
            x="comparison",
            y="difference_concordance",
            text="difference_concordance",
            hover_data=[
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
                "union_difference_count"
            ],
            title=(
                "Across-Group Assembly-Reference "
                "Difference Concordance"
            ),
            labels={
                "comparison":
                    "Assembly comparison",

                "difference_concordance":
                    "Difference concordance (%)"
            }
        )

        fig.update_traces(
            texttemplate="%{text:.1f}%",
            textposition="inside"
        )

        fig.update_yaxes(
            range=[0, 100]
        )


        sections.append(
            "<section>"
            "<h2>Across-Group Assembly-Reference "
            "Difference Comparison</h2>"
            "<p>"
            "Merged assembly-reference differences are compared "
            "between independently generated assemblies for the "
            "same sample and gene. Differences are considered "
            "shared when chromosome, position, difference type, "
            "reference sequence, and assembly sequence are "
            "identical. One representative flank comparison is "
            "shown for each sample, gene, and pair of groups."
            "</p>"
            + plot_html(fig)
        )


        ########################################################################
        # Shared / group-specific difference percentages

        difference_counts = difference_data[
            [
                "comparison",
                "shared_difference_count",
                "group1_only_count",
                "group2_only_count",
                "union_difference_count"
            ]
        ].copy()

        difference_counts = difference_counts.melt(
            id_vars=[
                "comparison",
                "union_difference_count"
            ],
            var_name="difference_category",
            value_name="count"
        )

        difference_counts[
            "difference_category"
        ] = (
            difference_counts[
                "difference_category"
            ]
            .replace(
                {
                    "shared_difference_count":
                        "Shared",

                    "group1_only_count":
                        "Group 1 only",

                    "group2_only_count":
                        "Group 2 only"
                }
            )
        )


        ########################################################################
        # Convert raw counts to percentage of union

        difference_counts["percentage"] = (
            difference_counts["count"]
            / difference_counts[
                "union_difference_count"
            ]
            * 100
        )


        ########################################################################
        # Raw count labels
        #
        # All categories remain in the stacked bar.
        #
        # Shared:
        #   Always display the raw count.
        #
        # Group-specific:
        #   Display the raw count only when >= 100.

        difference_counts["count_label"] = (
            difference_counts.apply(
                lambda row: (
                    str(int(row["count"]))
                    if (
                        row["difference_category"] == "Shared"
                        or row["count"] >= 100
                    )
                    else ""
                ),
                axis=1
            )
        )


        ########################################################################
        # Plot percentage while displaying selected raw counts inside bars

        fig = px.bar(
            difference_counts,
            x="comparison",
            y="percentage",
            color="difference_category",
            text="count_label",
            barmode="stack",
            title=(
                "Shared and Discordant "
                "Assembly-Reference Differences"
            ),
            labels={
                "comparison":
                    "Assembly comparison",

                "percentage":
                    "Percentage of differences (%)",

                "difference_category":
                    "Difference category"
            },
            hover_data={
                "count": True,
                "percentage": ":.2f",
                "union_difference_count": True,
                "count_label": False
            }
        )

        fig.update_traces(
            texttemplate="%{text}",
            textposition="inside"
        )

        fig.update_yaxes(
            range=[0, 100]
        )


        sections.append(
            plot_html(fig)
        )


        ### Difference comparison table

        sections.append(
            table_html(
                difference_data
            )
        )

        sections.append(
            "</section>"
        )

    else:

        sections.append(
            """
            <section>

                <h2>
                    Across-Group Assembly-Reference Difference Comparison
                </h2>

                <p>
                No paired-contiguous assemblies were available for
                across-group assembly-reference difference comparison.
                </p>

            </section>
            """
        )


    ##############
    # HTML

    plotly_js = get_plotlyjs()

    html = f"""
    <!DOCTYPE html>

    <html lang="en">

    <head>

        <meta charset="UTF-8">

        <meta
            name="viewport"
            content="width=device-width, initial-scale=1.0"
        >

        <title>
            Localized Long-Read Assembly Report
        </title>

        <script>
        {plotly_js}
        </script>

        <style>

            body {{
                font-family: Arial, sans-serif;
                margin: 0;
                background: #f5f5f5;
                color: #222;
            }}

            header {{
                background: #ffffff;
                padding: 30px 5%;
                border-bottom: 1px solid #ddd;
            }}

            header h1 {{
                margin: 0;
            }}

            header p {{
                color: #666;
                margin-bottom: 0;
            }}

            main {{
                width: 90%;
                max-width: 1400px;
                margin: 30px auto;
            }}

            section {{
                background: white;
                padding: 25px;
                margin-bottom: 25px;
                border-radius: 8px;
                box-shadow: 0 1px 4px rgba(0,0,0,0.08);
                overflow-x: auto;
            }}

            .cards {{
                display: flex;
                flex-wrap: wrap;
                gap: 15px;
            }}

            .card {{
                min-width: 160px;
                padding: 20px;
                border: 1px solid #ddd;
                border-radius: 6px;
                background: #fafafa;
            }}

            .card-value {{
                font-size: 28px;
                font-weight: bold;
            }}

            .card-label {{
                margin-top: 5px;
                color: #666;
            }}

            .data-table {{
                width: 100%;
                border-collapse: collapse;
                font-size: 13px;
                margin-bottom: 20px;
            }}

            .data-table th,
            .data-table td {{
                padding: 8px;
                border: 1px solid #ddd;
                text-align: left;
            }}

            .data-table th {{
                background: #f0f0f0;
            }}

            h2 {{
                margin-top: 0;
            }}

            p {{
                line-height: 1.5;
            }}

        </style>

    </head>

    <body>

        <header>

            <h1>
                Localized Long-Read Assembly Report
            </h1>

            <p>
                Assembly genotype, sequencing QC,
                cross-group sequence concordance, and
                assembly-reference difference summary
            </p>

        </header>

        <main>

            {''.join(sections)}

        </main>

    </body>

    </html>
    """

    Path(
        args.output
    ).write_text(
        html,
        encoding="utf-8"
    )


if __name__ == "__main__":
    main()
