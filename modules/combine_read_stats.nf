/*
################################################################################
Module: COMBINE_READ_STATS
################################################################################

Combine individual regional read-statistic files into a single summary.

Input:  Per-sample READ_STATS TSV files

Output: read_stat_summary.tsv

################################################################################
*/

process COMBINE_READ_STATS {

    tag "combine read statistics"

    cpus 1
    memory '4 GB'

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:

    path read_stats

    output:

    path "read_stat_summary.tsv",
        emit: summary

    script:

    """
    ${projectDir}/bin/combine_read_stats.py \
        --input ${read_stats} \
        --output read_stat_summary.tsv

    test -s read_stat_summary.tsv
    """
}
