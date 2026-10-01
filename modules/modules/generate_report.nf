/*
#################################################################################
Module: GENERATE_REPORT
#################################################################################

Generate a self-contained HTML summary report for the workflow.

The report summarizes:
    - assembly genotype / outcome
    - coverage QC
    - regional read statistics
    - representative across-group sequence comparisons
    - representative across-group assembly-reference difference comparisons

#################################################################################
*/

process GENERATE_REPORT {

    tag "HTML workflow report"

    publishDir "${params.outdir}/report", mode: 'copy'


    input:

    path assembly_summary
    path coverage_qc
    path read_stats
    path comparison_summary
    path difference_comparison


    output:

    path "local_assembly_report.html",
        emit: report


    script:

    """
    ${projectDir}/bin/generate_report.py \
        --assembly-summary ${assembly_summary} \
        --coverage-qc ${coverage_qc} \
        --read-stats ${read_stats} \
        --comparison-summary ${comparison_summary} \
        --difference-comparison ${difference_comparison} \
        --output local_assembly_report.html

    test -s local_assembly_report.html
    """
}
