/*
###################################################################################
Module: COMBINE_QC
###################################################################################

Combine individual sample coverage QC results into a single summary TSV.

Input: All sample coverage QC TSV files

Output: coverage_qc_summary.tsv

###################################################################################
*/

process COMBINE_QC {

    tag "coverage QC summary"

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path qc_results

    output:
    path "coverage_qc_summary.tsv", emit: summary

    script:
    """
    ${projectDir}/bin/combine_qc.py \
        --input ${qc_results} \
        --output coverage_qc_summary.tsv

    test -s coverage_qc_summary.tsv
    """
}
