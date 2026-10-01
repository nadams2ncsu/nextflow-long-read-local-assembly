/*
####################################################################################
Module: COMBINE_RESULTS
####################################################################################

Combine all per-assembly classified summary TSV files into one final
workflow summary.

Input: sample alignment summary tsv files

Output:
    local_assembly_summary.tsv ==> all assembly alignment results
    contiguous_assembly_summary.tsv ==> only contiguous assembly alignment results

####################################################################################
*/

process COMBINE_RESULTS {

    tag "final summary"

    cpus 1
    memory '4 GB'

    publishDir "${params.outdir}/summary", mode: 'copy'

    input:

    path classified_results

    output:

    path "local_assembly_summary.tsv",
        emit: summary

    path "contiguous_assembly_summary.tsv",
         emit: contiguous_summary

    script:

    """
    ${projectDir}/bin/combine_results.py \
        --input ${classified_results} \
        --output local_assembly_summary.tsv \
	--contiguous-output contiguous_assembly_summary.tsv

    test -s local_assembly_summary.tsv
    test -s contiguous_assembly_summary.tsv
    """
}
