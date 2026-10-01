/*
################################################################################
Module: COMPARE_GROUP_DIFFERENCES
################################################################################

Compare merged assembly-reference difference BED files across groups.
	- For the same sample + gene, every available flank combination is compared
		between different groups.

Outputs: nucleotide differences between individuals sequenced in 2 groups (TSV & BED)

################################################################################
*/

process COMPARE_GROUP_DIFFERENCES {

    tag "across-group difference comparison"

    cpus 1
    memory '4 GB'

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:

    val comparison_metadata

    path difference_beds

    output:

    path "across_groups_difference_comparison.tsv",
        emit: comparison

    path "across_groups_difference_discordance.bed",
        emit: discordance

    script:
	
    def manifest_rows = comparison_metadata.join("\n")

    """
    printf "sample\\tgroup\\tgene\\tflank_kb\\tbed\\n" \
        > difference_manifest.tsv

    cat >> difference_manifest.tsv <<'EOF'
${manifest_rows}
EOF

    ${projectDir}/bin/compare_group_differences.py \
        --manifest difference_manifest.tsv \
        --summary-output across_groups_difference_comparison.tsv \
        --discordance-output across_groups_difference_discordance.bed

    test -f across_groups_difference_comparison.tsv
    test -f across_groups_difference_discordance.bed
    """
}
