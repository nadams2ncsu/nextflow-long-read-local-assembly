/*
################################################################################
Module: ACROSS_GROUPS_COMPARISON
################################################################################

Compare paired-contiguous haplotype assemblies across groups.

- All valid paired-contiguous flank combinations are compared between groups
  for the same sample + gene.

Input:
    - metadata for paired-contiguous assemblies
    - BAM files
    - BAI files

Output:
    across_groups_comparison.tsv

################################################################################
*/

process ACROSS_GROUPS_COMPARISON {

    tag "across-group comparison"

    cpus 2
    memory '8 GB'

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:

    val comparison_metadata

    path bams

    path bais

    output:

    path "across_groups_comparison.tsv",
        emit: comparison

    script:

    def manifest_rows = comparison_metadata.join("\n")

    """
    printf "sample\\tgroup\\tgene\\tflank_kb\\tassembly_haplotype\\tbam\\n" \
        > comparison_manifest.tsv

    cat >> comparison_manifest.tsv <<'EOF'
${manifest_rows}
EOF

    ${projectDir}/bin/across_groups_comparison.py \
        --manifest comparison_manifest.tsv \
        --similarity-threshold ${params.group_similarity_threshold} \
        --output across_groups_comparison.tsv

    test -f across_groups_comparison.tsv
    """
}
