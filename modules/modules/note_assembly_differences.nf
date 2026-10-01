/*
################################################################################
Module: NOTE_ASSEMBLY_DIFFERENCES
################################################################################

Identify sequence differences between a contiguous haplotype assembly and the
reference genome from the assembly-to-reference BAM alignment.

For loci with a haplotype configuration, evaluation_coordinates from the
haplotype YAML are used because the assembled contigs are aligned against
the common workflow reference.

Input: reference genome & aligned contiguous contigs

Output: assembly-specific BED files noting differences with the reference genome

################################################################################
*/

process NOTE_ASSEMBLY_DIFFERENCES {

    tag "${meta.group} | ${meta.sample} | ${meta.gene} | ${meta.run_name} | ${meta.haplotype}"

    cpus 1
    memory '4 GB'

    input:

    tuple val(meta),
          path(bam),
          path(bai)

    path reference


    output:

    tuple val(meta),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}_${meta.haplotype}.differences.bed"),
          emit: differences


    script:

    def prefix =
        "${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}_${meta.haplotype}"

    def haplotype_arg =
        meta.haplotype_config.toUpperCase() != 'NA'
            ? "--haplotype-config ${meta.haplotype_config}"
            : ""

    """
    ${projectDir}/bin/note_assembly_differences.py \
        --bam ${bam} \
        --reference ${reference} \
        --haplotype "${meta.haplotype}" \
        --coordinates "${meta.coordinates}" \
        ${haplotype_arg} \
        --large-indel-threshold ${params.large_indel_threshold} \
        --output ${prefix}.differences.bed

    test -f ${prefix}.differences.bed
    """
}
