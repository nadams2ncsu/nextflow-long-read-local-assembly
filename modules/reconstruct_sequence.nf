/*
###################################################################################
Module: RECONSTRUCT_SEQUENCE
###################################################################################

Reconstruct paired contiguous haplotype sequences from aligned BAM files and
trim them to the workflow-reference evaluation coordinates.

This process only receives sample + group + gene + flank combinations where
BOTH hap1 and hap2 were classified as contiguous.
###################################################################################
*/

process RECONSTRUCT_SEQUENCE {

    tag "${meta.sample} | ${meta.group} | ${meta.gene} | ${meta.run_name}"

    publishDir "${params.outdir}/fasta", mode: 'copy'

    input:

    tuple val(meta),
          path(hap1_bam),
          path(hap1_bai),
          path(hap2_bam),
          path(hap2_bai)

    output:

    tuple val(meta),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.hap1.fa"),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.hap2.fa"),
          emit: fasta

    script:

    def prefix =
        "${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}"

    def haplotype_arg =
        meta.haplotype_config.toUpperCase() != 'NA'
            ? "--haplotype-config ${meta.haplotype_config}"
            : ""

    """
    ${projectDir}/bin/reconstruct_sequence.py \
        --bam ${hap1_bam} \
        --coordinates "${meta.coordinates}" \
        ${haplotype_arg} \
        --output ${prefix}.hap1.fa

    ${projectDir}/bin/reconstruct_sequence.py \
        --bam ${hap2_bam} \
        --coordinates "${meta.coordinates}" \
        ${haplotype_arg} \
        --output ${prefix}.hap2.fa

    test -s ${prefix}.hap1.fa
    test -s ${prefix}.hap2.fa
    """
}
