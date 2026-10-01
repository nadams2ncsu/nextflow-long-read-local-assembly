/*
################################################################################
Module: MERGE_ASSEMBLY_DIFFERENCES
################################################################################

Merge hap1 and hap2 assembly-reference differences for each paired-contiguous
sample + group + gene + flank assembly.

################################################################################
*/

process MERGE_ASSEMBLY_DIFFERENCES {

    tag "${meta.group} | ${meta.sample} | ${meta.gene} | ${meta.run_name}"

    cpus 1
    memory '4 GB'

    publishDir "${params.outdir}/differences", mode: 'copy'

    input:

    tuple val(meta),
          path(hap1_bed),
          path(hap2_bed)

    output:

    tuple val(meta),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.differences.bed"),
          emit: differences

    script:

    def prefix =
        "${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}"

    """
    ${projectDir}/bin/merge_assembly_differences.py \
        --hap1 ${hap1_bed} \
        --hap2 ${hap2_bed} \
        --output ${prefix}.differences.bed

    test -f ${prefix}.differences.bed
    """
}
