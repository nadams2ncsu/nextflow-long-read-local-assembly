/*
##############################################################################################
Module: EXTRACT_READS
##############################################################################################

 Extract reads mapping to the sample-specific genomic region.

 Input: sample metadata (config/samples.tsv)

 Output: Generate BAM file subsets from user-supplied coordinate regions and flank sizes

############################################################################################
*/

process EXTRACT_READS {

    tag "${meta.sample} | ${meta.gene} | ${meta.run_name}"

    cpus 2
    memory '8 GB'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.bam"),
          emit: bam

    script:

    def prefix =
        "${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}"

    """
    samtools view \
        -@ ${task.cpus} \
        -b \
        ${bam} \
        "${meta.region}" \
        > ${prefix}.bam
    """
}
