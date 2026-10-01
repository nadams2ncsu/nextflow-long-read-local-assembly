/*
###################################################################################
Module: SORT_INDEX
###################################################################################

 Sort and index assembled contig alignments.

 Input: Unsorted, aligned contigs to the reference genome

 Output: Aligned, sorted contigs to the reference genome

###################################################################################
*/

process SORT_INDEX {

    tag "${meta.sample} | ${meta.gene} | ${meta.run_name} | ${meta.haplotype}"

    cpus 4
    memory '8 GB'

    input:

    tuple val(meta), path(bam)

    output:

    tuple val(meta),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.${meta.haplotype}.bam"),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.${meta.haplotype}.bam.bai"),
          emit: bam


    script:

    def prefix =
        "${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.${meta.haplotype}"

    """
    samtools sort \
        -@ ${task.cpus} \
        -o ${prefix}.bam \
        ${bam}

    samtools index \
        -@ ${task.cpus} \
        ${prefix}.bam

    samtools quickcheck ${prefix}.bam
    test -s ${prefix}.bam.bai
    """
}
