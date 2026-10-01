/*
###################################################################################
Module:  ALIGN_CONTIGS
###################################################################################

 Align assembled haplotype contigs to the reference genome using minimap2.

 Default preset:
   asm5

 The preset can be overridden using:
   --minimap2_preset PRESET

 Additional minimap2 arguments can be supplied using:
   --minimap2_args "ARGS"

 --eqx is always enabled because downstream SNV counting requires
 = and X CIGAR operations.

 Input: run name, contige haplotype label from hifiasm, reference genome
	OPTIONAL minimap alignment parameters

 Output: aligned, unsorted contigs to the reference genome in BAM file format

###################################################################################
*/

process ALIGN_CONTIGS {

    tag "${meta.sample} | ${meta.gene} | ${meta.run_name} | ${meta.haplotype}"

    cpus 8
    memory '16 GB'

    input:

    tuple val(meta), path(fasta)

    path reference


    output:

    tuple val(meta),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.${meta.haplotype}.unsorted.bam"),
          emit: bam


    script:

    def prefix =
        "${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.${meta.haplotype}"

    """
    minimap2 -t ${task.cpus} -ax ${params.minimap2_preset} --eqx ${reference} ${fasta} \
        | samtools view -@ ${task.cpus} -b -o ${prefix}.unsorted.bam

    samtools quickcheck ${prefix}.unsorted.bam
    """
}
