/*
###################################################################################
Module: BAM_TO_FASTQ
###################################################################################

 Convert the regional BAM file to FASTQ.

 Input: sample metadata (sample ID, BAM file, gene name, gene coordinates, flank size) 
	from config/samples.tsv

 Output: subset region-specific reads in fastq file format

###################################################################################
*/

process BAM2FASTQ {

    tag "${meta.sample} | ${meta.gene} | ${meta.run_name}"

    cpus 2
    memory '8 GB'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.fastq"),
          emit: fastq

    script:

    def prefix =
        "${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}"

    """
    samtools fastq -@ ${task.cpus} ${bam} > ${prefix}.fastq

    test -s ${prefix}.fastq
    """
}
