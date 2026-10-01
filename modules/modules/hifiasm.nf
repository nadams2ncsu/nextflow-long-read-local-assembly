/*
###################################################################################
Module: HIFIASM
###################################################################################

 Perform localized long-read assembly using hifiasm.

 PacBio HiFi: Default hifiasm execution
 Oxford Nanopore: Adds the --ont option

 Additional hifiasm arguments can be suplied using: --hifiasm_args "ARGS"

 Input: Gene & flank size-specific fastq files for each sample

 Output: Gene & flank size-specific assembly graphs for each sample haplotype

###################################################################################
*/

process HIFIASM {

    tag "${meta.sample} | ${meta.gene} | ${meta.run_name} | ${meta.data_type}"

    cpus 16
    memory '64 GB'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.bp.hap1.p_ctg.gfa"),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.bp.hap2.p_ctg.gfa"),
          emit: gfa

    script:

    def prefix =
        "${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}"

    def technologyOption =
        meta.data_type == 'ONT' ? '--ont' : ''

    // Convert bp to hifiasm genome-size format
    def hgSize = meta.hg_size >= 1000000 ?
        "${meta.hg_size / 1000000.0}m" :
        "${meta.hg_size / 1000.0}k"

    """
    hifiasm \
        ${technologyOption} \
        -t ${task.cpus} \
        --hg-size ${hgSize} \
	${params.hifiasm_args} \
        -o ${prefix} \
        ${fastq}

    test -s ${prefix}.bp.hap1.p_ctg.gfa
    test -s ${prefix}.bp.hap2.p_ctg.gfa
    """
}
