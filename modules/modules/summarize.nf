/*
###################################################################################
Module: SUMMARIZE
###################################################################################

Generate assembly/alignment summary statistics.

Input: Sample configuration file; aligned, sorted contigs to the reference genome
       BAM file and its index

Output: Summary statistics on aligned contigs

###################################################################################
*/

process SUMMARIZE {

    tag "${meta.sample} | ${meta.gene} | ${meta.run_name} | ${meta.haplotype}"

    cpus 1
    memory '4 GB'

    input:

    tuple val(meta),
          path(bam),
          path(bai),
          val(haplotype_config)


    output:

    tuple val(meta),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.${meta.haplotype}.summary.tsv"),
          emit: summary


    script:

    def prefix =
        "${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.${meta.haplotype}"

    def haplotype_arg =
        haplotype_config.toUpperCase() != 'NA'
            ? "--haplotype-config ${haplotype_config}"
            : ""

    """
    ${projectDir}/bin/summarize_alignment.py \
        --bam ${bam} \
        --coordinates "${meta.coordinates}" \
        ${haplotype_arg} \
        --gene "${meta.gene}" \
        --sample "${meta.sample}" \
        --group "${meta.group}" \
        --flank ${meta.flank_kb} \
        --assembly-haplotype "${meta.haplotype}" \
        --large-indel-threshold ${params.large_indel_threshold} \
        --min-region-coverage ${params.min_region_coverage} \
        --output ${prefix}.summary.tsv

    test -s ${prefix}.summary.tsv
    """
}
