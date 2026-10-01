/*
################################################################################
Module: READ_STATS
################################################################################

Calculate regional read-length statistics from the original aligned BAM.
        - This process runs once per sample 

Default QC threshold: median read length <= 5000 bp ==> WARNING
        - The threshold can be overridden with --min_median_read_length

Input: sample configuration file

Output: per-sample gene read statistics tsv

################################################################################
*/

process READ_STATS {

    tag "${meta.group} | ${meta.sample} | ${meta.gene}"

    cpus 1
    memory '8 GB'

    input:

    tuple val(meta),
          path(bam),
          path(bai)

    output:

    tuple val(meta),
          path("${meta.group}_${meta.sample}_${meta.gene}.read_stats.tsv"),
          emit: stats

    script:

    def prefix =
        "${meta.group}_${meta.sample}_${meta.gene}"

    """
    ${projectDir}/bin/read_stats.py \
        --bam ${bam} \
        --coordinates "${meta.coordinates}" \
        --sample "${meta.sample}" \
        --group "${meta.group}" \
        --gene "${meta.gene}" \
        --min-median-read-length ${params.min_median_read_length} \
        --output ${prefix}.read_stats.tsv

    test -s ${prefix}.read_stats.tsv
    """
}
