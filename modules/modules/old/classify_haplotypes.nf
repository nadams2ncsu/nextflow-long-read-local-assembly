/*
###################################################################################
Module:  CLASSIFY_HAPLOTYPE
###################################################################################

Assign a first-pass structural haplotype using structural events identified
by summarize_alignment.py and haplotype definitions from a user-supplied YAML.

Input: sample metadata (config/samples.tsv) && known haplotype coordinates against
	input reference genome (config/haplotypes/human_FCGR2_3.yaml)

Output:
    tuple(meta, classified summary.tsv)
###################################################################################
*/

process CLASSIFY_HAPLOTYPE {

    tag "${meta.sample} | ${meta.gene} | ${meta.run_name} | ${meta.haplotype}"

    cpus 1
    memory '4 GB'

    input:
    tuple val(meta),
          path(summary)

    path haplotype_config

    output:
    tuple val(meta),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.${meta.haplotype}.haplotype.tsv"),
          emit: classified

    script:

    def prefix =
        "${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.${meta.haplotype}"

    """
    ${projectDir}/bin/classify_haplotypes.py \
        --summary ${summary} \
        --config ${haplotype_config} \
        --output ${prefix}.haplotype.tsv

    test -s ${prefix}.haplotype.tsv
    """
}
