/*
###################################################################################
Module: GFA2FASTA
###################################################################################

Convert phased hifiasm GFA assembly graphs into FASTA.

Input: Hifiasm assembly graphs for each haplotype (.gfa)

Output: Haplotype sequences in FASTA file format

###################################################################################
*/

process GFA2FASTA {

    tag "${meta.sample} | ${meta.gene} | ${meta.run_name}"

    input:
    tuple val(meta), path(hap1_gfa), path(hap2_gfa)

    output:
    tuple val(meta),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.hap1.fa"),
          emit: hap1

    tuple val(meta),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.hap2.fa"),
          emit: hap2

    script:
    def prefix = "${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}"

    """
    gfa2fasta --input ${hap1_gfa} --output ${prefix}.hap1.fa
    gfa2fasta --input ${hap2_gfa} --output ${prefix}.hap2.fa

    test -s ${prefix}.hap1.fa
    test -s ${prefix}.hap2.fa
    """
}
