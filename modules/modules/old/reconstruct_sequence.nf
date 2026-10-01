/*
###################################################################################
Module: RECONSTRUCT_SEQUENCE
###################################################################################

Reconstruct paired contiguous haplotype sequences from aligned BAM files and
trim them to the exact gene coordinates supplied in config/samples.tsv.
	- This process should only receive sample + group + gene + flank combinations
		where BOTH hap1 and hap2 were classified as contiguous.

Input: sample metadata information & aligned, sorted contig BAM files

Output: trimmed haplotype sequence per sample

###################################################################################
*/

process RECONSTRUCT_SEQUENCE {

    tag "${meta.sample} | ${meta.group} | ${meta.gene} | ${meta.run_name}"

    publishDir "${params.outdir}/fasta", mode: 'copy'

    input:

    tuple val(meta),
          path(hap1_bam),
          path(hap1_bai),
          path(hap2_bam),
          path(hap2_bai)

    output:

    tuple val(meta),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.hap1.fa"),
          path("${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}.hap2.fa"),
          emit: fasta

    script:

    def prefix =
        "${meta.group}_${meta.sample}_${meta.gene}_${meta.run_name}"

    """
    ${projectDir}/bin/reconstruct_sequence.py \
        --bam ${hap1_bam} \
        --coordinates "${meta.coordinates}" \
        --output ${prefix}.hap1.fa

    ${projectDir}/bin/reconstruct_sequence.py \
        --bam ${hap2_bam} \
        --coordinates "${meta.coordinates}" \
        --output ${prefix}.hap2.fa

    test -s ${prefix}.hap1.fa
    test -s ${prefix}.hap2.fa
    """
}
