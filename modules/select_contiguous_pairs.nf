/*
###################################################################################
Module: SELECT_CONTIGUOUS_PAIRS
###################################################################################

Identify sample + group + gene + flank combinations where BOTH assembly
haplotypes were classified as contiguous.

Input: Per-haplotype summary TSV files

Output: Manifest of paired-contiguous assemblies

###################################################################################
*/

process SELECT_CONTIGUOUS_PAIRS {

    tag "select paired-contiguous assemblies"

    input:

    path summaries

    output:

    path "paired_contiguous.tsv",
        emit: pairs

    script:

    """
    ${projectDir}/bin/select_contiguous_pairs.py \
        --input ${summaries} \
        --output paired_contiguous.tsv

    test -f paired_contiguous.tsv
    """
}
