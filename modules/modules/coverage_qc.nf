/*
###################################################################################
Module: COVERAGE_QC
###################################################################################

Calculate whole-genome and gene-level coverage QC from the original aligned BAM.

The QC script:
    1. Calculates average whole-genome depth.
    2. Selects gene coordinates from the user-supplied YAML based on group label.
	- Can be simplified based just use one set of fix coordinates
    3. Calculates mean depth for each configured gene interval.
    4. Sums interval depths when multiple intervals are defined for a gene.
    5. Normalizes gene depth by average whole-genome depth.

Input:
    tuple(meta, BAM)
    QC configuration YAML

Output:
    tuple(meta, coverage QC TSV)

###################################################################################
*/

process COVERAGE_QC {

    tag "${meta.group} | ${meta.sample} | ${meta.gene}"

    cpus 2
    memory '8 GB'
	
    input:

    tuple val(meta),
          path(bam),
	  path(bai),
	  path(qc_config)

    output:

    tuple val(meta),
          path("${meta.group}_${meta.sample}_${meta.gene}.coverage_qc.tsv"),
          emit: qc

    script:

    def prefix =
        "${meta.group}_${meta.sample}_${meta.gene}"

    """
    ${projectDir}/bin/coverage_qc.py \\
        --bam ${bam} \\
        --config ${qc_config} \\
        --sample "${meta.sample}" \\
        --group "${meta.group}" \\
        --gene "${meta.gene}" \\
        --output ${prefix}.coverage_qc.tsv

    test -s ${prefix}.coverage_qc.tsv
    """
}
