/*
################################################################################
Module: VALIDATE_INPUTS
################################################################################

Check all user-supplied inputs in samples.tsv file

################################################################################
*/

process VALIDATE_INPUTS {

    tag "workflow input validation"

    input:
    path samples
    path reference
    path reference_fai
    val base_dir

    output:
    path 'validation.ok', emit: done

    script:
    """
    ${projectDir}/bin/validate_inputs.py \
        --samples ${samples} \
        --reference ${reference} \
	--reference-fai ${reference_fai} \
        --base-dir '${base_dir}' \
        --output validation.ok
    """
}
