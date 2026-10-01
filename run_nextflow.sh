#!/usr/bin/env bash

set -euo pipefail

##########################################################################################
# Run Nextflow long-read local assembly workflow on an HPC compute node

############################################
# Tools

module load nextflow
module load apptainer

# Compute nodes do not have internet access
export NXF_DISABLE_CHECK_LATEST=true


############################################
# Confirm required environment variables

: "${SAMPLES:?ERROR: SAMPLES environment variable is not set.}"
: "${REFERENCE:?ERROR: REFERENCE environment variable is not set.}"


############################################
# Run information

echo
echo "Nextflow controller"
echo "-------------------"

echo "SLURM job:                  ${SLURM_JOB_ID:-NA}"
echo "Running on:                 $(hostname)"
echo
echo "Samples:                    $SAMPLES"
echo "Reference:                  $REFERENCE"

echo

if [[ -n "${HIFIASM_ARGS:-}" ]]; then
    echo "Hifiasm arguments:          $HIFIASM_ARGS"
else
    echo "Hifiasm arguments:          workflow defaults (-t --hg-size -o)"
fi


if [[ -n "${MINIMAP2_ARGS:-}" ]]; then
    echo "Minimap2 arguments:         $MINIMAP2_ARGS"
else
    echo "Minimap2 arguments:         workflow defaults (-ax asm5 --eqx --secondary=no)"

fi


echo "Large INDEL threshold:      ${LARGE_INDEL_THRESHOLD:-40000}"
echo "Min region coverage:        ${MIN_REGION_COVERAGE:-0.80}"
echo "Min median read length:     ${MIN_MEDIAN_READ_LENGTH:-5000}"
echo "Group similarity threshold: ${GROUP_SIMILARITY_THRESHOLD:-95.0}"
echo "Nextflow queue size:        ${QUEUE_SIZE:-18}"
echo "Output directory:           ${OUTDIR:-results}"

echo


############################################
# Build Nextflow command

CMD=(

    nextflow run main.nf
    -c nextflow.config
    --samples "$SAMPLES"
    --reference "$REFERENCE"
    --large_indel_threshold
    "${LARGE_INDEL_THRESHOLD:-40000}"
    --min_region_coverage
    "${MIN_REGION_COVERAGE:-0.80}"
    --min_median_read_length
    "${MIN_MEDIAN_READ_LENGTH:-5000}"
    --group_similarity_threshold
    "${GROUP_SIMILARITY_THRESHOLD:-95.0}"
    --queue_size
    "${QUEUE_SIZE:-18}"
    --outdir
    "${OUTDIR:-results}"
)


############################################
# Optional hifiasm arguments

if [[ -n "${HIFIASM_ARGS:-}" ]]; then

    CMD+=(
        --hifiasm_args
        "$HIFIASM_ARGS"
    )

fi


############################################
# Optional minimap2 arguments

if [[ -n "${MINIMAP2_ARGS:-}" ]]; then

    CMD+=(
        --minimap2_args
        "$MINIMAP2_ARGS"
    )

fi


############################################
# Print command

echo "Running:"

printf ' %q' "${CMD[@]}"

echo
echo


############################################
# Run workflow

"${CMD[@]}"
