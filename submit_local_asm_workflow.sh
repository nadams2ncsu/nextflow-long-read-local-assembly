#!/usr/bin/env bash

set -euo pipefail

##########################################################################################
# Submit Nextflow long-read local assembly workflow to SLURM
#
# Usage:
#
# ./submit_local_asm_workflow.sh \
#     <samples.tsv> \
#     <reference.fasta> \
#     [options]
#
# Required positional arguments:
#
#   samples.tsv
#       Sample configuration table. Locus-specific QC and haplotype
#       configuration files are specified within this file.
#
#   reference.fasta
#       Reference genome used for downstream assembly alignment.
#
# Optional arguments:
#
#   --hifiasm_args "ARGS"
#   --minimap2_args "ARGS"
#   --large_indel_threshold INT
#   --min_region_coverage FLOAT
#   --min_median_read_length INT
#   --group_similarity_threshold FLOAT
#   --queue_size INT
#   --outdir DIR
#


############################################
# I/O

mkdir -p logs


############################################
# Check required arguments

if [[ $# -lt 2 ]]; then

    echo "Usage:"
    echo "  $0 <samples.tsv> <reference.fasta> [options]"
    echo
    echo "Options:"
    echo '  --hifiasm_args "ARGS"'
    echo '  --minimap2_args "ARGS"'
    echo "  --large_indel_threshold INT"
    echo "  --min_region_coverage FLOAT"
    echo "  --min_median_read_length INT"
    echo "  --group_similarity_threshold FLOAT"
    echo "  --queue_size INT"
    echo "  --outdir DIR"
    echo

    exit 1
fi


############################################
# Required positional arguments

SAMPLES="$(realpath "$1")"
REFERENCE="$(realpath "$2")"

shift 2


############################################
# Default parameters

# Hifiasm
HIFIASM_ARGS=""

# Minimap2
MINIMAP2_ARGS=""

# Alignment analysis
LARGE_INDEL_THRESHOLD="40000"
MIN_REGION_COVERAGE="0.80"

# Read statistics QC
MIN_MEDIAN_READ_LENGTH="5000"

# Across-group comparison QC
GROUP_SIMILARITY_THRESHOLD="90.0"

# Maximum number of Nextflow tasks submitted to SLURM
QUEUE_SIZE="18"

# Output
OUTDIR="results"


############################################
# Parse optional workflow arguments

while [[ $# -gt 0 ]]; do
    case "$1" in
        --hifiasm_args)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --hifiasm_args requires a value."
                exit 1
            }

            HIFIASM_ARGS="$2"
            shift 2
            ;;


        --minimap2_args)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --minimap2_args requires a value."
                exit 1
            }

            MINIMAP2_ARGS="$2"
            shift 2
            ;;


        --large_indel_threshold)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --large_indel_threshold requires a value."
                exit 1
            }
            LARGE_INDEL_THRESHOLD="$2"
            shift 2
            ;;

        --min_region_coverage)

            [[ $# -ge 2 ]] || {
                echo "ERROR: --min_region_coverage requires a value."
                exit 1
            }

            MIN_REGION_COVERAGE="$2"
            shift 2
            ;;


        --min_median_read_length)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --min_median_read_length requires a value."
                exit 1
            }

            MIN_MEDIAN_READ_LENGTH="$2"
            shift 2
            ;;


        --group_similarity_threshold)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --group_similarity_threshold requires a value."
                exit 1
            }

            GROUP_SIMILARITY_THRESHOLD="$2"
            shift 2
            ;;


        --queue_size)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --queue_size requires a value."
                exit 1
            }

            QUEUE_SIZE="$2"
            shift 2
            ;;


        --outdir)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --outdir requires a value."
                exit 1
            }
            OUTDIR="$2"
            shift 2
            ;;


        *)
            echo "ERROR: Unknown argument: $1"
            exit 1
            ;;
    esac
done


############################################
# Validate input files

[[ -f "$SAMPLES" ]] || {
    echo "ERROR: samples file not found: $SAMPLES"
    exit 1
}

[[ -f "$REFERENCE" ]] || {
    echo "ERROR: reference file not found: $REFERENCE"
    exit 1
}


############################################
# Validate numeric parameters

if ! [[ "$LARGE_INDEL_THRESHOLD" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --large_indel_threshold must be an integer."
    exit 1

fi


if ! [[ "$MIN_REGION_COVERAGE" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]]; then
    echo "ERROR: --min_region_coverage must be between 0 and 1."
    exit 1

fi


if ! [[ "$MIN_MEDIAN_READ_LENGTH" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --min_median_read_length must be an integer."
    exit 1

fi


if ! [[ "$GROUP_SIMILARITY_THRESHOLD" =~ ^([0-9]+([.][0-9]+)?)$ ]]; then
    echo "ERROR: --group_similarity_threshold must be numeric."
    exit 1

fi


if ! [[ "$QUEUE_SIZE" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --queue_size must be a positive integer."
    exit 1
fi


############################################
# Print configuration

echo
echo "Submitting local assembly workflow"
echo "---------------------------------"
echo "Samples:                    $SAMPLES"
echo "Reference:                  $REFERENCE"

echo

if [[ -n "$HIFIASM_ARGS" ]]; then
    echo "Hifiasm arguments:          $HIFIASM_ARGS"

else
    echo "Hifiasm arguments:          workflow default (-t --hg-size -o)"
fi


if [[ -n "$MINIMAP2_ARGS" ]]; then
    echo "Minimap2 arguments:         $MINIMAP2_ARGS"
else
    echo "Minimap2 arguments:         workflow defaults (-ax asm5 --eqx --secondary=no)"
fi


echo "Large INDEL threshold:      $LARGE_INDEL_THRESHOLD"
echo "Min region coverage:        $MIN_REGION_COVERAGE"
echo "Min median read length:     $MIN_MEDIAN_READ_LENGTH"
echo "Group similarity threshold: $GROUP_SIMILARITY_THRESHOLD"
echo "Nextflow queue size:        $QUEUE_SIZE"
echo "Output directory:           $OUTDIR"

echo


############################################
# Submit Nextflow controller to SLURM

SBATCH_ARGS=(

    --job-name=local_asm_nf
    --cpus-per-task=1
    --mem=4G
    --time=12:00:00
    --output=logs/nextflow_%j.out
    --error=logs/nextflow_%j.err

    --export=ALL,SAMPLES="$SAMPLES",REFERENCE="$REFERENCE",HIFIASM_ARGS="$HIFIASM_ARGS",MINIMAP2_ARGS="$MINIMAP2_ARGS",LARGE_INDEL_THRESHOLD="$LARGE_INDEL_THRESHOLD",MIN_REGION_COVERAGE="$MIN_REGION_COVERAGE",MIN_MEDIAN_READ_LENGTH="$MIN_MEDIAN_READ_LENGTH",GROUP_SIMILARITY_THRESHOLD="$GROUP_SIMILARITY_THRESHOLD",QUEUE_SIZE="$QUEUE_SIZE",OUTDIR="$OUTDIR"

)


sbatch "${SBATCH_ARGS[@]}" run_nextflow.sh
