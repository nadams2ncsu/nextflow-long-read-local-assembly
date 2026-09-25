#!/usr/bin/env bash

set -euo pipefail

###################################################################################################################
# Suggest a conservative Nextflow queue size for SLURM.
#   - Not all systems are the same, use suggest_queue_size.sh to determine what is best for
#        your system.
#   - Update queue size with job submission:
#        ./submit_local_assembly.sh config/samples.tsv /path/to/reference/genome/fasta --queue_size INT
#

DEFAULT=10
CAP=25

echo "#################################"
echo "# Nextflow queue-size helper"
echo "#"
echo

if ! command -v squeue >/dev/null 2>&1; then
    echo "SLURM commands were not found."
    echo "Suggested queue size: $DEFAULT"
    exit 0
fi

suggested=$DEFAULT
current_jobs=$(squeue -h -u "$USER" 2>/dev/null | wc -l | tr -d ' ')

echo "User:                  $USER"
echo "Current SLURM jobs:    $current_jobs"

max_submit=""

if command -v sacctmgr >/dev/null 2>&1; then
    max_submit=$(
        sacctmgr -n -P show assoc where user="$USER" \
            format=MaxSubmitJobs 2>/dev/null |
        awk -F'|' '
            $1 ~ /^[0-9]+$/ && $1 > 0 {
                if (min == "" || $1 < min) min=$1
            }
            END { if (min != "") print min }
        '
    )
fi

if [[ -n "$max_submit" ]]; then
    echo "Max submitted jobs:    $max_submit"
    available=$((max_submit - current_jobs))

    if (( available < 1 )); then
        suggested=1
    elif (( available < suggested )); then
        suggested=$available
    fi
else
    echo "Max submitted jobs:    not reported"
fi

if (( suggested > CAP )); then
    suggested=$CAP
fi

echo
echo "Suggested queue size:  $suggested"
echo
echo "Use:"
echo "  --queue_size $suggested"
echo
echo "This is a conservative starting value."
echo "SLURM still determines how many submitted tasks can actually run."
