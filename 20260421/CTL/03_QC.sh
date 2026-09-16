#!/bin/bash

set -euo pipefail

scripts=(
03a_FastQC.sh
03b_FastQC.sh
)

for script in "${scripts[@]}"; do
    [ -f "$script" ] || { echo "$script is missing; please check carefully..."; exit 1; }
done

prev_job=""

for script in "${scripts[@]}"; do
    echo "Submitting Slurm job for $script"

    if [ -z "$prev_job" ]; then
        job_id=$(sbatch "$script" | awk '{print $4}')
    else
        job_id=$(sbatch --dependency=afterok:$prev_job "$script" | awk '{print $4}')
    fi

    echo " -> Slurm job ID: $job_id"
    prev_job=$job_id
done

