#!/bin/bash

set -euo pipefail

scripts=(
01_FastQC.sh
02_Cutadapt.sh
03_fastp.sh
04_FastQC.sh
# 05_STAR.sh
06_STAR.sh
07_SAMtools.sh
08_UMICollapse.sh
09_SAMtools.sh
10_CLIPper.sh
11_SAMtools.sh
12_merge_peaks.sh
13_BEDtools.sh
14_deepTools.sh
15_deepTools.sh
16_UCSC.sh
17_BEDtools.sh
18_UCSC.sh
# 19_deepTools.sh
20_deepTools.sh
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

