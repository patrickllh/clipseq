#!/bin/bash

#SBATCH --partition=cluster_short
#SBATCH --job-name=11_SAMtools
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=%x_%j.log
#SBATCH --mem=0
#SBATCH --time=4:00:00

set -euo pipefail
source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
ulimit -n 65536
conda activate skipper
conda list

mkdir -p output

mappings_sorted_deduplicated=(
SNRNP70-Sorbitol-rep1_L1_R1.genome.Aligned.sort.dedup.bam
SNRNP70-Sorbitol-rep2_L2_R1.genome.Aligned.sort.dedup.bam
mVenus-Sorbitol-rep5_L1_R1.genome.Aligned.sort.dedup.bam
mVenus-Sorbitol-rep6_L2_R1.genome.Aligned.sort.dedup.bam
)

for BAM in "${mappings_sorted_deduplicated[@]}"; do
    FILE=~/clipseq/20260421/OSM/output/${BAM}
    samtools view -c -F 4 "${FILE}" > "output/${BAM}.readnum"
done
