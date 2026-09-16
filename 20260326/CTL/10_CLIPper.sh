#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=10_CLIPper
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --output=%x_%A_%a.log
#SBATCH --mem=0
#SBATCH --time=100:00:00
#SBATCH --array=0-3

source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
ulimit -n 65536
conda activate clipper3
set -euo pipefail # set this after environment
conda list

mkdir -p output

mappings_sorted_deduplicated=(
SNRNP70-Con-rep1_L1_R1.genome.Aligned.sort.dedup.bam
SNRNP70-Con-rep2_L2_R1.genome.Aligned.sort.dedup.bam
mVenus-Con-rep5_L1_R1.genome.Aligned.sort.dedup.bam
mVenus-Con-rep6_L2_R1.genome.Aligned.sort.dedup.bam
)

peaks=(
SNRNP70-Con-rep1_L1_R1.peaks.bed
SNRNP70-Con-rep2_L2_R1.peaks.bed
mVenus-Con-rep5_L1_R1.peaks.bed
mVenus-Con-rep6_L2_R1.peaks.bed
)

clipper --processors 52 -v -s GRCh38v35noalt \
-b output/${mappings_sorted_deduplicated[$SLURM_ARRAY_TASK_ID]} \
-o output/${peaks[$SLURM_ARRAY_TASK_ID]}
