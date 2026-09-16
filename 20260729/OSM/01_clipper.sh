#!/bin/bash

#SBATCH --partition=cluster_short
#SBATCH --job-name=01_clipper
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=%x_%A_%a.log
#SBATCH --mem=0
#SBATCH --time=4:00:00
#SBATCH --array=0-3

source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
ulimit -n 65536
conda activate clipper3
set -euo pipefail # set this after environment
mkdir -p output

mappings_sorted_deduplicated=(
SNRNP70-Sorbitol-rep1_L1_R1.genome.Aligned.sort.dedup.bam
SNRNP70-Sorbitol-rep2_L2_R1.genome.Aligned.sort.dedup.bam
mVenus-Sorbitol-rep5_L1_R1.genome.Aligned.sort.dedup.bam
mVenus-Sorbitol-rep6_L2_R1.genome.Aligned.sort.dedup.bam
)

peaks=(
SNRNP70-Sorbitol-rep1_L1_R1.peaks.bed
SNRNP70-Sorbitol-rep2_L2_R1.peaks.bed
mVenus-Sorbitol-rep5_L1_R1.peaks.bed
mVenus-Sorbitol-rep6_L2_R1.peaks.bed
)

clipper --processors 48 -v -s GRCh38v35noalt \
-b ~/clipseq/20260421/OSM/output/${mappings_sorted_deduplicated[$SLURM_ARRAY_TASK_ID]} \
-o output/${peaks[$SLURM_ARRAY_TASK_ID]}
