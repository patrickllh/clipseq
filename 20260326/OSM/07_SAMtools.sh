#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=07_SAMtools
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=52
#SBATCH --output=%x_%A_%a.log
#SBATCH --mem=0
#SBATCH --time=100:00:00
#SBATCH --array=0-3

set -euo pipefail
source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
ulimit -n 65536
conda activate skipper
conda list

mkdir -p output

mappings=(
SNRNP70-Sorbitol-rep1_L1_R1.genome.Aligned.out.bam
SNRNP70-Sorbitol-rep2_L2_R1.genome.Aligned.out.bam
mVenus-Sorbitol-rep5_L1_R1.genome.Aligned.out.bam
mVenus-Sorbitol-rep6_L2_R1.genome.Aligned.out.bam
)

BAM=output/${mappings[$SLURM_ARRAY_TASK_ID]}
sample=$(basename ${BAM} .genome.Aligned.out.bam)

samtools sort -T output/${sample} -@ 52 -o output/${sample}.genome.Aligned.sort.bam ${BAM}
samtools index -@ 52 output/${sample}.genome.Aligned.sort.bam
