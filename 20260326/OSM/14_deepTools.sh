#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=14_deepTools
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
conda activate deeptools
conda list

mkdir -p output

mappings_sorted_deduplicated=(
SNRNP70-Sorbitol-rep1_L1_R1.genome.Aligned.sort.dedup.bam
SNRNP70-Sorbitol-rep2_L2_R1.genome.Aligned.sort.dedup.bam
mVenus-Sorbitol-rep5_L1_R1.genome.Aligned.sort.dedup.bam
mVenus-Sorbitol-rep6_L2_R1.genome.Aligned.sort.dedup.bam
)

BAM=output/${mappings_sorted_deduplicated[$SLURM_ARRAY_TASK_ID]}
sample=$(basename "${BAM}" .genome.Aligned.sort.dedup.bam)

bamCoverage --normalizeUsing CPM -bs 1 -p 52 -v -b "${BAM}" -o "output/${sample}.cpm.bw"
