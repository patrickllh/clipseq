#!/bin/bash

#SBATCH --partition=cluster_short
#SBATCH --job-name=04_deepTools
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=%x_%A_%a.log
#SBATCH --mem=0
#SBATCH --time=4:00:00
#SBATCH --array=0-3

set -euo pipefail
source /work/pa-hilario/anaconda3/etc/profile.d/conda.sh
ulimit -n 65536
conda activate deeptools
conda list

mkdir -p output

mappings_sorted_deduplicated=(
SNRNP70-Con-rep1_L1_R1.genome.Aligned.sort.dedup.bam
SNRNP70-Con-rep2_L2_R1.genome.Aligned.sort.dedup.bam
mVenus-Con-rep5_L1_R1.genome.Aligned.sort.dedup.bam
mVenus-Con-rep6_L2_R1.genome.Aligned.sort.dedup.bam
)

BAM=~/clipseq/20260421/CTL/output/${mappings_sorted_deduplicated[$SLURM_ARRAY_TASK_ID]}
sample=$(basename "${BAM}" .genome.Aligned.sort.dedup.bam)

bamCoverage --normalizeUsing CPM -bs 1 -p 48 -v -b "${BAM}" -o "output/${sample}.cpm.bw"
