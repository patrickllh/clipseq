#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=08_UMICollapse
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

mappings_sorted=(
SNRNP70-Con-rep1_L1_R1.genome.Aligned.sort.bam
SNRNP70-Con-rep2_L2_R1.genome.Aligned.sort.bam
mVenus-Con-rep5_L1_R1.genome.Aligned.sort.bam
mVenus-Con-rep6_L2_R1.genome.Aligned.sort.bam
)

BAM=output/${mappings_sorted[$SLURM_ARRAY_TASK_ID]}
sample=$(basename ${BAM} .genome.Aligned.sort.bam)

java -server -Xms8G -Xmx8G -Xss20M \
    -jar ~/skipper/installation/UMICollapse-1.0.0/umicollapse.jar bam \
    -i ${BAM} \
    -o output/${sample}.genome.Aligned.sort.dedup.bam \
    --umi-sep : \
    --two-pass
