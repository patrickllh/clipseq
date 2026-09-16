#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=06_STAR
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

STAR_DIR=output/genome_ref

reads_trimmed_UMIextracted=(
SNRNP70-Sorbitol-rep1_L1_R1.trimmed.umi.fq.gz
SNRNP70-Sorbitol-rep2_L2_R1.trimmed.umi.fq.gz
mVenus-Sorbitol-rep5_L1_R1.trimmed.umi.fq.gz
mVenus-Sorbitol-rep6_L2_R1.trimmed.umi.fq.gz
)

fq=output/${reads_trimmed_UMIextracted[$SLURM_ARRAY_TASK_ID]}
sample=$(basename ${fq} .trimmed.umi.fq.gz)

STAR \
--alignEndsType EndToEnd \
--genomeDir $STAR_DIR \
--genomeLoad NoSharedMemory \
--outBAMcompression 10 \
--outFileNamePrefix output/${sample}.genome. \
--winAnchorMultimapNmax 100 \
--outFilterMultimapNmax 100 \
--outFilterMultimapScoreRange 1 \
--outSAMmultNmax 1 \
--outMultimapperOrder Random \
--outFilterScoreMin 10 \
--outFilterType BySJout \
--limitOutSJcollapsed 5000000 \
--outReadsUnmapped None \
--outSAMattrRGline ID:${sample} \
--outSAMattributes All \
--outSAMmode Full \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--readFilesCommand zcat \
--outStd Log \
--readFilesIn $fq \
--runMode alignReads \
--runThreadN 52
