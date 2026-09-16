#!/bin/bash

#SBATCH --partition=cluster_long
#SBATCH --job-name=03b_FastQC
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

reads_trimmed_UMIextracted=(
SNRNP70-Con-rep1_L1_R1.trimmed.umi.fq.gz
SNRNP70-Con-rep2_L2_R1.trimmed.umi.fq.gz
mVenus-Con-rep5_L1_R1.trimmed.umi.fq.gz
mVenus-Con-rep6_L2_R1.trimmed.umi.fq.gz
)

FASTQ=output/${reads_trimmed_UMIextracted[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$(basename ${FASTQ} .fq.gz)

zcat ${FASTQ} | fastqc stdin:${SAMPLE} --extract --outdir output -t 1
